#include "supregraph.hpp"
#include "assembly_graph/ag_algorithms.hpp"
#include "common/omp_utils.hpp"

using namespace ag;

namespace spg {
    bool Supregraph::segmentsFitEdge(ag::Edge &e) const {
        Vertex &u = e.getStart();
        Vertex &v = e.getFinish();
        const size_t t = v.size() - e.truncSize();
        const bool u_last_at_end = !u.coverage_info.empty() && u.coverage_info.back().finish() == u.size();
        if (v.coverage_info.empty())
            return u_last_at_end;
//        A kpomer chunk stands for each of its (k+1)-mers, so its first segment ends k + 1 after its start.
        ag::CoverageSamples::SampleView first = v.coverage_info.front();
        const size_t v_first = first.sample.type == ag::CoverageSamples::SampleType::kpomer ?
                               first.start() + getK() + 1 : first.finish();
        return v_first == t + 1 || (u_last_at_end && v_first == t);
    }

    bool Supregraph::checkConsistency(std::ostream &os) {
        bool ok = AssemblyGraph::checkConsistency(os);
        for (ag::Edge &e : edges()) {
            if (e.getStart().outDeg() != 1 && e.getFinish().inDeg() != 1) {
                os << "Supregraph inconsistency: edge " << e.getId() << " " << e.getStart().getId() << " -> "
                   << e.getFinish().getId() << " is neither virtually suffix (out-degree " << e.getStart().outDeg()
                   << ") nor virtually prefix (in-degree " << e.getFinish().inDeg() << ")\n";
                ok = false;
            }
            if (e.getVirtualType() == ag::VirtualType::unset) {
                os << "Supregraph inconsistency: edge " << e.getId() << " has no virtual type\n";
                ok = false;
            } else if (((e.isVirtualSuffix() || e.isMerge()) && !e.isSuffix() && e.getStart().outDeg() != 1) ||
                       ((e.isVirtualPrefix() || e.isMerge()) && !e.isPrefix() && e.getFinish().inDeg() != 1)) {
                os << "Supregraph inconsistency: edge " << e.getId() << " has virtual type " << int(e.getVirtualType())
                   << " but start out-degree " << e.getStart().outDeg() << " and finish in-degree "
                   << e.getFinish().inDeg() << "\n";
                ok = false;
            }
            if (e.getStart() == e.getFinish()) {
                bool alone = true;
                for (ag::Vertex &x : ThisAndRC(e.getStart()))
                    for (ag::Edge &other : x)
                        alone = alone && other.getFinish() == x;
                if (!alone) {
                    os << "Supregraph inconsistency: edge " << e.getId() << " goes from vertex " << e.getStart().getId()
                       << " to itself, but the vertex has other edges\n";
                    ok = false;
                }
            }
            if (e.getStart() != e.getFinish() && !segmentsFitEdge(e)) {
                ag::Vertex &u = e.getStart();
                ag::Vertex &v = e.getFinish();
                os << "Supregraph inconsistency: sample segments at edge " << e.getId() << " " << u.getId() << " -> "
                   << v.getId() << " (overlap " << __int64_t(v.size()) - __int64_t(e.truncSize())
                   << ") do not meet at the overlap: last segment of the start (size " << u.size() << ") ";
                if (u.coverage_info.empty())
                    os << "absent";
                else
                    os << "ends at " << u.coverage_info.back().finish();
                os << ", first segment of the finish ";
                if (v.coverage_info.empty())
                    os << "absent";
                else
                    os << "starts at " << v.coverage_info.front().start() << " (" <<
                       (v.coverage_info.front().sample.type == ag::CoverageSamples::SampleType::kpomer ? "kpomer" : "resolution")
                       << ")";
                os << "\n";
                ok = false;
            }
        }
        return ok;
    }


    bool CanFlip(const ag::Edge &inc) {
        const Vertex &core = inc.getFinish();
        if (!inc.isSuffix() || core.inDeg() != 1)
            return false;
        for (const Edge &out : core)
            if (!out.isPrefix())
                return false;
        return true;
    }

    GraphPath MergeablePath(const GraphPath &path) {
        if (!path.getStart().isJunction()) {
            VERIFY(path.getStart() == path.getFinish());
            return path;
        }
        const bool keep_start = path.getStart().outDeg() == 1;
        const bool keep_finish = path.getFinish().inDeg() == 1;
        PathPosition from = path.firstPosition();
        PathPosition to = path.lastPosition();
        if (!keep_start)
            ++from;
        if (path.calculateSize() + 1 <= 1 + !keep_start + !keep_finish)
            return GraphPath(from.getVertex());
        if (!keep_finish)
            --to;
        return path.subPath(from, to);
    }

    static Vertex &MergeLockFree(const GraphPath &path, Supregraph &graph, Vertex::id_type id = 0) {
        return path.calculateSize() == 0 ? path.getStart() : graph.mergePathLockFree(path, id);
    }

//    path comes from MergeablePath. Locks nothing: MergeAllSPG merges only paths that touch disjoint parts of the graph
//    in parallel.
    static Vertex &MergeMergeable(const GraphPath &path, Supregraph &graph, Vertex::id_type id) {
        if (path.calculateSize() == 0 || path.getStart() != path.getFinish())
            return MergeLockFree(path, graph, id);
//        An isolated cycle. If it also passes through its start's rc, it is a self-rc path start -> start.rc() followed
//        by a self-rc path back: the first one is merged into a self-rc vertex, then the inner part of the second one.
//        Otherwise the cycle becomes a single vertex with a loop.
        PathPosition split = path.firstPosition();
        for (PathPosition pp = path.firstPosition() + 1; pp < path.lastPosition(); ++pp)
            if (pp.getVertex() == path.getStart().rc()) {
                split = pp;
                break;
            }
        if (split == path.firstPosition())
            return MergeLockFree(path, graph, id);
        GraphPath back = path.subPath(split, path.lastPosition());
        const bool inner_edges = back.calculateSize() >= 3;
        GraphPath back_inner = inner_edges ? back.subPath(back.firstPosition() + 1, back.lastPosition() - 1) : back;
        Vertex &res = MergeLockFree(path.subPath(path.firstPosition(), split), graph, id);
        if (inner_edges)
            MergeLockFree(back_inner, graph);
        return res;
    }

    Vertex &MergePathOrLoop(const GraphPath &path, Supregraph &graph, Vertex::id_type id) {
        return MergeMergeable(MergeablePath(path), graph, id);
    }

    Vertex &MergeUnbranching(logging::Logger &logger, Supregraph &graph, Vertex &vertex) {
//        A junction lies on no unbranching path through it: nothing to merge.
        if (vertex.isJunction())
            return vertex;
        ag::GraphPath path = ag::PathHelper::WalkForward(vertex.front());
        ag::GraphPath rpath = ag::PathHelper::WalkForward(vertex.rc().front()).RC();
        if (path.getFinish().isJunction()) {
            if (rpath.getStart().isJunction()) {
                path = rpath + path;
            } else {
                VERIFY(rpath.getStart() == rpath.getFinish().rc());
                path = path.RC() + rpath + path;
            }
        } else {
            if (path.getFinish() == path.getStart().rc()) {
                path = rpath + path;
            } else {
                VERIFY(path.getStart() == path.getFinish());
            }
        }
        logger.trace() << "Merging path if necessary: " << path.str() << std::endl;
        return MergePathOrLoop(path, graph);
    }

//    The paths are trimmed before any of them is merged: then they never share vertices and no merge removes a vertex of
//    another one. Merging one changes the edges of its neighbours, so only the paths of one independent bucket are merged
//    in parallel, without locks.
    void MergeAllSPG(logging::Logger &logger, size_t threads, Supregraph &graph) {
        omp_set_num_threads(threads);
        const std::vector<GraphPath> unbranching_paths = AllUnbranchingPaths(logger, threads, graph);
        const std::vector<std::vector<size_t>> buckets = ag::IndependentPathBuckets(logger, threads, unbranching_paths);
        std::vector<GraphPath> mergeable(unbranching_paths.size());
#pragma omp parallel for schedule(dynamic, 100)
        for (size_t i = 0; i < unbranching_paths.size(); i++)
            mergeable[i] = MergeablePath(unbranching_paths[i]);
        logger.trace() << "Merging unbranching paths in " << buckets.size() << " buckets" << std::endl;
        //TODO: improve this system. make new vertex creation take not just id but an object that knows about the
        //reserved package of ids. Inside this package there is an indexed structure that distributes ids.
        //When this container is destroyed, contents spill into the general list in the assembly graph. Spilling should
        //also be done carefully to avoid changing the order. Also something should be done about edge ids.
        int first_id = graph.reserveVertexIds(unbranching_paths.size());
        for (const std::vector<size_t> &bucket : buckets) {
#pragma omp parallel for schedule(dynamic, 10)
            for (size_t b = 0; b < bucket.size(); b++)
                MergeMergeable(mergeable[bucket[b]], graph, first_id + bucket[b]);
        }
        logger.trace() << "Removing isolated vertices" << std::endl;
        graph.removeMarked();
        logger.trace() << "Finished merging unbranching paths" << std::endl;
    }

    static bool OuterVertexHasSupport(Vertex &vertex) {
        size_t res = vertex.getRawIntSPGCoverage();
        if (vertex.inDeg() == 1 && vertex.incFrontVertex().outDeg() == 1)
            res += vertex.incFrontVertex().getRawIntSPGCoverage();
        if (vertex.outDeg() == 1 && vertex.frontVertex().inDeg() == 1)
            res += vertex.frontVertex().getRawIntSPGCoverage();
        return res > 0;
    }

    size_t SimpleRemoveUncovered(logging::Logger &logger, size_t threads, Supregraph &spg) {
        logger.trace() << "Removing completely uncovered edges" << std::endl;
        omp_set_num_threads(threads);
        std::vector<VertexId> vertices_to_delete;
        for (Vertex &v : spg.verticesUnique())
            if (v.hasCoverageInfo() && v.coverage_info.raw_support == 0)
                vertices_to_delete.push_back(v.getId());
        size_t res = vertices_to_delete.size();
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(vertices_to_delete, spg)
        for(size_t i = 0; i < vertices_to_delete.size(); i++) {
            spg.isolateAndMark(*vertices_to_delete[i]);
        }
        vertices_to_delete.clear();
        for (Vertex &v : spg.verticesUnique())
            if (!v.hasCoverageInfo() && v.outDeg() == 0 && v.inDeg() == 0)
                vertices_to_delete.push_back(v.getId());
        res += vertices_to_delete.size();
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(vertices_to_delete, spg)
        for(size_t i = 0; i < vertices_to_delete.size(); i++) {
            spg.isolateAndMark(*vertices_to_delete[i]);
        }
        spg.removeMarked();

//         std::vector<VertexId> vertices_to_check;
//         for (Vertex &vertex: spg.verticesUnique()) {
//             VERIFY(!vertex.marked());
//             if (vertex.isOuter())
//                 vertices_to_check.emplace_back(vertex.getId());
//         }
//         size_t res = 0;
//         while (!vertices_to_check.empty()) {
//             std::vector<VertexId> vertices_to_delete;
//             for(VertexId vid : vertices_to_check) {
// //                A vertex isolated in an earlier round is picked up again as a neighbour; it is already gone.
//                 if (!vid->marked() && vid->isOuter() && !OuterVertexHasSupport(*vid))
//                     vertices_to_delete.emplace_back(vid);
//             }
//             vertices_to_check.clear();
//             for(VertexId vid : vertices_to_delete) {
//                 for (Vertex & vdir : ThisAndRC(*vid))
//                     for (Edge &edge : vdir)
//                         vertices_to_check.emplace_back(edge.getFinish().getCanonical().getId());
//             }
//             std::sort(vertices_to_check.begin(), vertices_to_check.end());
//             vertices_to_check.erase(std::unique(vertices_to_check.begin(), vertices_to_check.end()), vertices_to_check.end());
// #pragma omp parallel for default(none) schedule(dynamic, 100) shared(vertices_to_delete, spg)
//             for(size_t i = 0; i < vertices_to_delete.size(); i++) {
//                 spg.isolateAndMark(*vertices_to_delete[i]);
//             }
//             logger.trace() << "Removed " << vertices_to_delete.size() << " vertices." << std::endl;
//             res += vertices_to_delete.size();
//             vertices_to_delete.clear();
//         }
        logger.trace() << "Finished removing completely uncovered vertices." << std::endl;
        MergeAllSPG(logger, threads, spg);
        return res;
    }
}