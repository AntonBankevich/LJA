#include "supregraph.hpp"
#include "assembly_graph/ag_algorithms.hpp"
#include "common/omp_utils.hpp"

using namespace ag;

namespace spg {
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
            } else if ((e.isVirtualSuffix() && !e.isSuffix() && e.getStart().outDeg() != 1) ||
                       (e.isVirtualPrefix() && !e.isPrefix() && e.getFinish().inDeg() != 1)) {
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
        }
        return ok;
    }


    Vertex &Supregraph::mergePathSPG(const GraphPath &path, Vertex::id_type id) {
        VertexId start = path.getStart().getId();
        VertexId finish = path.getFinish().getId();
        Vertex &res = mergePath(path, id);
//        res itself must survive: it is the start/finish when the path is all suffix/all prefix edges.
        auto can_flip = [&res](Vertex &v) {
            return v != res && v != res.rc() && !v.marked() && v.inDeg() == 1 && v.incFront().isSuffix() &&
                   v.outDeg() != 1;
        };
        if (can_flip(*finish))
            flipVertexLockFree(finish->incFront());
        if (can_flip(start->rc()))
            flipVertexLockFree(start->rc().incFront());
        return res;
    }

    Vertex &MergePathOrLoop(const GraphPath &path, Supregraph &graph, Vertex::id_type id) {
        if(path.getStart().isJunction()) {
//            The start stays in the merge only if it has out-degree 1, the finish only if it has in-degree 1: other
//            vertices keep their connections. With a single vertex left there is nothing to merge.
            const bool keep_start = path.getStart().outDeg() == 1;
            const bool keep_finish = path.getFinish().inDeg() == 1;
            PathPosition from = path.firstPosition();
            PathPosition to = path.lastPosition();
            if (!keep_start)
                ++from;
            if (path.calculateSize() + 1 <= 1 + !keep_start + !keep_finish)
                return from.getVertex();
            if (!keep_finish)
                --to;
            return graph.mergePath(path.subPath(from, to), id);
        } else {
            VERIFY(path.getStart() == path.getFinish());
            // TODO: Create real loop merging here
            // graph.mergeLoop(path);
            PathPosition split = path.firstPosition();
            for(PathPosition pp = path.firstPosition() + 1; pp < path.lastPosition(); ++pp) {
                if (pp.getVertex().size() < split.getVertex().size()) {
                    split = pp;
                    break;
                }
            }
            GraphPath new_path = path.subPath(split, path.lastPosition()) + path.subPath(path.firstPosition(), split);

            split = new_path.firstPosition();
            for(PathPosition pp = new_path.firstPosition() + 1; pp < new_path.lastPosition(); ++pp) {
                if (pp.getVertex() == new_path.getStart().rc()) {
                    split = pp;
                    break;
                }
            }
            if (split != new_path.firstPosition()) {
                GraphPath p1 = new_path.subPath(new_path.firstPosition(), split);
                GraphPath p2 = new_path.subPath(split, new_path.lastPosition());
                if (p1.calculateSize() > 2)
                    graph.mergePath(p1);
                if (p2.calculateSize() > 2)
                    graph.mergePath(p2);
                return new_path.getStart().frontVertex();
            } else {
                return graph.mergePath(new_path, id);
            }
        }
    }

//    The finish of a suffix edge is flipped into its start when the edge is its only incoming edge.
    static bool FlippableAlong(Edge &edge) {
        Vertex &core = edge.getFinish();
        return edge.isSuffix() && core.inDeg() == 1 && core != edge.getStart() && core != edge.getStart().rc();
    }

//    Flips the single out-neighbour of a merged vertex, and its single in-neighbour on the rc strand, into it wherever
//    FlippableAlong allows, so that the merged vertex replaces the path ends that have no other connections.
    static void FlipNeighboursInto(Supregraph &graph, Vertex &res) {
        if (res.marked())
            return;
        if (res.outDeg() == 1 && FlippableAlong(res.front()))
            graph.flipVertexLockFree(res.front());
        if (res.inDeg() == 1 && FlippableAlong(res.incFront().rc()))
            graph.flipVertexLockFree(res.incFront().rc());
    }

    Vertex &MergeUnbranching(logging::Logger &logger, Supregraph &graph, Vertex &vertex) {
        VERIFY(!vertex.isJunction());
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
        Vertex &res = MergePathOrLoop(path, graph);
        FlipNeighboursInto(graph, res);
        return res;
    }

//    Returns the vertex each path was merged into.
    static std::vector<VertexId> MergePathsSPG(logging::Logger &logger, size_t threads, Supregraph &graph, const std::vector<GraphPath> &paths) {
        logger.trace() << "Merging unbranching paths" << std::endl;
        //TODO: improve this system. make new vertex creation take not just id but an object that knows about the
        //reserved package of ids. Inside this package there is an indexed structure that distributes ids.
        //When this container is destroyed, contents spill into the general list in the assembly graph. Spilling should
        //also be done carefully to avoid changing the order. Also something should be done about edge ids.
        int first_id = graph.reserveVertexIds(paths.size());
        ParallelRecordCollector<VertexId> merged(threads);
        std::function<void(size_t, const GraphPath &)> task =
                [&graph, &merged, first_id](size_t pos, const GraphPath &path) {
                    merged.emplace_back(MergePathOrLoop(path, graph, first_id + pos).getId());
                };
        ParallelProcessor<const GraphPath>(task, logger, threads).processObjects(paths.begin(), paths.end());
        return merged.collectUnique();
    }


    void MergeAllSPG(logging::Logger &logger, size_t threads, Supregraph &graph) {
        // graph.resetMarkers();
        omp_set_num_threads(threads);
        auto unbranching_paths = AllUnbranchingPaths(logger, threads, graph);
        std::vector<VertexId> merged = MergePathsSPG(logger, threads, graph, unbranching_paths);
        logger.trace() << "Flipping vertices" << std::endl;
        for (VertexId vid : merged)
            FlipNeighboursInto(graph, *vid);
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