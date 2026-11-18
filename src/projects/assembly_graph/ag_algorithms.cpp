#include "ag_algorithms.hpp"
#include "common/graph_coloring.hpp"
namespace ag {

    GraphPath PathHelper::WalkForward(Edge &start) {
        GraphPath res(start);
        Vertex *next = &start.getFinish();
        VERIFY(next != nullptr);
        while (*next != start.getStart() && *next != start.getStart().rc() && !next->isJunction()) {
            VERIFY(next != nullptr);
            res += next->front();
            next = &res.getFinish();
        }
        return res;
    }

    void processLoop(ParallelRecordCollector<GraphPath> &result, Vertex &start) {
        GraphPath to_merge = PathHelper::WalkForward(start.front());
        GraphPath second_part;
        VERIFY(to_merge.getFinish() == start || to_merge.getFinish() == start.rc());
        if (to_merge.getFinish() != start && to_merge.getFinish() == start.rc()) {
            second_part = PathHelper::WalkForward(start.rc().front());
        }
        bool ok = true;
        for (Vertex &v: to_merge.vertices()) {
            if (v < start || v.rc() < start) {
                ok = false;
                break;
            }
        }
        for (Vertex &v: second_part.vertices()) {
            if (v < start || v.rc() < start) {
                ok = false;
                break;
            }
        }
        if (ok) {
            if (!to_merge.isSingleton())
                result.emplace_back(to_merge);
            if (!second_part.empty() && !second_part.isSingleton())
                result.emplace_back(second_part);
        }
    }

    std::vector<GraphPath> ConstructUnbranchingLoops(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
        ParallelRecordCollector<GraphPath> result(threads);
        std::function<void(size_t, Vertex &)> task =
                [&result](size_t pos, Vertex &start) {
                    if (start.isJunction())
                        return;
                    processLoop(result, start);
                };
        processObjects(graph.vertices().begin(), graph.vertices().end(), logger, threads, task);
        return result.collect();
    }

    std::vector<GraphPath>
    AllUnbranchingPaths(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
        logger.trace() << "Collecting linear unbranching paths" << std::endl;
        ParallelRecordCollector<GraphPath> result(threads);
        std::function<void(size_t, Edge &)> pathTask =
                [&result](size_t pos, Edge &start) {
                    if (!start.getStart().isJunction())
                        return;
                    GraphPath to_merge = PathHelper::WalkForward(start);
                    for (auto &v: to_merge.innerVertices()) {
                        v.mark();
                    }
                    if (to_merge.frontEdge() > to_merge.backEdge().rc())
                        return;
                    if (to_merge.calculateSize() == 2 && to_merge.frontEdge().isPrefix() && to_merge.backEdge().isSuffix())
                        return;
                    result.emplace_back(to_merge);
                };
        processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, pathTask);
        logger.trace() << "Collecting circular unbranching paths" << std::endl;
        std::function<void(size_t, Vertex &)> loopTask =
                [&result](size_t pos, Vertex &start) {
                    if (start.isJunction() || start.marked()) {
                        start.unmark();
                        return;
                    }
                    processLoop(result, start);
                };
        processObjects(graph.vertices().begin(), graph.vertices().end(), logger, threads, loopTask);
        return result.collect();
    }

//    An unbranching path ends at junctions and its inner vertices belong to it alone (an isolated cycle touches no other
//    path), so paths conflict only through their ends. Pairs (canonical end, path) are sorted in parallel; paths with the
//    same end conflict. Once the paths are trimmed for merging, paths of one bucket are not even joined by an edge: an
//    edge u -> v between ends of two others would need u's out-degree and v's in-degree above 1, which the supregraph
//    invariant forbids. Buckets are the color classes of the conflict graph. With vertex degrees at most 8 a path
//    conflicts with at most 30 others.
    std::vector<std::vector<size_t>>
    IndependentPathBuckets(logging::Logger &logger, size_t threads, const std::vector<GraphPath> &paths) {
        logger.trace() << "Splitting " << paths.size() << " paths into independent buckets" << std::endl;
        omp_set_num_threads(threads);
        const size_t n = paths.size();
        auto canonical = [](const Vertex &v) {return v.isCanonical() ? v.getId() : v.rc().getId();};
        std::vector<std::pair<ConstVertexId, size_t>> ends(2 * n);
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n; i++) {
            ends[2 * i] = {canonical(paths[i].getStart()), i};
            ends[2 * i + 1] = {canonical(paths[i].getFinish()), i};
        }
        __gnu_parallel::sort(ends.begin(), ends.end());
        std::vector<std::vector<size_t>> conflicts(n);
        for (size_t from = 0, to = 0; from < ends.size(); from = to) {
            while (to < ends.size() && ends[to].first == ends[from].first)
                to++;
            for (size_t a = from; a < to; a++)
                for (size_t b = from; b < to; b++)
                    if (ends[a].second != ends[b].second)
                        conflicts[ends[a].second].push_back(ends[b].second);
        }
        std::vector<std::vector<size_t>> buckets = ParallelColorClasses(threads, conflicts);
        logger.trace() << "Split paths into " << buckets.size() << " buckets" << std::endl;
        return buckets;
    }

    void
    MergePathsToEdges(logging::Logger &logger, size_t threads, AssemblyGraph &graph, const std::vector<GraphPath> &paths) {
        logger.trace() << "Merging unbranching paths" << std::endl;
        std::function<void(size_t, const GraphPath &)> task =
                [&graph](size_t pos, const GraphPath &path) {
                    graph.mergePathToEdge(path);
                };
        ParallelProcessor<const GraphPath>(task, logger, threads).processObjects(paths.begin(), paths.end());
    }

    void MergeAllToEdges(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
        graph.resetMarkers();
        auto unbranching_paths = AllUnbranchingPaths(logger, threads, graph);
        MergePathsToEdges(logger, threads, graph, unbranching_paths);
        logger.trace() << "Removing isolated vertices" << std::endl;
        graph.removeMarked();
        graph.removeIsolated();
        logger.trace() << "Finished merging unbranching paths" << std::endl;
    }

    AssemblyGraph
    LoadSupregraphFromGFA(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path,
                          size_t k) {
        logger.info() << "Loading graph from fasta" << std::endl;
        std::vector<std::tuple<Edge::id_type, Edge::id_type>> edges;
        std::vector<std::tuple<Sequence, Vertex::id_type>> vertices(threads);
        std::ifstream is(path);
        std::string line;
        std::unordered_map<Vertex::id_type, VertexId> vmap;
        AssemblyGraph res(k);
        while (std::getline(is, line)) {
            if (line.empty() || line[0] == '#') continue;  // skip empty/comment
            auto fields = split(line, "\t");
            if (fields.empty()) continue;
            char type = fields[0][0];
            if (type == 'S' && fields.size() >= 3) {
                auto vid = Parse<Vertex::id_type>(fields[1]);
                Sequence seq(fields[2]);
                Vertex &v = res.addVertex(seq, vid);
                vmap[v.getInnerId()] = v.getId();
                vmap[v.rc().getInnerId()] = v.rc().getId();
            } else if (type == 'L' && fields.size() >= 7) {
                Vertex::id_type from =
                        fields[2][0] == '+' ? Parse<Vertex::id_type>(fields[1]) : -Parse<Vertex::id_type>(fields[1]);
                Vertex::id_type to =
                        fields[4][0] == '+' ? Parse<Vertex::id_type>(fields[3]) : -Parse<Vertex::id_type>(fields[3]);
                size_t overlap = Parse<int>(fields[5], 0, fields[5].size());

                VERIFY_MSG(startsWith(fields[6], "ID:Z:"), line);
                ag::EdgeSaveLabel eids = Parse<ag::EdgeSaveLabel>(fields[6], 5, fields[6].size());
                Vertex &vfrom = *vmap[from];
                Vertex &vto = *vmap[to];
                VERIFY(overlap <= vfrom.size() && overlap <= vto.size());
                Edge &edge = res.addEdge(vfrom, vto, vto.getSeq().Subseq(overlap),vfrom.rc().getSeq().Subseq(overlap), eids.fId, eids.rcId);
            }
        }
        is.close();
        return res;
    }
}
