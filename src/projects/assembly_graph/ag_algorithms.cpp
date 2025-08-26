#include "ag_algorithms.hpp"
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
        return std::move(res);
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
                    if (to_merge.isSingleton())
                        return;
                    for (auto &v: to_merge.innerVertices()) {
                        v.mark();
                    }
                    if (to_merge.getStart() < to_merge.getFinish().rc() ||
                        (to_merge.getStart() == to_merge.getFinish().rc()) &&
                        !to_merge.isSingleton() &&
                        (to_merge.calculateSize() > 2 || (to_merge.calculateSize() == 2 && (!to_merge.frontEdge().isPrefix() || !to_merge.backEdge().isSuffix()))) &&
                        to_merge.frontEdge().truncSeq() <= to_merge.backEdge().rc().truncSeq()) {
                        result.emplace_back(to_merge);
                    }
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

    void
    MergePathsToEdges(logging::Logger &logger, size_t threads, AssemblyGraph &graph, const std::vector<GraphPath> &paths) {
        logger.trace() << "Merging unbranching paths" << std::endl;
        std::function<void(size_t, const GraphPath &)> task =
                [&graph](size_t pos, const GraphPath &path) {
                    graph.mergePathToEdge(path);
                };
        ParallelProcessor<const GraphPath>(task, logger, threads).processObjects(paths.begin(), paths.end());
    }

    Vertex &MergePathSPG(const GraphPath &path, AssemblyGraph &graph) {
        if(path.getStart().isJunction())
            return graph.mergePath(path);
        else {
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
                return new_path.getStart().front().getFinish();
            } else {
                return graph.mergePath(new_path);
            }
        }
    }

    void
    MergePathsSPG(logging::Logger &logger, AssemblyGraph &graph, const std::vector<GraphPath> &paths) {
        logger.trace() << "Merging unbranching paths" << std::endl;
        std::function<void(size_t, const GraphPath &)> task =
                [&graph](size_t pos, const GraphPath &path) {
                    MergePathSPG(path, graph);
                };
        ParallelProcessor<const GraphPath>(task, logger, 1).processObjects(paths.begin(), paths.end());
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

    // TODO: make this work in parallel after concurrent adding vertices to the graph is implemented.
    void MergeAllSPG(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
        graph.resetMarkers();
        auto unbranching_paths = AllUnbranchingPaths(logger, threads, graph);
        MergePathsSPG(logger, graph, unbranching_paths);
        logger.trace() << "Removing isolated vertices" << std::endl;
        graph.removeMarked();
        graph.removeIsolated();
        logger.trace() << "Finished merging unbranching paths" << std::endl;
    }

    AssemblyGraph
    LoadSupregraphFromGFA(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path) {
        logger.info() << "Loading graph from fasta" << std::endl;
        std::vector<std::tuple<Edge::id_type, Edge::id_type>> edges;
        std::vector<std::tuple<Sequence, Vertex::id_type>> vertices(threads);
        std::ifstream is(path);
        std::string line;
        std::unordered_map<Vertex::id_type, VertexId> vmap;
        AssemblyGraph res;
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
                VERIFY(overlap == std::min(vfrom.size(), vto.size()));
                Edge &edge = res.addEdge(vfrom, vto, vto.getSeq().Subseq(overlap),vfrom.rc().getSeq().Subseq(overlap), eids.fId, eids.rcId);
            }
        }
        is.close();
        return std::move(res);
    }

}
