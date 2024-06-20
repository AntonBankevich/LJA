#pragma once
#include "assembly_graph.hpp"
#include "common/omp_utils.hpp"
namespace ag {
    template<class Traits>
    class PathHelper {
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Vertex::VertexId VertexId;
        static GraphPath<Traits> WalkForward(Edge &start) {
            GraphPath<Traits> res(start);
            Vertex *next = &start.getFinish();
            VERIFY(next != nullptr);
            while (*next != start.getStart() && *next != start.getStart().rc() && !next->isJunction()) {
                VERIFY(next != nullptr);
                res += next->front();
                next = &res.getFinish();
            }
            return std::move(res);
        }
    };

    template<class Traits>
    void processLoop(ParallelRecordCollector<GraphPath<Traits>> &result, typename Traits::Vertex &start) {
        GraphPath<Traits> to_merge = PathHelper<Traits>::WalkForward(start.front());
        GraphPath<Traits> second_part;
        VERIFY(to_merge.getFinish() == start || to_merge.getFinish() == start.rc());
        if(to_merge.getFinish() != start && to_merge.getFinish() == start.rc()) {
            second_part = PathHelper<Traits>::WalkForward(start.rc().front());
        }
        bool ok = true;
        for(typename Traits::Vertex &v : to_merge.vertices()) {
            if(v < start || v.rc() < start) {
                ok = false;
                break;
            }
        }
        for(typename Traits::Vertex &v : second_part.vertices()) {
            if(v < start || v.rc() < start) {
                ok = false;
                break;
            }
        }
        if(ok) {
            if(!to_merge.isSingleton())
                result.emplace_back(to_merge);
            if (!second_part.empty() && !second_part.isSingleton())
                result.emplace_back(second_part);
        }
    }


    //    This method can only be invoked after all linear paths are collapsed.
    template<class Traits>
    std::vector<GraphPath<Traits>> ConstructUnbranchingLoops(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        ParallelRecordCollector<GraphPath<Traits>> result;
        std::function<void(size_t, typename Traits::Vertex &)> task =
                [&result](size_t pos, typename Traits::Vertex &start) {
                    if (start.isJunction())
                        return;
                    processLoop(result, start);
                };
        processObjects(graph.vertices().begin(), graph.vertices().end(), logger, threads, task);
        return result.collect();
    }


    template<class Traits>
    std::vector<GraphPath<Traits>> AllUnbranchingPaths(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph, bool include_trivial = true) {
        logger.trace() << "Collecting linear unbranching paths" << std::endl;
        ParallelRecordCollector<GraphPath<Traits>> result(threads);
        std::function<void(size_t, typename Traits::Edge &)> pathTask =
                [&result, include_trivial](size_t pos, typename Traits::Edge &start) {
                    if (!start.getStart().isJunction())
                        return;
                    GraphPath<Traits> to_merge = PathHelper<Traits>::WalkForward(start);
                    if(to_merge.isSingleton())
                        return;
                    for(auto &v : to_merge.innerVertices()) {
                        v.mark();
                    }
                    if(to_merge.getStart() < to_merge.getFinish().rc() || (to_merge.getStart() == to_merge.getFinish().rc()) &&
                                (include_trivial || !to_merge.isSingleton()) &&
                                to_merge.frontEdge().truncSeq() <= to_merge.backEdge().rc().truncSeq()) {
                        result.emplace_back(to_merge);
                    }
                };
        processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, pathTask);
        logger.trace() << "Collecting circular unbranching paths" << std::endl;
        std::function<void(size_t, typename Traits::Vertex &)> loopTask =
                [&result](size_t pos, typename Traits::Vertex &start) {
                    if (start.isJunction() || start.marked()) {
                        start.unmark();
                        return;
                    }
                    processLoop(result, start);
                };
        processObjects(graph.vertices().begin(), graph.vertices().end(), logger, threads, loopTask);
        return result.collect();
    }


    template<class Traits>
    void MergePaths(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph, const std::vector<GraphPath<Traits>> &paths) {
        logger.trace() << "Merging unbranching paths" << std::endl;
        std::function<void(size_t, const GraphPath<Traits> &)> task =
                                                                       [&graph](size_t pos, const GraphPath<Traits> &path) {
                                                                           graph.mergePathToEdge(path);
                                                                       };
        ParallelProcessor<const GraphPath<Traits>>(task, logger, threads).processObjects(paths.begin(), paths.end());
    }

    template<class Traits>
    void MergeAll(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        graph.resetMarkers();
        auto linear_paths = AllUnbranchingPaths(logger, threads, graph, false);
        MergePaths(logger, threads, graph, linear_paths);
        logger.trace() << "Removing isolated vertices" << std::endl;
        graph.removeMarked();
        graph.removeIsolated();
        logger.trace() << "Finished merging unbranching paths" << std::endl;
    }
}
