#pragma once
#include "assembly_graph.hpp"
#include "common/omp_utils.hpp"
namespace ag {
    template<class Traits>
    std::vector<CompactPath<Traits>> ConstructUnbranchingPaths(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        ParallelRecordCollector<CompactPath<Traits>> result;
        std::function<void(size_t, typename Traits::Vertex &)> task =
                [&result](size_t pos, typename Traits::Edge &start) {
                    if (!start.getStart().isJunction())
                        return;
                    GraphPath<Traits> to_merge = GraphPath<Traits>::WalkForward(start);
                    if(to_merge.start() < to_merge.finish().rc() || (to_merge.start() == to_merge.finish().rc()) &&
                                                                    to_merge.frontEdge().truncSeq() < to_merge.backEdge().rc().truncSeq()) {
                        result.emplace_back(to_merge);
                    }
                };
        processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, task);
        return result.collect();
    }

    template<class Traits>
    void processLoop(ParallelRecordCollector<GraphPath<Traits>> &result, typename Traits::Vertex &start) {
        GraphPath<Traits> to_merge = GraphPath<Traits>::WalkForward(start.front());
        GraphPath<Traits> second_part;
        VERIFY(to_merge.finish() == start || to_merge.finish() == start.rc());
        if(to_merge.finish() != start && to_merge.finish() == start.rc()) {
            second_part = GraphPath<Traits>::WalkForward(start.rc().front());
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
            if(to_merge.size() > 1)
                result.emplace_back(to_merge);
            if (!second_part.empty() && second_part.size() > 1)
                result.emplace_back(second_part);
        }
    }


    //    This method can only be invoked after all linear paths are collapsed.
    template<class Traits>
    std::vector<CompactPath<Traits>> ConstructUnbranchingLoops(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        ParallelRecordCollector<CompactPath<Traits>> result;
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
    std::vector<GraphPath<Traits>> AllUnbranchingPaths(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        logger.trace() << "Collecting linear unbranching paths" << std::endl;
        ParallelRecordCollector<GraphPath<Traits>> result(threads);
        std::function<void(size_t, typename Traits::Edge &)> pathTask =
                [&result](size_t pos, typename Traits::Edge &start) {
                    if (!start.getStart().isJunction())
                        return;
                    GraphPath<Traits> to_merge = GraphPath<Traits>::WalkForward(start);
                    if(to_merge.size() == 1)
                        return;
                    for(size_t i = 1; i < to_merge.size(); i++) {
                        to_merge.getVertex(i).mark();
                    }
                    if(to_merge.start() < to_merge.finish().rc() || (to_merge.start() == to_merge.finish().rc()) &&
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
        auto linear_paths = AllUnbranchingPaths(logger, threads, graph);
        MergePaths(logger, threads, graph, linear_paths);
        logger.trace() << "Removing isolated vertices" << std::endl;
        graph.removeMarked();
        graph.removeIsolated();
        logger.trace() << "Finished merging unbranching paths" << std::endl;
    }
}
