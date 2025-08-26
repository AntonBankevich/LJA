#pragma once
#include "assembly_graph.hpp"
#include "common/omp_utils.hpp"
namespace ag {
    class PathHelper {
    public:
        static GraphPath WalkForward(Edge &start);
    };
    void processLoop(ParallelRecordCollector<GraphPath> &result, Vertex &start);
    //    This method can only be invoked after all linear paths are collapsed.
    std::vector<GraphPath> ConstructUnbranchingLoops(logging::Logger &logger, size_t threads, AssemblyGraph &graph);
    std::vector<GraphPath> AllUnbranchingPaths(logging::Logger &logger, size_t threads, AssemblyGraph &graph);
    void MergePathsToEdges(logging::Logger &logger, size_t threads, AssemblyGraph &graph, const std::vector<GraphPath> &paths);
    void MergeAllToEdges(logging::Logger &logger, size_t threads, AssemblyGraph &graph);
    Vertex &MergePathSPG(const GraphPath &path, AssemblyGraph &graph);
    void MergeAllSPG(logging::Logger &logger, size_t threads, AssemblyGraph &graph);
    AssemblyGraph LoadSupregraphFromGFA(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path);
}
