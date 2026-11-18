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
//    Splits unbranching paths (as from AllUnbranchingPaths; given by their indices) into buckets such that no two paths
//    in a bucket share a vertex, up to rc: in a supregraph the paths of one bucket can be merged in parallel without
//    locks. Each bucket lists its indices in increasing order, and the split does not depend on the number of threads.
//    Vertex degrees must be at most 8; there are then at most 31 buckets.
    std::vector<std::vector<size_t>>
    IndependentPathBuckets(logging::Logger &logger, size_t threads, const std::vector<GraphPath> &paths);
    void MergePathsToEdges(logging::Logger &logger, size_t threads, AssemblyGraph &graph, const std::vector<GraphPath> &paths);
    void MergeAllToEdges(logging::Logger &logger, size_t threads, AssemblyGraph &graph);
//    GFA stores neither k nor sample segments (see the TODO at AlignedReadStorage::Load), so the graph gets the k given
//    here; 0 if no code that depends on sample segments will run on it.
    AssemblyGraph LoadSupregraphFromGFA(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path,
                                        size_t k = 0);
}
