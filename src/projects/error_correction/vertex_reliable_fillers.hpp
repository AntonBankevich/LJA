#pragma once

#include "assembly_graph/assembly_graph.hpp"

namespace ag {
//    Coverage-threshold analogue of dbg's CoverageReliableFiller (see reliable_fillers.hpp), but for
//    the vertex-level VertexReliability marker consumed by AbstractReliableVertexSplittingCorrectionAlgorithm
//    (see reliable_vertex_correction.hpp). Coverage below threshold only means "not enough evidence
//    yet" (VertexReliability::unknown), not "known bad" (VertexReliability::unreliable), so this
//    filler only ever produces the reliable mark, same as the edge filler it mirrors.
    class VertexCoverageReliableFiller {
    private:
        double threshold;
        size_t normalizeReliability(AssemblyGraph &graph) const ;

    public:
        explicit VertexCoverageReliableFiller(double threshold) : threshold(threshold) {}

//        Marks vertices with coverage >= threshold as reliable; leaves everything else untouched.
        size_t fill(AssemblyGraph &graph) const;

//        Resets every vertex to VertexReliability::unknown first, then fills -- use this between
//        passes of an iterative pipeline where coverage estimates keep changing.
        size_t refill(AssemblyGraph &graph) const;
        size_t loggedRefill(logging::Logger &logger, AssemblyGraph &graph) const;
    };
}
