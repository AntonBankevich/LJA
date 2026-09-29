#pragma once

#include "assembly_graph/graph_paths.hpp"
#include <functional>
#include <unordered_map>

namespace ag {
//    Generic (graph-agnostic) building blocks for finding alternative paths through a region of the
//    graph and choosing between candidate corrections purely by sequence similarity. `isTraversable`
//    decides which edges the search may walk through -- callers that trust edge coverage plug in a
//    coverage-based predicate, callers that only trust the vertex-level VertexReliability marker (see
//    reliable_vertex_correction.hpp) plug in a predicate based on that instead.
    std::unordered_map<Vertex *, size_t> findReachable(Vertex &start,
                                                        const std::function<bool(const Edge &)> &isTraversable,
                                                        size_t max_dist);

//    Alternative paths from path.getStart() to path.getFinish(), each within max_diff of path's length.
//    path itself is excluded from the result.
    std::vector<GraphPath> FindAlternativeSegments(const GraphPath &path, size_t max_diff,
                                                    const std::function<bool(const Edge &)> &isTraversable);

//    Alternative extensions from path.getStart() of roughly path's length (open-ended, unlike
//    FindAlternativeSegments which additionally fixes the finish vertex). Results that start with path
//    itself are excluded.
    std::vector<GraphPath> FindAlternativeTips(const GraphPath &path, size_t max_diff,
                                                const std::function<bool(const Edge &)> &isTraversable);

//    Index of the candidate closest (by edit distance) to original, or size_t(-1) if no candidate is
//    unambiguously closer than all the others.
    size_t tournament(const Sequence &original, const std::vector<Sequence> &candidates);

//    Truncates al to the prefix that best aligns to seq -- used to find the correct cutoff for a tip
//    alternative, whose natural length rarely matches the original tip's length exactly.
    std::pair<GraphPath, size_t> bestAlignmentPrefix(const GraphPath &al, const Sequence &seq, size_t max_diff);
}
