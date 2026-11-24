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

//    An alternative found by FindAlignedTipAlternatives / FindAlignedBulgeAlternatives and the edit distance of its
//    truncSeq() to the original's.
    struct AlignedAlternative {
        GraphPath path;
        size_t cost;
    };

//    Limit on the alignment cost of an explored path after L letters of it: base + rate * L (default: none).
    struct DivergenceLimit {
        size_t base = size_t(-1) / 2;
        double rate = 0;

        size_t at(size_t len) const { return base + size_t(rate * double(len)); }
    };

//    Tuning of FindAlignedTipAlternatives / FindAlignedBulgeAlternatives (see below).
    struct AlignedSearchParams {
        size_t band_slack = 30;
        DivergenceLimit divergence;
        size_t max_choices = 5;
        size_t max_steps = 100000;
    };

//    Read-guided searches for alternatives to a tip or a bulge of a read, starting at its first vertex (no left cut).
//    Paths from that vertex are explored best-first (Dijkstra on the tree of paths, following edges allowed by
//    isTraversable). Each explored path carries a band: its costs of alignment to the query prefixes around its best
//    end, within band_slack of the band's minimum and within divergence.at(L) for the path's length L (checked after
//    every letter; a path whose band empties is dropped). Extending a path by an edge runs the edit-distance DP along
//    the edge; the band minimum, the path's priority, never decreases along a path.
//    Many paths reach the same vertex with the same best query position. Of those, a path is extended only if it is
//    among the max_choices with the lowest band minimum, or if at some query position it is better than all the
//    others kept: near-best paths still matter, as they can make the tournament between alternatives ambiguous.
//    The search stops when the next path costs more than max_cost or the best alternative + band_slack, and returns
//    the alternatives within both limits.
//    Cost 0 (the original itself, or a path spelling exactly its sequence) is never the best and is not returned.
//    Both return nothing if they give up after max_steps paths.

//    Tip: the query is tip.truncSeq(). An alternative may end anywhere and is cut where its alignment to the whole query
//    ends. A path and its extensions count as one route (its cheapest alignment).
    std::vector<AlignedAlternative> FindAlignedTipAlternatives(const GraphPath &tip, size_t max_cost,
                                                               const std::function<bool(const Edge &)> &isTraversable,
                                                               const AlignedSearchParams &params = {});

//    Bulge (closed, ending at vertex F): past the point where F's sequence starts, every route to F just spells F, so
//    only the parts before it (inner parts) are aligned. An alternative ends where F's sequence starts, inside an edge
//    found by walking back from F, and a fixed route then leads to F. If the start overlaps F along a route, the route
//    spells a suffix of F, and the cost is the length difference. The divergence limit applies to the alignment of
//    the inner parts.
    std::vector<AlignedAlternative> FindAlignedBulgeAlternatives(const GraphPath &bulge, size_t max_cost,
                                                                 const std::function<bool(const Edge &)> &isTraversable,
                                                                 const AlignedSearchParams &params = {});

//    Chooses among candidate corrections of an original sequence by edit distance. The candidate closest to the
//    original (the winner, at distance dist[w]) wins if dist[w] is below the distance limit limit.at(|original|)
//    and every other candidate i below the limit is behind it:
//        dist[w] <= dist[i] - diff(w, i) * factor,   factor = max(min_diff_factor, 1 - divergence_factor * d),
//    where diff(w, i) is the distance between the two candidates and d = dist[w] / |original| estimates the original's
//    error rate. With factor 1 the winner has to lie exactly between the original and i (dist[i] = dist[w] + diff),
//    the strictest possible rule. For a noisy original, errors in it can match some of the differences between the
//    candidates by chance, so a fraction (about divergence_factor * d) of diff is allowed to go unexplained, but the
//    factor never drops below min_diff_factor. divergence_factor = 0 gives the strict rule.
    class Tournament {
    public:
        explicit Tournament(double divergence_factor = 2, double min_diff_factor = 0.75,
                            DivergenceLimit limit = {20, 0.01});

//        Index of the winner among candidates, or size_t(-1) if there is none.
        size_t tournament(const Sequence &original, const std::vector<Sequence> &candidates) const;

//        The same for alternatives found by the aligned searches: their costs are used as the distances to the
//        original (of size original_size); only the distances between the winner and the others are computed.
        size_t tournament(size_t original_size, const std::vector<AlignedAlternative> &alternatives) const;

    private:
        double divergence_factor;
        double min_diff_factor;
        DivergenceLimit limit;

        size_t maxDist(size_t original_size) const;
        size_t tournamentOnDistances(size_t original_size, const std::vector<size_t> &dists,
                                     const std::vector<Sequence> &candidates) const;
    };

//    Truncates al to the prefix that best aligns to seq -- used to find the correct cutoff for a tip
//    alternative, whose natural length rarely matches the original tip's length exactly.
    std::pair<GraphPath, size_t> bestAlignmentPrefix(const GraphPath &al, const Sequence &seq, size_t max_diff);
}
