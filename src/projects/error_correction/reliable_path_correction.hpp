#pragma once

#include "reliable_vertex_correction.hpp"
#include "path_alternatives.hpp"

namespace ag {
//    Corrects the tips/segments produced by AbstractReliableVertexSplittingCorrectionAlgorithm by
//    searching the graph itself for alternative routes aligned to the read's sequence
//    (FindAlignedTipAlternatives / FindAlignedBulgeAlternatives, tuned by search_params) --
//    it never looks at coverage, only at whether a candidate edge leads to a vertex that has been
//    explicitly marked VertexReliability::unreliable (vertices marked reliable or left unknown are
//    both fair game, since "unknown" is the normal state of the interior of a bulge/tip). If exactly
//    one such alternative exists, it is applied directly; if several are found, the one closest (by
//    edit distance) to the original sequence is chosen by the tournament (see Tournament; relaxed for
//    noisy reads by default), with the distances the search computed.
//    How divergent a correction may be is set by search_params.divergence: a tip or bulge whose query
//    (sequence) has length L is only corrected to alternatives within divergence.at(L) of it, and the
//    search never goes further. Callers can loosen it over time (a new corrector per pass), so that
//    reads of low quality are only corrected once the graph has fewer options left.
    class ReliablePathCorrector : public AbstractReliableVertexSplittingCorrectionAlgorithm {
    private:
        AlignedSearchParams search_params;
        Tournament tournament;
    protected:
        std::pair<ag::GraphPath, std::string> correctSegment(const ag::GraphPath &segment) override;
        std::pair<ag::GraphPath, std::string> correctTip(const ag::GraphPath &tip) override;

    public:
//        The tournament's distance limit is the search's divergence limit.
        explicit ReliablePathCorrector(AlignedSearchParams search_params = {30, {30, 0.03}, 5, 100000},
                                       double divergence_factor = 2, double min_diff_factor = 0.75) :
                AbstractReliableVertexSplittingCorrectionAlgorithm("ReliablePathCorrector"),
                search_params(search_params),
                tournament(divergence_factor, min_diff_factor, search_params.divergence) {}
    };
}
