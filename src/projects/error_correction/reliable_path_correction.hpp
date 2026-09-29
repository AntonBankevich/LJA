#pragma once

#include "reliable_vertex_correction.hpp"

namespace ag {
//    Corrects the tips/segments produced by AbstractReliableVertexSplittingCorrectionAlgorithm by
//    searching the graph itself for an alternative route between the same (reliable) endpoints --
//    it never looks at coverage, only at whether a candidate edge leads to a vertex that has been
//    explicitly marked VertexReliability::unreliable (vertices marked reliable or left unknown are
//    both fair game, since "unknown" is the normal state of the interior of a bulge/tip). If exactly
//    one such alternative exists, it is applied directly; if several are found, the one closest (by
//    edit distance) to the original sequence is used, the same tournament(...) comparison
//    tournament_correction.cpp uses to disambiguate coverage-based alternatives. Tips additionally
//    need their candidate's length trimmed down to the best-aligning prefix, since an open-ended
//    search rarely produces a candidate of exactly the right length.
//    Whatever candidate is picked (uniquely found, or chosen by tournament) is only applied if it is
//    within max_edit_distance of the original sequence -- this rejects "corrections" that are simply
//    wrong even though they were the only/best candidate the graph search turned up.
    class ReliablePathCorrector : public AbstractReliableVertexSplittingCorrectionAlgorithm {
    private:
        size_t max_edit_distance;
    protected:
        std::pair<ag::GraphPath, std::string> correctSegment(const ag::GraphPath &segment) override;
        std::pair<ag::GraphPath, std::string> correctTip(const ag::GraphPath &tip) override;

    public:
        explicit ReliablePathCorrector(size_t max_edit_distance = 100) :
                AbstractReliableVertexSplittingCorrectionAlgorithm("ReliablePathCorrector"),
                max_edit_distance(max_edit_distance) {}
    };
}
