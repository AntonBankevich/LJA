#include "reliable_path_correction.hpp"
#include "path_alternatives.hpp"
#include "sequences/edit_distance.hpp"

namespace ag {
    namespace {
//        The only thing this predicate consults is the vertex-level reliability marker: edges into a
//        vertex explicitly marked unreliable are off-limits, edges into a reliable or still-unknown
//        vertex are fair game (unknown is the normal state of a bulge/tip interior).
        bool isTraversable(const Edge &edge) {
            return edge.getFinish().reliability != VertexReliability::unreliable && edge.getStart().reliability != VertexReliability::unreliable;
        }

//        Picks the single alternative if there is only one, otherwise disambiguates by sequence
//        similarity to the original (same rule tournament_correction.cpp uses for coverage-based
//        alternatives); gives up (keeps the original) if the search found nothing, the tournament is
//        inconclusive, or the picked candidate still differs from the original by max_edit_distance
//        or more (a graph-search match is not necessarily a correct one).
        std::pair<GraphPath, std::string> chooseBySequence(const GraphPath &original,
                                                             std::vector<GraphPath> &alternatives,
                                                             size_t max_edit_distance) {
            if (alternatives.empty())
                return {original, ""};
            Sequence old = original.truncSeq();
            GraphPath *chosen;
            std::string message;
            if (alternatives.size() == 1) {
                chosen = &alternatives[0];
                message = "s";
            } else {
                std::vector<Sequence> candidates;
                for (GraphPath &al: alternatives)
                    candidates.push_back(al.truncSeq());
                size_t winner = tournament(old, candidates);
                if (winner == size_t(-1))
                    return {original, ""};
                chosen = &alternatives[winner];
                message = "m";
            }
            if (edit_distance(old, chosen->truncSeq(), max_edit_distance) >= max_edit_distance)
                return {original, ""};
            return {std::move(*chosen), message};
        }
    }

    std::pair<GraphPath, std::string> ReliablePathCorrector::correctSegment(const GraphPath &segment) {
        VERIFY(segment.startClosed());
        VERIFY(segment.endClosed());
        size_t inner_len = std::max(segment.truncLen(), segment.getFinish().size()) - segment.getFinish().size();
        size_t max_diff = std::max<size_t>(30, inner_len * 10 /100);
        std::vector<GraphPath> alternatives = FindAlternativeSegments(segment, max_diff, isTraversable);
        if (alternatives.empty() || (alternatives.size() == 1 && alternatives.front() == segment))
            return {segment, ""};
        return chooseBySequence(segment, alternatives, max_edit_distance);
    }

    std::pair<GraphPath, std::string> ReliablePathCorrector::correctTip(const GraphPath &tip) {
        if (tip.getFinish().getInnerId() == 690616 || tip.getFinish().size() == 6772) {
            std::cout << "found" << std::endl;
        }
        size_t max_diff = std::max<size_t>(100, tip.truncLen() * 3 / 100);
        std::vector<GraphPath> alternatives = FindAlternativeTips(tip, max_diff, isTraversable);
        Sequence old = tip.truncSeq();
        std::vector<GraphPath> trunc_alternatives;
        for (GraphPath &al: alternatives) {
//            An open-ended tip search rarely lands on a candidate of exactly tip's length, so the
//            cutoff has to be found by aligning the candidate to the original tip sequence.
            std::pair<GraphPath, size_t> cut = bestAlignmentPrefix(al, old, max_diff);
            if (cut.second < max_diff)
                trunc_alternatives.emplace_back(std::move(cut.first));
        }
        if (alternatives.empty() || (alternatives.size() == 1 && alternatives.front() == tip))
            return {tip, ""};
        return chooseBySequence(tip, trunc_alternatives, max_edit_distance);
    }
}
