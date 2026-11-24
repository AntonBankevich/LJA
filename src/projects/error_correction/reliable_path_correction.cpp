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
//        alternatives, on the distances the search computed); keeps the original if the search found
//        nothing or the tournament is inconclusive.
        std::pair<GraphPath, std::string> chooseAlternative(const GraphPath &original,
                                                             std::vector<AlignedAlternative> &alternatives,
                                                             const Tournament &tournament) {
            if (alternatives.empty())
                return {original, ""};
            size_t winner = 0;
            std::string message = "s";
            if (alternatives.size() > 1) {
                winner = tournament.tournament(original.truncLen(), alternatives);
                if (winner == size_t(-1))
                    return {original, ""};
                message = "m";
            }
            return {std::move(alternatives[winner].path), message};
        }
    }

    std::pair<GraphPath, std::string> ReliablePathCorrector::correctSegment(const GraphPath &segment) {
        VERIFY(segment.startClosed());
        VERIFY(segment.endClosed());
        const size_t inner_len = std::max(segment.truncLen(), segment.getFinish().size()) - segment.getFinish().size();
        std::vector<AlignedAlternative> alternatives = FindAlignedBulgeAlternatives(
                segment, search_params.divergence.at(inner_len), isTraversable, search_params);
        return chooseAlternative(segment, alternatives, tournament);
    }

    std::pair<GraphPath, std::string> ReliablePathCorrector::correctTip(const GraphPath &tip) {
        if (tip.getFinish().getInnerId() == 690616 || tip.getFinish().size() == 6772) {
            std::cout << "found" << std::endl;
        }
        std::vector<AlignedAlternative> alternatives = FindAlignedTipAlternatives(
                tip, search_params.divergence.at(tip.truncLen()), isTraversable, search_params);
        return chooseAlternative(tip, alternatives, tournament);
    }
}
