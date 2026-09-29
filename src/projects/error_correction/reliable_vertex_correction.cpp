#include "reliable_vertex_correction.hpp"
#include "common/string_utils.hpp"

namespace ag {
    namespace {
        std::vector<std::pair<PathPosition, PathPosition>> findUnreliableRegions(const GraphPath &read_path) {
            std::vector<std::pair<PathPosition, PathPosition>> res;
            PathPosition pos = read_path.firstPosition();
            while(pos != read_path.endPosition()) {
                if (pos.getVertex().reliability == VertexReliability::unreliable) {
                    PathPosition left = pos;
                    PathPosition right = pos;
                    while(left != read_path.firstPosition() && (left-1).getVertex().reliability != VertexReliability::reliable) {
                        --left;
                    }
                    while(right != read_path.lastPosition() && (right+1).getVertex().reliability != VertexReliability::reliable) {
                        ++right;
                    }
                    res.emplace_back(left, right);
                    pos = right;
                }
                ++pos;
            }
            return res;
        }
    }

    std::string AbstractReliableVertexSplittingCorrectionAlgorithm::correctRead(const std::string &name,
                                                                                 ag::GraphPath &read_path) {
        std::vector<std::pair<PathPosition, PathPosition>> unreliable_regions = findUnreliableRegions(read_path);
        std::vector<std::string> messages;
        GraphPath corrected;
        if (unreliable_regions.empty() || (unreliable_regions.size() == 1 &&
            unreliable_regions.front().first == read_path.firstPosition() &&
            unreliable_regions.front().second == read_path.lastPosition())) {
            return "";
        }
        PathPosition prev_pos = read_path.firstPosition();
        for (std::pair<PathPosition, PathPosition> region : unreliable_regions) {
            if (region.first != read_path.firstPosition() && region.first - 1 != prev_pos) {
                corrected += read_path.subPath(prev_pos, region.first - 1);
            }
            prev_pos = region.second + 1;
            if (region.first == read_path.firstPosition()) {
                GraphPath tip = read_path.subPath(region.first, region.second + 1).RC();
                auto [tc, tip_message] = correctTip(tip);
                VERIFY(tc.getStart() == tip.getStart());
                if(!tip_message.empty()) {
                    messages.emplace_back("i" + tip_message);
                    VERIFY(tip != tc);
                }
                VERIFY(!corrected.valid());
                corrected += tc.RC();
            } else if (region.second == read_path.lastPosition()) {
                GraphPath tip = read_path.subPath(region.first - 1, region.second);
                auto [tc, tip_message] = correctTip(tip);
                VERIFY(tc.getStart() == tip.getStart());
                if(!tip_message.empty()) {
                    messages.emplace_back("o" + tip_message);
                    VERIFY(tip != tc);
                }
                corrected += tc;
            } else {
                GraphPath segment = read_path.subPath(region.first - 1, region.second + 1);
                VERIFY(segment.getStart().reliability == VertexReliability::reliable);
                VERIFY(segment.getFinish().reliability == VertexReliability::reliable);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    auto [sc, segment_message] = correctSegment(segment);
                VERIFY(sc.getStart() == segment.getStart() && sc.getFinish() == segment.getFinish());
                if(!segment_message.empty()) {
                    messages.emplace_back("b" + segment_message);
                    VERIFY(segment != sc);
                }
                corrected += sc;
            }
        }
        if (prev_pos != read_path.endPosition()) {
            corrected += read_path.subPath(prev_pos, read_path.lastPosition());
        }
        if(messages.empty())
            return "";
        VERIFY(!corrected.frontEdge().isPrefix());
        VERIFY(!corrected.backEdge().isSuffix());
        read_path = std::move(corrected);
        return join("_", messages);
    }
}
