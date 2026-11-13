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

//        A tip is its anchor (the reliable vertex it starts at) and the read's sequence after the anchor. path runs from
//        the anchor to the read's end; the result is the minimal path from the anchor covering that sequence, with no
//        cut at the anchor (sample segments are not considered). Empty if no sequence remains after the anchor.
        GraphPath TipAfterAnchor(GraphPath path) {
            path.setCutLeft(0);
            size_t cut = path.rightCut();
            path.setCutRight(0);
            while (!path.empty() && path.backEdge().truncSize() <= cut) {
                cut -= path.backEdge().truncSize();
                path.pop_back();
            }
            path.setCutRight(cut);
            return path;
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
                VERIFY(!corrected.valid());
                GraphPath tip = TipAfterAnchor(read_path.subPath(region.first, region.second + 1).RC());
                if (tip.empty()) {
//                    Nothing of the read before the anchor: kept as it is. The read's right cut, if the anchor is its last
//                    vertex, comes back with that vertex after the loop.
                    GraphPath kept = read_path.subPath(region.first, region.second + 1);
                    kept.setCutRight(0);
                    corrected += kept;
                } else {
                    auto [tc, tip_message] = correctTip(tip);
                    VERIFY(tc.getStart() == tip.getStart());
                    if(!tip_message.empty()) {
                        messages.emplace_back("i" + tip_message);
                        VERIFY(tip != tc);
                    }
                    corrected += tc.RC();
                }
            } else if (region.second == read_path.lastPosition()) {
                GraphPath tail = read_path.subPath(region.first - 1, region.second);
                GraphPath tip = TipAfterAnchor(tail);
                if (tip.empty()) {
//                    Nothing of the read after the anchor: kept as it is.
                    corrected += tail;
                } else {
                    auto [tc, tip_message] = correctTip(tip);
                    VERIFY(tc.getStart() == tip.getStart());
                    if(!tip_message.empty()) {
                        messages.emplace_back("o" + tip_message);
                        VERIFY(tip != tc);
                    }
//                    If the anchor is the read's first vertex, the read's left cut lies in it.
                    tc.setCutLeft(tail.leftCut());
                    corrected += tc;
                }
            } else {
                GraphPath segment = read_path.subPath(region.first - 1, region.second + 1);
                VERIFY(segment.getStart().reliability == VertexReliability::reliable);
                VERIFY(segment.getFinish().reliability == VertexReliability::reliable);
//                In a supregraph a read may start or end inside a reliable vertex, so a segment reaching the read's end
//                inherits its cut. The segment is corrected closed (it keeps its end vertices); the left cut is put back,
//                the right one comes back with the read's last vertex, appended after the loop.
                GraphPath closed = segment;
                closed.setCutLeft(0);
                closed.setCutRight(0);
                auto [sc, segment_message] = correctSegment(closed);
                VERIFY(sc.getStart() == segment.getStart() && sc.getFinish() == segment.getFinish());
                if(!segment_message.empty()) {
                    messages.emplace_back("b" + segment_message);
                    VERIFY(closed != sc);
                }
                sc.setCutLeft(segment.leftCut());
                corrected += sc;
            }
        }
        if (prev_pos != read_path.endPosition()) {
            corrected += read_path.subPath(prev_pos, read_path.lastPosition());
        }
        if(messages.empty())
            return "";
        read_path = std::move(corrected);
        return join("_", messages);
    }
}
