#include "vertex_reduction.hpp"

#include "assembly_graph/assembly_graph.hpp"

using namespace ag;

Segment<Vertex> ExtendToSize(Segment<Vertex> seg, size_t min_size) {
    if(seg.size() < min_size) {
        size_t left_ext = seg.cutLeft();
        size_t right_ext = seg.cutRight();
        if(left_ext *2 < min_size + 1 - seg.size()) {
            right_ext = min_size - left_ext;
        } else if(left_ext *2 < min_size + 1 - seg.size()) {
            left_ext = min_size - right_ext;
        } else {
            left_ext = (min_size + 1 -seg.size()) / 2;
            right_ext = (min_size + 1 - seg.size()) / 2;
        }
        VERIFY(seg.size() + left_ext + right_ext >= min_size);
        return {seg.contig(), seg.left - std::min(seg.left, left_ext), std::min(seg.right + right_ext, seg.contig().size())};
    } else {
        return seg;
    }
}

std::vector<Segment<Vertex>> ExtendSegments(std::vector<Segment<Vertex>> &segs, size_t min_alignment) {
    std::vector<Segment<Vertex>> extended_segs;
    for(const Segment<Vertex> &seg : segs) {
        extended_segs.emplace_back(ExtendToSize(seg, min_alignment));
    }
    return extended_segs;
}

void ensureOverlap(std::unordered_map<VertexId, std::pair<size_t, size_t>> reduction, Edge &edge) {
    if (edge.isPrefix()) {
        reduction[edge.getFinish().getId()].second = std::max(reduction[edge.getFinish().getId()].second,reduction[edge.getFinish().getId()].second);
        reduction[edge.getStart().getId()].first = std::min(reduction[edge.getStart().getId()].first, reduction[edge.getFinish().getId()].first);
    } else if (edge.isSuffix()) {
        reduction[edge.getStart().getId()].first = std::min(reduction[edge.getStart().getId()].first, reduction[edge.getFinish().getId()].first + edge.rc().truncSize());
        reduction[edge.getFinish().getId()].second = std::max(reduction[edge.getFinish().getId()].second, reduction[edge.getStart().getId()].second - edge.rc().truncSize());
    }
}

void CollapseCoreVertexToPrefix(Vertex &v, size_t min_overlap, std::unordered_map<VertexId, std::pair<size_t, size_t>> &reduction) {
    reduction[v.getId()] = {0, min_overlap};
    reduction[v.rc().getId()] = {v.size() - min_overlap, v.size()};
    for (Edge &edge : v) {
        reduction[edge.getFinish().getId()].first = 0;
        reduction[edge.getFinish().rc().getId()].second = edge.getFinish().size();
    }
}

std::unordered_map<ag::VertexId, Segment<ag::Vertex>> ConstructReduction(ag::AssemblyGraph &graph, size_t min_overlap, size_t max_repeat){
    std::unordered_map<VertexId, std::pair<size_t, size_t>> reduction;
    std::vector<VertexId> list = oneline::map(graph.verticesUnique().begin(), graph.verticesUnique().end(),
                                              IdTransformer<Vertex>());
    std::sort(list.begin(), list.end(), [](const VertexId &vid1, const VertexId &vid2) {return vid1->size() > vid2->size() || (vid1->size() == vid2->size() && vid1 > vid2);});
    VERIFY(list.empty() || list.front()->size() >= list.back()->size());
    for(Vertex &vertex : graph.vertices()) {
        reduction[vertex.getId()] = {0, vertex.size()};
    }
    for (Edge &edge : graph.edges()) {
        if (edge.isPrefix()) {
            VERIFY(edge.getStart().size() >= min_overlap);
            reduction[edge.getFinish().getId()].first = edge.getStart().size() - min_overlap;
        } else if (edge.isSuffix()) {
            reduction[edge.getStart().getId()].second = edge.rc().truncSize() + min_overlap;
        }
    }
    for (Vertex &v: graph.vertices()) {
        if (v.isCore() && v.outDeg() == 1 && v.inDeg() > 1 && v.front().isPrefix()) {
            VERIFY(reduction[v.getId()].first == 0 && reduction[v.getId()].second == v.size());
            CollapseCoreVertexToPrefix(v, min_overlap, reduction);
        }
    }
    for (Vertex &v: graph.verticesUnique()) {
        if (v.isCore() && v.inDeg() >= 2 && v.outDeg() >= 2 && v.size() < max_repeat) {
            size_t left = std::min(v.rc().front().getFinish().size(), v.rc().back().getFinish().size());
            size_t right = std::min(v.front().getFinish().size(), v.back().getFinish().size());
            if (left > right && left > 2 * v.size()) {
                CollapseCoreVertexToPrefix(v.rc(), min_overlap, reduction);
            } else if (right >= left && right > 2 * v.size()) {
                CollapseCoreVertexToPrefix(v, min_overlap, reduction);
            }
        }
    }
    for (Vertex &vertex : graph.verticesUnique()) {
        if (!vertex.isJunction() && !vertex.front().getFinish().isJunction()) {
            reduction[vertex.getId()] = {0, vertex.size()};
        }
    }
    for (VertexId vid : list) {
        for (Edge &edge : *vid)
            ensureOverlap(reduction, edge);
    }
    for (Edge &edge : graph.edges()) {
        auto l = reduction[edge.getStart().getId()];
        auto r = reduction[edge.getFinish().getId()];
        VERIFY(l.first <= r.first + edge.rc().truncSize());
        VERIFY(l.second <= r.second + edge.rc().truncSize());
        VERIFY(l.second >= r.first + edge.rc().truncSize() + min_overlap);
        VERIFY(l.second <= r.first + edge.rc().truncSize() + 30000);
    }

    for (Vertex &v : graph.vertices()) {
        auto l = reduction[v.getId()];
        auto r = reduction[v.rc().getId()];
        VERIFY(l.first + r.second == v.size());
        VERIFY(l.second + r.first == v.size());
    }

    // for(VertexId vid : list) {
    //     if(vid->isCore()) {
    //         if (vid->inDeg() == 1) {
    //             reduction[vid] = {*vid, vid->size(), vid->size()};
    //             reduction[vid->rc().getId()] = {vid->rc(), 0, 0};
    //         }
    //         continue;
    //     }
    //     size_t left_cut = vid->size();
    //     size_t right_cut = vid->size();
    //     for(Edge &edge : *vid) {
    //         if(edge.isPrefix()) {
    //             Segment<Vertex> seg = reduction.at(edge.getFinish().getId());
    //             left_cut = std::min(left_cut, seg.left);
    //         } else {
    //             VERIFY(edge.isSuffix());
    //             right_cut = std::min(right_cut, edge.getFinish().size());
    //         }
    //     }
    //     for(Edge &edge : vid->rc()) {
    //         if(edge.isPrefix()) {
    //             Segment<Vertex> seg = reduction.at(edge.getFinish().getId());
    //             right_cut = std::min(right_cut, seg.left);
    //         } else {
    //             VERIFY(edge.isSuffix());
    //             left_cut = std::min(left_cut, edge.getFinish().size());
    //         }
    //     }
    //     if(left_cut == vid->size()) left_cut = 0;
    //     if(right_cut == vid->size()) right_cut = 0;
    //     size_t left = left_cut;
    //     size_t right = vid->size() - right_cut;
    //     VERIFY(left < right || vid->isOuter());
    //     Segment<Vertex> res(*vid, std::min(left, right), std::max(left, right));
    //     res = ExtendToSize(res, min_size);
    //     reduction[vid] = res;
    //     reduction[vid->rc().getId()] = res.RC();
    // }
    std::unordered_map<VertexId, Segment<Vertex>> res;
    for (Vertex &vertex : graph.vertices())
        res[vertex.getId()] = {vertex, reduction[vertex.getId()].first, reduction[vertex.getId()].second};
    return std::move(res);
}
