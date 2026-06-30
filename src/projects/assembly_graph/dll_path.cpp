#include "dll_path.hpp"

#include "dll_path_storage.hpp"

ag::AlignmentFragment ag::AlignmentFragment::InnerFragment(Segment<Contig> seg, Vertex &vertex, size_t from) {
    return {seg, vertex.getId(), {}, from, vertex.size() - from - seg.size()};
}

ag::AlignmentFragment ag::AlignmentFragment::EdgeSegment(Segment<Contig> seg, Edge &edge, size_t from) {
    return {seg, {}, edge.getId(), from, edge.fullSize() - from - seg.size()};
}

bool ag::AlignmentFragment::isLinked(const AlignmentFragment &other) const {
    if (!valid() || !other.valid() || isInnerFragment() || other.isInnerFragment()) return false;
    VERIFY(isEdgeSegment());
    VERIFY(other.isEdgeSegment());
    return edge->getFinish() == other.edge->getStart() && cut_right == 0 && other.cut_left == 0 &&
           seg.contig() == other.seg.contig() && seg.right == other.seg.left + edge->getFinish().size();
}

bool ag::AlignmentFragment::checkOverlap(const AlignmentFragment &other) const {
    return seg.inter(other.seg) && vertex == other.vertex && edge == other.edge && cut_left - seg.left == other.cut_left - other.seg.left;;
}

bool ag::AlignmentFragment::checkSeqMatch() const {
    if (isEdgeSegment())
        return seg.contig().getSeq().Subseq(seg.left, seg.right) == edge->fullSeq().Subseq(cut_left, edge->fullSize() - cut_right);
    else if (isInnerFragment()) {
        return seg.contig().getSeq().Subseq(seg.left, seg.right) == vertex->getSeq().Subseq(cut_left, vertex->size() - cut_right);
    } else {
        VERIFY(false);
    }
    return true;
}

ag::AlignmentFragment ag::AlignmentFragment::merge(const AlignmentFragment &other) const {
    return {seg.unite(other.seg), vertex, edge, std::min(cut_left, other.cut_left), std::min(cut_right, other.cut_right)};
}

ag::DLLAlignmentPath::DLLAlignmentPath(const std::vector<AlignmentChain<Contig, Edge>> &als) {
    for (const AlignmentChain<Contig, Edge> &al : als) {
        // if (!empty() && !(back().edge.valid() && back().edge->getFinish() == al.seg_to.contig().getStart() && back().cut_right == 0 && al.seg_to.cutLeft() == 0)) {
        //     size_t prev_end = back().seg.right;
        //     size_t next_start = al.seg_from.left;
        //     push_back(AlignmentFragment::Gap({al.seg_from.contig(), std::min(prev_end, next_start), std::max(prev_end, next_start)}));
        // }
        push_back(AlignmentFragment::EdgeSegment(al.seg_from.extendRight(al.seg_to.contig().getStart().size()),
            al.seg_to.contig(), al.seg_to.cutLeft()));
    }
}

ag::DLLAlignmentPath::DLLAlignmentPath(Contig &contig, const GraphPath &path) {
    size_t shift = 0;
    for (PathPosition pp = path.firstPosition(); pp < path.lastPosition(); ++pp) {
        size_t left = pp == path.firstPosition() ? path.leftCut() : 0;
        size_t right = pp == path.lastPosition() - 1 ? path.rightCut() : 0;
        size_t size = pp.nextEdge().fullSize() - left - right;
        push_back(AlignmentFragment::EdgeSegment(contig.segment(shift, shift + size), pp.nextEdge(), left));
        shift += pp.nextEdge().rc().truncSize();
    }
}

