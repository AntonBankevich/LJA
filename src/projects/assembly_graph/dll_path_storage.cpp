#include "dll_path_storage.hpp"

#include "dll_path.hpp"
#include "sequences/contigs.hpp"

namespace ag {
    class Edge;
}

void ag::DLLAlignmentStorage::addContig(Contig &contig, DLLAlignmentPath &&p) {
    DLLAlignmentPath &path = (alignments[contig.getInnerId()] = std::move(p));
    for (DLLPosition pos = path.begin(); pos != path.end(); ++pos) {
        if (pos->isEdgeSegment()) {
            edge_map[pos->edge].emplace_back(pos);
        } else if (pos->isInnerFragment()) {
            vertex_map[(*pos).vertex].emplace_back(pos);
        } else {
            VERIFY(false);
        }
    }
}

ag::DLLAlignmentStorage::DLLPosition ag::DLLAlignmentStorage::insertBefore(Edge &edge, DLLPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.edge== edge.getId());
    DLLPosition new_pos = pos.insertBefore(fragment);
    edge_map[edge.getId()].emplace_back(new_pos);
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPosition ag::DLLAlignmentStorage::insertAfter(Edge &edge, DLLPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.edge== edge.getId());
    DLLPosition new_pos = pos.insertAfter(fragment);
    edge_map[edge.getId()].emplace_back(new_pos);
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPosition ag::DLLAlignmentStorage::insertBefore(Vertex &vertex, DLLPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.vertex== vertex.getId());
    DLLPosition new_pos = pos.insertBefore(fragment);
    vertex_map[vertex.getId()].emplace_back(new_pos);
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPosition ag::DLLAlignmentStorage::insertAfter(Vertex &vertex, DLLPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.vertex== vertex.getId());
    DLLPosition new_pos = pos.insertAfter(fragment);
    vertex_map[vertex.getId()].emplace_back(new_pos);
    return new_pos;
}

ag::DLLAlignmentPath & ag::DLLAlignmentStorage::addContig(Contig new_contig,
                                                          std::vector<AlignmentChain<Contig, Edge>> &als) {
    contigs.push_back(std::move(new_contig));
    Contig &contig = contigs.back();
    contigs.push_back(contig.RC());
    Contig &rc_contig = contigs.back();
    std::vector<AlignmentChain<Contig, Edge>> forward;
    std::vector<AlignmentChain<Contig, Edge>> backward;
    for (AlignmentChain<Contig, Edge> &chain : als) {
        size_t k = chain.seg_to.contig().getStart().size();
        forward.emplace_back(contig, chain.seg_to.contig(), chain.seg_from.left, chain.seg_to.left, chain.seg_from.size());
        backward.emplace_back(rc_contig, chain.seg_to.contig().rc(), rc_contig.fullSize() - chain.seg_from.right - k, chain.seg_to.contig().truncSize() - chain.seg_to.right, chain.seg_from.size());
    }
    backward = {backward.rbegin(), backward.rend()};
    addContig(contig, std::move(forward));
    addContig(rc_contig, std::move(backward));
    return alignments[contig.getInnerId()];
}

void ag::DLLAlignmentStorage::fireDeleteVertex(Vertex &v) {
    if (!hasRecords(v))
        return;
    for (DLLPosition &pos : vertex_map.at(v.getId())) {
        VERIFY(pos.valid());
        VERIFY(pos->isInnerFragment());
        pos.erase();
    }
    vertex_map.erase(v.getId());
}

void ag::DLLAlignmentStorage::fireDeleteEdge(Edge &e) {
    if (!hasRecords(e))
        return;
    for (DLLPosition &pos : edge_map.at(e.getId())) {
        if (!pos.valid()) continue;
        if (!pos.extracted()) {
            if (pos->edge->isSuffix()) {
                DLLPosition cur_pos = pos;
                while (cur_pos.prev()->isLinked(*cur_pos) && cur_pos->edge->isSuffix()) {
                    cur_pos = cur_pos.prev();
                    cur_pos.next().extract();
                }
                if (!cur_pos.prev()->isLinked(*cur_pos) && cur_pos->edge->isSuffix()) {
                    insertBefore(cur_pos->edge->getStart(), cur_pos,
                        AlignmentFragment::InnerFragment(cur_pos->seg, cur_pos->edge->getStart(), cur_pos->cut_left));
                    cur_pos.extract();
                }
            } else if (pos->edge->isPrefix()) {
                DLLPosition cur_pos = pos;
                while (cur_pos->isLinked(*cur_pos.next()) && cur_pos->edge->isPrefix()) {
                    cur_pos = cur_pos.next();
                    cur_pos.prev().extract();
                }
                if (cur_pos->isLinked(*cur_pos.next()) && cur_pos->edge->isPrefix()) {
                    insertAfter(cur_pos->edge->getFinish(), cur_pos,
                        AlignmentFragment::InnerFragment(cur_pos->seg, cur_pos->edge->getFinish(), cur_pos->cut_right));
                    cur_pos.extract();
                }
            }
        }
        pos.erase();
    }
    edge_map.erase(e.getId());
}

void ag::DLLAlignmentStorage::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    if (!hasRecords(e)) return;
    for (DLLPosition &pos : edge_map.at(e.getId())) {
        bool link_left = pos.prev()->isLinked(*pos);
        bool link_right = pos->isLinked(*pos.next());
        if (!link_left && !link_right) {
            insertBefore(v, pos, AlignmentFragment::InnerFragment(pos->seg, v, pos->cut_left));
        } else {
            if (link_left) {
                insertBefore(v.rc().front().rc(), pos, AlignmentFragment::EdgeSegment(pos->seg, v.rc().front().rc(), pos->cut_left));
            }
            if (link_right) {
                insertAfter(v.front(), pos, AlignmentFragment::EdgeSegment(pos->seg, v.front(), pos->cut_left));
            }
        }
        pos.extract();
    }
}

void ag::DLLAlignmentStorage::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    size_t shift = 0;
    for (Edge &edge: path.edges()) {
        if (!hasRecords(edge))
            continue;
        for (DLLPosition &pos : edge_map.at(edge.getId())) {
            if (!pos.valid() || pos.extracted())
                continue;
            size_t cut_left = shift + pos->cut_left;
            DLLPosition next_pos = pos.next();
            DLLPosition prev_pos = pos.prev();
            size_t len = pos->edge->fullSize() - pos->cut_left - pos->cut_right;
            Segment<Contig> seg = pos->seg;
            for (; next_pos.prev()->isLinked(*next_pos) && next_pos->edge->getStart() != path.getFinish(); next_pos = next_pos.next()) {
                len += next_pos->edge->truncSize() - next_pos->cut_right;
                seg = seg.unite(next_pos->seg);
            }
            while (prev_pos.next() != next_pos)
                prev_pos.next().extract();
            VERIFY(len == seg.size());
            size_t cut_right = new_vertex.size() - cut_left - len;
            bool linked = false;
            if (cut_left == 0) {
                AlignmentFragment new_al = AlignmentFragment::EdgeSegment(seg, new_vertex.rc().front().rc(), 0);
                if (prev_pos->isLinked(new_al)) {
                    VERIFY(shift == 0);
                    insertAfter(new_vertex.rc().front().rc(), prev_pos, new_al);
                    linked = true;
                }
            }
            if (cut_right == 0) {
                AlignmentFragment new_al = AlignmentFragment::EdgeSegment(seg, new_vertex.front(), cut_left);
                if (new_al.isLinked(*next_pos)) {
                    insertBefore(new_vertex.front(), next_pos, new_al);
                    linked = true;
                }
            }
            if (!linked) {
                insertBefore(new_vertex, next_pos, AlignmentFragment::InnerFragment(seg, new_vertex, cut_left));
            }
            VERIFY(pos.extracted())
        }
        shift += edge.rc().truncSize();
    }

    shift = 0;
    for (Vertex &v : path.innerVertices()) {
        if (hasRecords(v))
            for (DLLPosition &pos : vertex_map.at(v.getId())) {
                if (!pos.valid())
                    continue;
                insertBefore(new_vertex, pos, AlignmentFragment::InnerFragment(pos->seg, new_vertex, shift + pos->cut_left));
                pos.extract();
            }
        shift += v.rc().front().truncSize();
    }
}

void ag::DLLAlignmentStorage::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    size_t shift = 0;
    for (Edge &edge: path.edges()) {
        if (!hasRecords(edge))
            continue;
        for (DLLPosition &pos : edge_map.at(edge.getId())) {
            if (!pos.valid() || pos.extracted())
                continue;
            size_t extra_shift = 0;
            AlignmentFragment new_fragment = AlignmentFragment::EdgeSegment(pos->seg, new_edge, shift + pos->cut_left);
            DLLPosition last_pos = pos.next();
            for (; last_pos.prev()->isLinked(*last_pos); last_pos = last_pos.next()) {
                AlignmentFragment edge_fragment = AlignmentFragment::EdgeSegment(last_pos->seg, new_edge, shift + extra_shift + last_pos->cut_left);
                new_fragment = new_fragment.merge(edge_fragment);
                extra_shift += last_pos->edge->rc().truncSize();
                last_pos.prev().extract();
            }
            insertBefore(new_edge, last_pos, new_fragment);
        }
        shift += edge.rc().truncSize();
    }
    shift = 0;
    Vertex &new_vertex = new_edge.isSuffix() ? new_edge.getStart() : new_edge.getFinish();
    for (Vertex &v : path.innerVertices()) {
        if (hasRecords(v))
            for (DLLPosition &pos : vertex_map.at(v.getId())) {
                if (!pos.valid())
                    continue;
                insertBefore(new_vertex, pos, AlignmentFragment::InnerFragment(pos->seg, new_vertex, shift + pos->cut_left));
                pos.extract();
            }
        shift += v.rc().front().truncSize();
    }
}

void ag::DLLAlignmentStorage::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &left_al,
    const AlignmentForm &right_al) {
    if (hasRecords(left)) {
        size_t match_size = 0;
        size_t start_pos = left.truncSize() - left_al.queryLength();
        auto col = left_al.columns().begin();
        for (auto col = left_al.columns().begin(); col != left_al.columns().end() &&
                                                   match_size < left_al.queryLength() && (*col).event == CigarEvent::M &&
                                                   left.truncSeq()[start_pos + match_size]==new_edge.truncSeq()[start_pos + match_size]; ++col) {
            match_size++;
        }
        size_t min_right_cut = left_al.queryLength() - match_size;
        for (DLLPosition &pos : edge_map.at(left.getId())) {
            if (!pos.valid()) continue;
            VERIFY(!pos.extracted());
            if (left.truncSize() > pos->cut_left + min_right_cut) {
                size_t extra_cut = std::max(min_right_cut, pos->cut_right) - pos->cut_right;
                AlignmentFragment new_fragment = AlignmentFragment::EdgeSegment(pos->seg.shrinkRightBy(extra_cut), new_edge, pos->cut_left);
                insertBefore(new_edge, pos, new_fragment);
            }
            pos.extract();
        }
        for (DLLPosition &pos : edge_map.at(left.rc().getId())) {
            if (!pos.valid()) continue;
            VERIFY(!pos.extracted());
            if (left.rc().truncSize() > pos->cut_right + min_right_cut) {
                size_t extra_cut = std::max(min_right_cut, pos->cut_left) - pos->cut_left;
                AlignmentFragment new_fragment = AlignmentFragment::EdgeSegment(pos->seg.shrinkLeftBy(extra_cut),
                    new_edge.rc(), new_edge.fullSize() - pos->cut_left - pos->seg.size() + extra_cut);
                insertAfter(new_edge.rc(), pos, new_fragment);
            }
            pos.extract();
        }
    }
}

void ag::DLLAlignmentStorage::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
    if (!hasRecords(edge))
        return;
    for (DLLPosition &pos : edge_map.at(edge.getId())) {
        size_t e_left_cut = 0;
        size_t e_right_cut = edge.truncSize();
        for (Edge &e : split.edges()) {
            e_right_cut -= e.truncSize();
            if (pos->cut_left < edge.fullSize() - e_right_cut - e.getFinish().size() &&
                pos->cut_right > edge.fullSize() - e_left_cut - e.getStart().size()) {
                size_t sub_left_cut = std::max(e_left_cut, pos->cut_left);
                size_t sub_right_cut = std::max(e_right_cut, pos->cut_right);
                size_t len = edge.fullSize() - sub_left_cut - sub_right_cut;
                Segment<Contig> subseg(pos->seg.contig(), pos->seg.left + sub_left_cut - pos->cut_left,
                                       pos->seg.left + sub_left_cut - pos->cut_left + len);
                insertBefore(e, pos, AlignmentFragment::EdgeSegment(subseg, e, sub_left_cut - e_left_cut));
            }
            e_left_cut += e.rc().truncSize();
        }
        pos.extract();
    }
    edge_map.erase(edge.getId());
}

void ag::DLLAlignmentStorage::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    for (Edge &inc: core.incoming()) {
        if (!hasRecords(inc))
            continue;
        for (DLLPosition &pos : edge_map.at(inc.getId())) {
            VERIFY(pos.next()->isEdgeSegment());
            VERIFY(pos.next()->edge->getStart() == core);
            DLLPosition left = pos;
            DLLPosition right = pos.next();
            VERIFY(left->isLinked(*right));
            if (resolution.contains(*left->edge, *right->edge)) {
                Vertex &new_vertex = resolution.get(*left->edge, *right->edge);
                bool is_linked = false;
                if (left.prev()->isLinked(*left)) {
                    AlignmentFragment al = AlignmentFragment::EdgeSegment(left->seg.unite(right->seg),
                                                                          new_vertex.rc().front().rc(), left->cut_left);
                    insertBefore(new_vertex.rc().front().rc(), left, al);
                    is_linked = true;
                }
                if (right->isLinked(*right.next())) {
                    AlignmentFragment al = AlignmentFragment::EdgeSegment(left->seg.unite(right->seg),
                                                                          new_vertex.front(), left->cut_left);
                    insertAfter(new_vertex.front(), right, al);
                    is_linked = true;
                }
                if (!is_linked) {
                    insertBefore(new_vertex, right, AlignmentFragment::InnerFragment(left->seg.unite(right->seg), new_vertex, left->cut_left));
                }
                left.extract();
                right.extract();
            }
        }
    }
}

std::function<std::string(const ag::Vertex &)> ag::DLLAlignmentStorage::getVertexTooltipper() const {
    std::function<std::string(const Vertex &)> res = [this](const Vertex &v)->std::string {
        if (!hasRecords(v)) return "";
        std::stringstream ss;
        for (DLLPosition pos : vertex_map.at(v.getId())) {
            AlignmentFragment f = *pos;
            ss << f << "\n";
        }
        return ss.str();
    };
    return res;
}

std::function<std::string(const ag::Edge &)> ag::DLLAlignmentStorage::getEdgeTooltipper() const {
    std::function<std::string(const Edge &)> res = [this](const Edge &e)->std::string {
        if (!hasRecords(e)) return "";
        std::stringstream ss;
        for (DLLPosition pos : edge_map.at(e.getId())) {
            AlignmentFragment f = *pos;
            ss << f << "\n";
        }
        return ss.str();
    };
    return res;
}

std::function<std::string(const ag::Edge &)> ag::DLLAlignmentStorage::getEdgeColorer() const {
    std::function<std::string(const Edge &)> res = [this](const Edge &e)->std::string {
        if (!hasRecords(e)) return "";
        return "blue";
    };
    return res;
}

std::function<std::string(const ag::Vertex &)> ag::DLLAlignmentStorage::getVertexColorer() const {
    std::function<std::string(const Vertex &)> res = [this](const Vertex &e)->std::string {
        if (!hasRecords(e)) return "";
        return "blue";
    };
    return res;
}

ag::VertexInfo ag::DLLAlignmentStorage::getVertexInfo() const {
    return VertexInfo::Colorer(getVertexColorer()) + VertexInfo::Tooltiper(getVertexTooltipper());
}

ag::EdgeInfo ag::DLLAlignmentStorage::getEdgeInfo() const {
    return EdgeInfo::Colorer(getEdgeColorer()) + EdgeInfo::Tooltiper(getEdgeTooltipper());
}

ag::Printer ag::DLLAlignmentStorage::getPrinter() const {
    return Printer(getVertexInfo(), getEdgeInfo());
}

void ag::DLLAlignmentStorage::print(std::ostream &out) {
    out << "Contigs" << std::endl;
    for (Contig &c: contigs) {
        out << c.getInnerId() << std::endl;
        DLLAlignmentPath &path = alignments.at(c.getInnerId());
        for (AlignmentFragment f: path) {
            out << f << std::endl;
        }
    }
    out << "Vertices" << std::endl;
    for (auto &p: vertex_map) {
        out << p.first << std::endl;
        for (DLLPosition pos: p.second) {
            out << *pos << std::endl;
        }
    }
    out << "Edges" << std::endl;
    for (auto &p: edge_map) {
        out << p.first << std::endl;
        for (DLLPosition pos: p.second) {
            out << *pos << std::endl;
        }
    }
}
