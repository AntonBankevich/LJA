#include "dll_path_storage.hpp"

#include "dll_path.hpp"
#include "sequences/contigs.hpp"
#include <unordered_set>

namespace ag {
    class Edge;
}

void ag::PathDrawer::draw(Contig &contig, const std::vector<VertexId> &vertices, const std::string &event_tag) {
    if (vertices.empty())
        return;
    const std::string &contig_name = contig.getInnerId();
    if (startsWith(contig_name, "-"))
        return;
    std::unordered_set<ConstVertexId> vertex_set(vertices.begin(), vertices.end());
    std::function<std::string(const Vertex &)> colorer = [vertex_set](const Vertex &v) -> std::string {
        return vertex_set.count(v.getId()) ? "orange" : "";
    };
    Printer p = printer + VertexInfo::Colorer(colorer);
    std::experimental::filesystem::path fname = figuresFor(contig_name).nextFile(event_tag);
    p.printDot(fname, Component::neighbourhood(*graph, vertices, radius, max_size));
}

void ag::DLLAlignmentStorage::logPath(Contig &contig, const std::string &source) {
    if (!log_stream)
        return;
    *log_stream << "Change in observed path during " << source << ": " << contig.getInnerId() << ":";
    for (AlignmentFragment &f : alignments.at(contig.getInnerId()))
        *log_stream << " " << f;
    *log_stream << std::endl;
}

void ag::DLLAlignmentStorage::recordTouch(std::unordered_map<Contig *, std::unordered_set<VertexId>> &touched,
                                           const AlignmentFragment &fragment) {
    Contig &contig = fragment.seg.contig();
    std::unordered_set<VertexId> &verts = touched[&contig];
    if (fragment.isInnerFragment()) {
        verts.insert(fragment.vertex);
    } else if (fragment.isEdgeSegment()) {
        verts.insert(fragment.edge->getStart().getId());
        verts.insert(fragment.edge->getFinish().getId());
    }
}

void ag::DLLAlignmentStorage::addContig(Contig &contig, DLLAlignmentPath &&p) {
    DLLAlignmentPath &path = (alignments[contig.getInnerId()] = std::move(p));
    for (DLLPathPosition pos = path.begin(); pos != path.end(); ++pos) {
        if (pos->isEdgeSegment()) {
            edge_map[pos->edge].emplace_back(pos);
        } else if (pos->isInnerFragment()) {
            vertex_map[(*pos).vertex].emplace_back(pos);
        } else {
            VERIFY(false);
        }
    }
}

ag::DLLAlignmentStorage::DLLPathPosition ag::DLLAlignmentStorage::insertBefore(Edge &edge, DLLPathPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.edge== edge.getId());
    DLLPathPosition new_pos = pos.insertBefore(fragment);
    edge_map[edge.getId()].emplace_back(new_pos);
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPathPosition ag::DLLAlignmentStorage::insertAfter(Edge &edge, DLLPathPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.edge== edge.getId());
    DLLPathPosition new_pos = pos.insertAfter(fragment);
    edge_map[edge.getId()].emplace_back(new_pos);
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPathPosition ag::DLLAlignmentStorage::insertBefore(Vertex &vertex, DLLPathPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.vertex== vertex.getId());
    DLLPathPosition new_pos = pos.insertBefore(fragment);
    vertex_map[vertex.getId()].emplace_back(new_pos);
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPathPosition ag::DLLAlignmentStorage::insertAfter(Vertex &vertex, DLLPathPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.vertex== vertex.getId());
    DLLPathPosition new_pos = pos.insertAfter(fragment);
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
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    for (DLLPathPosition &pos : vertex_map.at(v.getId())) {
        VERIFY(pos.valid());
        VERIFY(pos->isInnerFragment());
        recordTouch(touched, *pos);
        pos.erase();
    }
    vertex_map.erase(v.getId());
    for (auto &p : touched)
        logPath(*p.first, "fireDeleteVertex");
}

void ag::DLLAlignmentStorage::fireDeleteEdge(Edge &e) {
    if (!hasRecords(e))
        return;
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    for (DLLPathPosition &pos : edge_map.at(e.getId())) {
        if (!pos.valid()) continue;
        recordTouch(touched, *pos);
        if (!pos.extracted()) {
            if (pos->edge->isSuffix()) {
                DLLPathPosition cur_pos = pos;
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
                DLLPathPosition cur_pos = pos;
                while (cur_pos->isLinked(*cur_pos.next()) && cur_pos->edge->isPrefix()) {
                    cur_pos = cur_pos.next();
                    cur_pos.prev().extract();
                }
                if (!cur_pos->isLinked(*cur_pos.next()) && cur_pos->edge->isPrefix()) {
                    insertAfter(cur_pos->edge->getFinish(), cur_pos,
                        AlignmentFragment::InnerFragment(cur_pos->seg, cur_pos->edge->getFinish(), cur_pos->cut_right));
                    cur_pos.extract();
                }
            }
        }
        pos.erase();
    }
    edge_map.erase(e.getId());
    for (auto &p : touched)
        logPath(*p.first, "fireDeleteEdge");
}

void ag::DLLAlignmentStorage::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    if (!hasRecords(e)) return;
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    for (DLLPathPosition &pos : edge_map.at(e.getId())) {
        recordTouch(touched, *pos);
        bool link_left = pos.prev()->isLinked(*pos);
        bool link_right = pos->isLinked(*pos.next());
        if (!link_left && !link_right) {
            insertBefore(v, pos, AlignmentFragment::InnerFragment(pos->seg, v, pos->cut_left));
        } else {
            if (link_left) {
                insertBefore(v.incFront(), pos, AlignmentFragment::EdgeSegment(pos->seg, v.incFront(), pos->cut_left));
            }
            if (link_right) {
                insertAfter(v.front(), pos, AlignmentFragment::EdgeSegment(pos->seg, v.front(), pos->cut_left));
            }
        }
        pos.extract();
    }
    for (auto &p : touched) {
        logPath(*p.first, "fireEdgeToSupreVertex");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireEdgeToSupreVertex");
    }
}

void ag::DLLAlignmentStorage::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    size_t shift = 0;
    for (Edge &edge: path.edges()) {
        if (hasRecords(edge)) {
            for (DLLPathPosition &pos : edge_map.at(edge.getId())) {
                if (!pos.valid() || pos.extracted())
                    continue;
                recordTouch(touched, *pos);
                size_t cut_left = shift + pos->cut_left;
                DLLPathPosition next_pos = pos.next();
                DLLPathPosition prev_pos = pos.prev();
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
                    AlignmentFragment new_al = AlignmentFragment::EdgeSegment(seg, new_vertex.incFront(), 0);
                    if (prev_pos->isLinked(new_al)) {
                        VERIFY(shift == 0);
                        insertAfter(new_vertex.incFront(), prev_pos, new_al);
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
        }
        shift += edge.rc().truncSize();
    }

    shift = 0;
    for (Vertex &v : path.innerVertices()) {
        shift += v.rc().front().truncSize();
        if (hasRecords(v))
            for (DLLPathPosition &pos : vertex_map.at(v.getId())) {
                if (!pos.valid())
                    continue;
                recordTouch(touched, *pos);
                insertBefore(new_vertex, pos, AlignmentFragment::InnerFragment(pos->seg, new_vertex, shift + pos->cut_left));
                pos.extract();
            }
    }
    for (auto &p : touched) {
        logPath(*p.first, "fireMergePath");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireMergePath");
    }
}

void ag::DLLAlignmentStorage::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    size_t shift = 0;
    for (Edge &edge: path.edges()) {
        if (hasRecords(edge)) {
            for (DLLPathPosition &pos : edge_map.at(edge.getId())) {
                if (!pos.valid() || pos.extracted())
                    continue;
                recordTouch(touched, *pos);
                size_t extra_shift = 0;
                AlignmentFragment new_fragment = AlignmentFragment::EdgeSegment(pos->seg, new_edge, shift + pos->cut_left);
                DLLPathPosition last_pos = pos.next();
                for (; last_pos.prev()->isLinked(*last_pos) && last_pos->edge->getStart() != path.getFinish(); last_pos = last_pos.next()) {
                    extra_shift += last_pos.prev()->edge->rc().truncSize();
                    AlignmentFragment edge_fragment = AlignmentFragment::EdgeSegment(last_pos->seg, new_edge, shift + extra_shift + last_pos->cut_left);
                    new_fragment = new_fragment.merge(edge_fragment);
                    last_pos.prev().extract();
                }
                last_pos.prev().extract();
                insertBefore(new_edge, last_pos, new_fragment);
            }
        }
        shift += edge.rc().truncSize();
    }
    shift = 0;
    Vertex &new_vertex = new_edge.isSuffix() ? new_edge.getStart() : new_edge.getFinish();
    for (Vertex &v : path.innerVertices()) {
        shift += v.rc().front().truncSize();
        if (hasRecords(v))
            for (DLLPathPosition &pos : vertex_map.at(v.getId())) {
                if (!pos.valid())
                    continue;
                recordTouch(touched, *pos);
                insertBefore(new_vertex, pos, AlignmentFragment::InnerFragment(pos->seg, new_vertex, shift + pos->cut_left));
                pos.extract();
            }
    }
    for (auto &p : touched) {
        logPath(*p.first, "fireMergePathToEdge");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireMergePathToEdge");
    }
}

void ag::DLLAlignmentStorage::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &left_al,
    const AlignmentForm &right_al) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
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
        for (DLLPathPosition &pos : edge_map.at(left.getId())) {
            if (!pos.valid()) continue;
            VERIFY(!pos.extracted());
            recordTouch(touched, *pos);
            if (left.truncSize() > pos->cut_left + min_right_cut) {
                size_t extra_cut = std::max(min_right_cut, pos->cut_right) - pos->cut_right;
                AlignmentFragment new_fragment = AlignmentFragment::EdgeSegment(pos->seg.shrinkRightBy(extra_cut), new_edge, pos->cut_left);
                insertBefore(new_edge, pos, new_fragment);
            }
            pos.extract();
        }
        for (DLLPathPosition &pos : edge_map.at(left.rc().getId())) {
            if (!pos.valid()) continue;
            VERIFY(!pos.extracted());
            recordTouch(touched, *pos);
            if (left.rc().truncSize() > pos->cut_right + min_right_cut) {
                size_t extra_cut = std::max(min_right_cut, pos->cut_left) - pos->cut_left;
                AlignmentFragment new_fragment = AlignmentFragment::EdgeSegment(pos->seg.shrinkLeftBy(extra_cut),
                    new_edge.rc(), new_edge.fullSize() - pos->cut_left - pos->seg.size() + extra_cut);
                insertAfter(new_edge.rc(), pos, new_fragment);
            }
            pos.extract();
        }
    }
    for (auto &p : touched) {
        logPath(*p.first, "fireMergeTipsToEdge");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireMergeTipsToEdge");
    }
}

void ag::DLLAlignmentStorage::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    if (hasRecords(edge)) {
        for (DLLPathPosition &pos : edge_map.at(edge.getId())) {
            recordTouch(touched, *pos);
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
    for (auto &p : touched) {
        logPath(*p.first, "fireSplitEdge");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireSplitEdge");
    }
}

void ag::DLLAlignmentStorage::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    for (Edge &inc: core.incoming()) {
        if (!hasRecords(inc))
            continue;
        for (DLLPathPosition &pos : edge_map.at(inc.getId())) {
            VERIFY(pos.next()->isEdgeSegment());
            VERIFY(pos.next()->edge->getStart() == core);
            DLLPathPosition left = pos;
            DLLPathPosition right = pos.next();
            VERIFY(left->isLinked(*right));
            if (resolution.contains(*left->edge, *right->edge)) {
                Vertex &new_vertex = resolution.get(*left->edge, *right->edge);
//                The interesting new structure is new_vertex itself, not just left/right's old edges.
                touched[&left->seg.contig()].insert(new_vertex.getId());
                bool is_linked = false;
                if (left.prev()->isLinked(*left)) {
                    AlignmentFragment al = AlignmentFragment::EdgeSegment(left->seg.unite(right->seg),
                                                                          new_vertex.incFront(), left->cut_left);
                    insertBefore(new_vertex.incFront(), left, al);
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
    for (auto &p : touched) {
        logPath(*p.first, "fireResolveVertex");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireResolveVertex");
    }
}

std::vector<std::string> ag::DLLAlignmentStorage::passingContigs(const Vertex &vertex) const {
    std::vector<std::string> res;
    if (!hasRecords(vertex))
        return res;
    std::unordered_set<std::string> seen;
    for (const DLLPathPosition &pos : vertex_map.at(vertex.getId())) {
        VERIFY(pos->isInnerFragment());
        const std::string &name = pos->seg.contig().getInnerId();
        if (seen.insert(name).second)
            res.push_back(name);
    }
    return res;
}

std::vector<std::string> ag::DLLAlignmentStorage::passingForwardContigs(const Vertex &vertex) const {
    std::vector<std::string> res;
    if (!hasRecords(vertex))
        return res;
    std::unordered_set<std::string> seen;
    for (const DLLPathPosition &pos : vertex_map.at(vertex.getId())) {
        VERIFY(pos->isInnerFragment());
        const std::string &name = pos->seg.contig().getInnerId();
        if (!startsWith(name, "-") && seen.insert(name).second)
            res.push_back(name);
    }
    return res;
}

std::function<std::string(const ag::Vertex &)> ag::DLLAlignmentStorage::getVertexTooltipper() const {
    std::function<std::string(const Vertex &)> res = [this](const Vertex &v)->std::string {
        if (!hasRecords(v)) return "";
        std::stringstream ss;
        for (DLLPathPosition pos : vertex_map.at(v.getId())) {
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
        for (DLLPathPosition pos : edge_map.at(e.getId())) {
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
        for (DLLPathPosition pos: p.second) {
            out << *pos << std::endl;
        }
    }
    out << "Edges" << std::endl;
    for (auto &p: edge_map) {
        out << p.first << std::endl;
        for (DLLPathPosition pos: p.second) {
            out << *pos << std::endl;
        }
    }
}
