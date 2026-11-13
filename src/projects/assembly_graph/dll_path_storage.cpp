#include "dll_path_storage.hpp"

#include "dll_path.hpp"
#include "sequences/contigs.hpp"
#include <unordered_set>
#include <set>

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

ag::DLLAlignmentStorage::PositionRecord &ag::DLLAlignmentStorage::getRecords(ConstVertexId vid) const {
    PositionRecord *res = nullptr;
    vertex_map.find_fn(vid, [&res](const std::unique_ptr<PositionRecord> &ptr) {res = ptr.get();});
    VERIFY(res != nullptr);
    return *res;
}

ag::DLLAlignmentStorage::PositionRecord &ag::DLLAlignmentStorage::getRecords(ConstEdgeId eid) const {
    PositionRecord *res = nullptr;
    edge_map.find_fn(eid, [&res](const std::unique_ptr<PositionRecord> &ptr) {res = ptr.get();});
    VERIFY(res != nullptr);
    return *res;
}

void ag::DLLAlignmentStorage::logPath(Contig &contig, const std::string &source) {
    if (!log_stream)
        return;
    std::lock_guard<std::mutex> guard(pathLock(contig));
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
    path_locks.emplace(&contig, std::make_unique<std::mutex>());
    DLLAlignmentPath &path = (alignments[contig.getInnerId()] = std::move(p));
    for (DLLPathPosition pos = path.begin(); pos != path.end(); ++pos) {
        if (pos->isEdgeSegment()) {
            edge_map.upsert(pos->edge, [&](std::unique_ptr<PositionRecord> &rec) {rec->emplace_back(pos);},
            std::make_unique<PositionRecord>(1, pos));
        } else if (pos->isInnerFragment()) {
            vertex_map.upsert((*pos).vertex, [&](std::unique_ptr<PositionRecord> &rec) {rec->emplace_back(pos);},
            std::make_unique<PositionRecord>(1, pos));
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
    edge_map.upsert(edge.getId(), [&](std::unique_ptr<PositionRecord> &rec) {rec->emplace_back(new_pos);},
            std::make_unique<PositionRecord>(1, new_pos));
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPathPosition ag::DLLAlignmentStorage::insertAfter(Edge &edge, DLLPathPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.edge== edge.getId());
    DLLPathPosition new_pos = pos.insertAfter(fragment);
    edge_map.upsert(edge.getId(), [&](std::unique_ptr<PositionRecord> &rec) {rec->emplace_back(new_pos);},
            std::make_unique<PositionRecord>(1, new_pos));
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPathPosition ag::DLLAlignmentStorage::insertBefore(Vertex &vertex, DLLPathPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.vertex== vertex.getId());
    DLLPathPosition new_pos = pos.insertBefore(fragment);
    vertex_map.upsert(vertex.getId(), [&](std::unique_ptr<PositionRecord> &rec) {rec->emplace_back(new_pos);},
            std::make_unique<PositionRecord>(1, new_pos));
    return new_pos;
}

ag::DLLAlignmentStorage::DLLPathPosition ag::DLLAlignmentStorage::insertAfter(Vertex &vertex, DLLPathPosition pos,
    AlignmentFragment fragment) {
    VERIFY(fragment.checkSeqMatch());
    VERIFY(fragment.vertex== vertex.getId());
    DLLPathPosition new_pos = pos.insertAfter(fragment);
    vertex_map.upsert(vertex.getId(), [&](std::unique_ptr<PositionRecord> &rec) {rec->emplace_back(new_pos);},
            std::make_unique<PositionRecord>(1, new_pos));
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
    for (DLLPathPosition &pos : getRecords(v.getId())) {
        VERIFY(pos.valid());
        VERIFY(pos->isInnerFragment());
        recordTouch(touched, *pos);
        pos.erase();
    }
    vertex_map.erase(v.getId());
    for (auto &p : touched)
        logPath(*p.first, "fireDeleteVertex");
}

void ag::DLLAlignmentStorage::keepStart(DLLPathPosition pos) {
    Vertex &start = pos->edge->getStart();
    Segment<Contig> start_seg(pos->seg.contig(), pos->seg.left, pos->seg.left + start.size());
    DLLPathPosition prev = pos.prev();
    bool covered = prev->isLinked(*pos) || (prev->isInnerFragment() && prev->vertex == start.getId() && prev->seg == start_seg);
    if (pos->cut_left == 0 && !covered)
        insertBefore(start, pos, AlignmentFragment::InnerFragment(start_seg, start, 0));
}

void ag::DLLAlignmentStorage::keepFinish(DLLPathPosition pos) {
    Vertex &finish = pos->edge->getFinish();
    Segment<Contig> finish_seg(pos->seg.contig(), pos->seg.right - finish.size(), pos->seg.right);
    DLLPathPosition next = pos.next();
    bool covered = pos->isLinked(*next) || (next->isInnerFragment() && next->vertex == finish.getId() && next->seg == finish_seg);
    if (pos->cut_right == 0 && !covered)
        insertAfter(finish, pos, AlignmentFragment::InnerFragment(finish_seg, finish, 0));
}

void ag::DLLAlignmentStorage::fireDeleteEdge(Edge &e) {
    if (!hasRecords(e))
        return;
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    for (DLLPathPosition &pos : getRecords(e.getId())) {
        if (!pos.valid()) continue;
//        Neighbouring fragments may be changed by threads deleting other edges of the same path.
        std::lock_guard<std::mutex> guard(pathLock(pos->seg.contig()));
        recordTouch(touched, *pos);
        if (!pos.extracted()) {
//            A suffix (prefix) edge spells its start (finish), so the fragment becomes an inner fragment there unless
//            the linked neighbour already ends (starts) at it. Only this fragment goes: the vertices of a linked run
//            of suffix (prefix) edges before (after) it are live proper infixes of the contig and stay in the path.
            if (pos->edge->isSuffix()) {
                keepFinish(pos);
                if (!pos.prev()->isLinked(*pos))
                    insertBefore(pos->edge->getStart(), pos,
                        AlignmentFragment::InnerFragment(pos->seg, pos->edge->getStart(), pos->cut_left));
            } else if (pos->edge->isPrefix()) {
                keepStart(pos);
                if (!pos->isLinked(*pos.next()))
                    insertAfter(pos->edge->getFinish(), pos,
                        AlignmentFragment::InnerFragment(pos->seg, pos->edge->getFinish(), pos->cut_left));
            } else {
                keepStart(pos);
                keepFinish(pos);
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
    for (DLLPathPosition &pos : getRecords(e.getId())) {
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
            for (DLLPathPosition &pos : getRecords(edge.getId())) {
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
//                The path's start and finish survive the merge: a fragment reaching one of them uncut stays an edge
//                segment on the new vertex's edge to it (both edges spell new_vertex), linked to its neighbour or not,
//                so that vertex stays in the path.
                bool on_edge = false;
                if (cut_left == 0) {
                    VERIFY(shift == 0);
                    insertAfter(new_vertex.incFront(), prev_pos, AlignmentFragment::EdgeSegment(seg, new_vertex.incFront(), 0));
                    on_edge = true;
                }
                if (cut_right == 0) {
                    insertBefore(new_vertex.front(), next_pos, AlignmentFragment::EdgeSegment(seg, new_vertex.front(), cut_left));
                    on_edge = true;
                }
                if (!on_edge) {
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
            for (DLLPathPosition &pos : getRecords(v.getId())) {
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
            for (DLLPathPosition &pos : getRecords(edge.getId())) {
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
//    An inner fragment on an interior vertex (at shift in new_edge) goes to new_edge itself. If new_edge is a prefix
//    (suffix) edge it spells its finish (start), so the fragment goes into that vertex at the same position.
    shift = 0;
    for (Vertex &v : path.innerVertices()) {
        shift += v.rc().front().truncSize();
        if (hasRecords(v))
            for (DLLPathPosition &pos : getRecords(v.getId()))
                if (pos.valid() && !pos.extracted()) {
                    recordTouch(touched, *pos);
                    if (new_edge.isPrefix() || new_edge.isSuffix()) {
                        Vertex &spelled = new_edge.isPrefix() ? new_edge.getFinish() : new_edge.getStart();
                        insertBefore(spelled, pos, AlignmentFragment::InnerFragment(pos->seg, spelled, shift + pos->cut_left));
                    } else {
                        insertBefore(new_edge, pos, AlignmentFragment::EdgeSegment(pos->seg, new_edge, shift + pos->cut_left));
                    }
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
        for (DLLPathPosition &pos : getRecords(left.getId())) {
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
        for (DLLPathPosition &pos : getRecords(left.rc().getId())) {
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
        for (DLLPathPosition &pos : getRecords(edge.getId())) {
            recordTouch(touched, *pos);
            size_t e_left_cut = 0;
            size_t e_right_cut = edge.truncSize();
            for (Edge &e : split.edges()) {
                e_right_cut -= e.truncSize();
                if (pos->cut_left < edge.fullSize() - e_right_cut - e.getFinish().size() &&
                    pos->cut_right < edge.fullSize() - e_left_cut - e.getStart().size()) {
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
    std::unordered_map<EdgeId, std::unordered_map<EdgeId, EdgeId>> direct;
    for (const std::pair<EdgeId, InOutEdgePair> &p : resolution.newEdges())
        direct[p.second.incomingId()][p.second.outgoingId()] = p.first;
//    Only a linked pair in -> out of a resolved connection is rerouted. Everything else stays on the old edges and is
//    handled by fireDeleteEdge: the contig keeps u and w but loses core.
    for (Edge &inc: core.incoming())
        if (hasRecords(inc))
            for (DLLPathPosition &left : getRecords(inc.getId())) {
                DLLPathPosition right = left.next();
                bool paired = left.valid() && !left.extracted() && left->isLinked(*right) && right->edge->getStart() == core;
                if (paired) {
                    Edge &out = *right->edge;
                    Contig &contig = left->seg.contig();
//                    The contig covers [L, R); core occurs at [c0, c1).
                    size_t L = left->seg.left, R = right->seg.right;
                    size_t c0 = right->seg.left, c1 = left->seg.right;
                    Segment<Contig> whole = left->seg.unite(right->seg);
                    auto dit = direct.find(inc.getId());
                    bool rerouted = true;
                    if (dit != direct.end() && dit->second.find(out.getId()) != dit->second.end()) {
//                        Direct edge u -> w spells u + out.truncSeq().
                        Edge &edge = *dit->second.at(out.getId());
                        touched[&contig].insert(edge.getStart().getId());
                        insertBefore(edge, left, AlignmentFragment::EdgeSegment(whole, edge, left->cut_left));
                    } else if (resolution.contains(inc, out)) {
//                        newv = a + core + b at [c0 - 1, c1 + 1): u -> newv spells u + b, newv -> w spells a + w.
                        Vertex &newv = resolution.get(inc, out);
                        touched[&contig].insert(newv.getId());
                        size_t nl = std::max(L, c0 - 1), nr = std::min(R, c1 + 1);
//                        A side gets its edge if the contig reaches beyond newv there, or covers u from its start
//                        (w to its end): that keeps u (w) in the path and linked to its neighbour.
                        bool to_u = L + 1 < c0 || left->cut_left == 0;
                        bool to_w = R > c1 + 1 || right->cut_right == 0;
                        if (to_u)
                            insertBefore(newv.incFront(), left, AlignmentFragment::EdgeSegment(
                                    Segment<Contig>(contig, L, nr), newv.incFront(), left->cut_left));
                        if (to_w)
                            insertBefore(newv.front(), left, AlignmentFragment::EdgeSegment(
                                    Segment<Contig>(contig, nl, R), newv.front(), nl - (c0 - 1)));
                        if (!to_u && !to_w)
                            insertBefore(newv, left, AlignmentFragment::InnerFragment(
                                    Segment<Contig>(contig, nl, nr), newv, nl - (c0 - 1)));
                    } else {
                        rerouted = false;
                    }
                    if (rerouted) {
                        right.extract();
                        left.extract();
                    }
                }
            }
    for (auto &p : touched) {
        logPath(*p.first, "fireResolveVertex");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireResolveVertex");
    }
}

void ag::DLLAlignmentStorage::fireFlipVertex(Vertex &core, const std::unordered_map<EdgeId, EdgeId> &extensions) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
//    Forward: u -inc-> core -e-> w becomes u -> w, spelling u + e.truncSeq(); core is a suffix of u.
//    Mirrored: w' -e-> core -inc-> u' becomes w' -> u'; core is a prefix of u'.
    const bool forward = core.outDeg() != 1;
    Edge &inc = forward ? core.incFront() : core.front();
    Vertex &u = forward ? inc.getStart() : inc.getFinish();
    const size_t shift = forward ? u.size() - core.size() : 0;
    std::vector<EdgeId> replaced;
    if (forward)
        for (Edge &e : core) replaced.push_back(e.getId());
    else
        for (Edge &e : core.incoming()) replaced.push_back(e.getId());
    for (EdgeId eid : replaced)
        if (hasRecords(*eid)) {
            Edge &new_edge = *extensions.at(eid);
            for (DLLPathPosition &pos : getRecords(eid))
                if (pos.valid() && !pos.extracted()) {
                    recordTouch(touched, *pos);
//                    A linked fragment on inc merges into the new fragment, which then stays linked to whatever was
//                    linked to it.
                    DLLPathPosition other = forward ? pos.prev() : pos.next();
                    bool merge = forward ? other->isLinked(*pos) && other->edge == inc.getId()
                                         : pos->isLinked(*other) && other->edge == inc.getId();
                    AlignmentFragment new_fragment = merge ?
                            AlignmentFragment::EdgeSegment(pos->seg.unite(other->seg), new_edge, forward ? other->cut_left : pos->cut_left) :
                            AlignmentFragment::EdgeSegment(pos->seg, new_edge, shift + pos->cut_left);
                    if (merge)
                        other.extract();
                    insertBefore(new_edge, pos, new_fragment);
                    pos.extract();
                }
        }
    if (hasRecords(core))
        for (DLLPathPosition &pos : getRecords(core.getId()))
            if (pos.valid() && !pos.extracted()) {
                recordTouch(touched, *pos);
                insertBefore(u, pos, AlignmentFragment::InnerFragment(pos->seg, u, shift + pos->cut_left));
                pos.extract();
            }
    for (auto &p : touched) {
        logPath(*p.first, "fireFlipVertex");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireFlipVertex");
    }
}

std::vector<std::string> ag::DLLAlignmentStorage::passingContigs(const Vertex &vertex) const {
    std::vector<std::string> res;
    if (!hasRecords(vertex))
        return res;
    std::unordered_set<std::string> seen;
    for (const DLLPathPosition &pos : getRecords(vertex.getId())) {
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
    for (const DLLPathPosition &pos : getRecords(vertex.getId())) {
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
        for (DLLPathPosition pos : getRecords(v.getId())) {
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
        for (DLLPathPosition pos : getRecords(e.getId())) {
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
    for (auto &p: vertex_map.lock_table()) {
        out << p.first << std::endl;
        for (DLLPathPosition pos: *p.second) {
            out << *pos << std::endl;
        }
    }
    out << "Edges" << std::endl;
    for (auto &p: edge_map.lock_table()) {
        out << p.first << std::endl;
        for (DLLPathPosition pos: *p.second) {
            out << *pos << std::endl;
        }
    }
}

namespace {
//    g (a fragment of the rc contig's path) is the mirror of f (a fragment of the contig's path).
    bool IsMirror(const ag::AlignmentFragment &f, const ag::AlignmentFragment &g, size_t contig_size) {
        bool same_object = f.isEdgeSegment() ? g.isEdgeSegment() && g.edge == f.edge->rc().getId()
                                             : g.isInnerFragment() && g.vertex == f.vertex->rc().getId();
        return same_object && g.cut_left == f.cut_right && g.cut_right == f.cut_left &&
               g.seg.left == contig_size - f.seg.right && g.seg.right == contig_size - f.seg.left;
    }
}

bool ag::DLLAlignmentStorage::checkConsistency(std::ostream &os) {
    AssemblyGraph &graph = getGraph();
    bool ok = true;
    std::unordered_set<ConstEdgeId> live_edges;
    std::unordered_set<ConstVertexId> live_vertices;
    for (Edge &e : graph.edges())
        live_edges.insert(e.getId());
    size_t min_size = size_t(-1);
    for (Vertex &v : graph.vertices()) {
        live_vertices.insert(v.getId());
        if (!v.getSeq().empty())
            min_size = std::min(min_size, v.size());
    }
//    1. Index side: every live entry is a fragment on its key; count them to compare with the paths.
    size_t indexed = 0;
    for (auto &p : edge_map.lock_table())
        for (const DLLPathPosition &pos : *p.second)
            if (pos.valid() && !pos.extracted()) {
                indexed++;
                if (!pos->isEdgeSegment() || pos->edge != p.first) {
                    os << "DLLAlignmentStorage: entry of edge " << p.first << " points to fragment " << *pos << std::endl;
                    ok = false;
                }
            }
    for (auto &p : vertex_map.lock_table())
        for (const DLLPathPosition &pos : *p.second)
            if (pos.valid() && !pos.extracted()) {
                indexed++;
                if (!pos->isInnerFragment() || pos->vertex != p.first) {
                    os << "DLLAlignmentStorage: entry of vertex " << p.first << " points to fragment " << *pos << std::endl;
                    ok = false;
                }
            }
//    Vertex index for 4: each vertex under its first min_size nucleotides.
    std::unordered_map<std::string, std::vector<std::pair<Vertex *, std::string>>> by_prefix;
    if (min_size != size_t(-1))
        for (Vertex &v : graph.vertices())
            if (!v.getSeq().empty()) {
                std::string s = v.getSeq().str();
                by_prefix[s.substr(0, min_size)].emplace_back(&v, std::move(s));
            }
    size_t fragments = 0;
    for (Contig &contig : contigs) {
        DLLAlignmentPath &path = alignments.at(contig.getInnerId());
        std::set<std::pair<const Vertex *, size_t>> present;
        const AlignmentFragment *prev = nullptr;
        for (const AlignmentFragment &f : path) {
            fragments++;
//            1. Path side.
            bool live = f.isEdgeSegment() ? live_edges.count(f.edge) > 0 : f.isInnerFragment() && live_vertices.count(f.vertex) > 0;
            if (!f.valid() || !live || &f.seg.contig() != &contig) {
                os << "DLLAlignmentStorage: contig " << contig.getInnerId() << " has a dead or foreign fragment" << std::endl;
                ok = false;
            } else {
                if (!f.checkSeqMatch()) {
                    os << "DLLAlignmentStorage: fragment " << f << " does not spell its contig segment" << std::endl;
                    ok = false;
                }
//                2.
                if (prev != nullptr && (f.seg.left < prev->seg.left || f.seg.right < prev->seg.right)) {
                    os << "DLLAlignmentStorage: fragments " << *prev << " and " << f << " are out of order" << std::endl;
                    ok = false;
                }
                if (f.isInnerFragment()) {
                    if (f.cut_left <= f.seg.left)
                        present.emplace(&*f.vertex, f.seg.left - f.cut_left);
                } else {
                    if (f.cut_left == 0)
                        present.emplace(&f.edge->getStart(), f.seg.left);
                    if (f.cut_right == 0)
                        present.emplace(&f.edge->getFinish(), f.seg.right - f.edge->getFinish().size());
                }
            }
            prev = &f;
        }
//        4.
        if (min_size != size_t(-1)) {
            std::string s = contig.getSeq().str();
            for (size_t p = 1; p + min_size < s.size(); p++) {
                auto cand = by_prefix.find(s.substr(p, min_size));
                if (cand != by_prefix.end())
                    for (const std::pair<Vertex *, std::string> &v : cand->second)
                        if (p + v.second.size() < s.size() && s.compare(p, v.second.size(), v.second) == 0 &&
                                present.count({v.first, p}) == 0) {
                            os << "DLLAlignmentStorage: vertex " << v.first->getId() << " occurs inside contig "
                               << contig.getInnerId() << " at " << p << " but is not in its path" << std::endl;
                            ok = false;
                        }
            }
        }
    }
//    3. addContig stores contigs in pairs: the contig, then its rc.
    for (auto it = contigs.begin(); it != contigs.end(); std::advance(it, 2)) {
        Contig &contig = *it;
        Contig &rc_contig = *std::next(it);
        std::vector<const AlignmentFragment *> fw, bw;
        for (const AlignmentFragment &f : alignments.at(contig.getInnerId())) fw.push_back(&f);
        for (const AlignmentFragment &f : alignments.at(rc_contig.getInnerId())) bw.push_back(&f);
        bool mirrored = fw.size() == bw.size();
        for (size_t i = 0; mirrored && i < fw.size(); i++)
            mirrored = IsMirror(*fw[i], *bw[bw.size() - 1 - i], contig.fullSize());
        if (!mirrored) {
            os << "DLLAlignmentStorage: path of " << rc_contig.getInnerId() << " does not mirror the path of "
               << contig.getInnerId() << std::endl;
            ok = false;
        }
    }
    if (fragments != indexed) {
        os << "DLLAlignmentStorage: " << fragments << " fragments in paths but " << indexed << " live index entries" << std::endl;
        ok = false;
    }
    return ok;
}
