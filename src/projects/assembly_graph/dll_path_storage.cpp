#include "dll_path_storage.hpp"

#include "dll_path.hpp"
#include "sequences/contigs.hpp"
#include <omp.h>
#include <sstream>
#include <unordered_set>
#include <set>

namespace ag {
    class Edge;
}

void ag::PathDrawer::draw(Contig &contig, const std::vector<VertexId> &vertices, const std::string &event_tag) {
    if (vertices.empty() || omp_in_parallel())
        return;
    const std::string &contig_name = contig.getInnerId();
    if (startsWith(contig_name, "-") || contig.truncSize() > max_contig_length)
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

void ag::DLLAlignmentStorage::logChange(Contig &contig, const std::unordered_set<VertexId> &touched,
                                        const std::string &source) {
    if (!log_stream)
        return;
    std::lock_guard<std::mutex> guard(pathLock(contig));
    DLLAlignmentPath &path = alignments.at(contig.getInnerId());
    std::vector<DLLPathPosition> positions;
    auto collect = [&](const PositionRecord &record) {
        for (const DLLPathPosition &pos : record) {
            if (!pos.valid() || pos.extracted() || &pos->seg.contig() != &contig)
                continue;
            positions.push_back(pos);
            if (pos != path.begin())
                positions.push_back(pos.prev());
            if (pos.next() != path.end())
                positions.push_back(pos.next());
        }
    };
    for (VertexId v : touched) {
        if (hasRecords(*v))
            collect(getRecords(v));
        for (Edge &e : *v)
            if (hasRecords(e))
                collect(getRecords(e.getId()));
        for (Edge &e : v->incoming())
            if (hasRecords(e))
                collect(getRecords(e.getId()));
    }
    std::sort(positions.begin(), positions.end(), [](const DLLPathPosition &a, const DLLPathPosition &b) {
        return std::make_tuple(a->seg.left, a->seg.right, &*a) < std::make_tuple(b->seg.left, b->seg.right, &*b);
    });
    positions.erase(std::unique(positions.begin(), positions.end(),
                                [](const DLLPathPosition &a, const DLLPathPosition &b) {return &*a == &*b;}),
                    positions.end());
    std::ostringstream message;
    message << "Change in observed path during " << source << ": " << contig.getInnerId() << ":";
    for (const DLLPathPosition &pos : positions)
        message << " " << *pos;
    message << "\n";
    std::lock_guard<std::mutex> log_guard(*log_lock);
    *log_stream << message.str() << std::flush;
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
        logChange(*p.first, p.second, "fireDeleteVertex");
}

//    A fragment that starts inside the edge's start (cut_left > 0) does not contain it, and its start segment may not even
//    fit in the contig.
void ag::DLLAlignmentStorage::keepStart(DLLPathPosition pos) {
    if (pos->cut_left != 0)
        return;
    Vertex &start = pos->edge->getStart();
    Segment<Contig> start_seg(pos->seg.contig(), pos->seg.left, pos->seg.left + start.size());
    DLLPathPosition prev = pos.prev();
    bool covered = prev->isLinked(*pos) || (prev->isInnerFragment() && prev->vertex == start.getId() && prev->seg == start_seg);
    if (!covered)
        insertBefore(start, pos, AlignmentFragment::InnerFragment(start_seg, start, 0));
}

void ag::DLLAlignmentStorage::keepFinish(DLLPathPosition pos) {
    if (pos->cut_right != 0)
        return;
    Vertex &finish = pos->edge->getFinish();
    Segment<Contig> finish_seg(pos->seg.contig(), pos->seg.right - finish.size(), pos->seg.right);
    DLLPathPosition next = pos.next();
    bool covered = pos->isLinked(*next) || (next->isInnerFragment() && next->vertex == finish.getId() && next->seg == finish_seg);
    if (!covered)
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
        logChange(*p.first, p.second, "fireDeleteEdge");
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
        logChange(*p.first, p.second, "fireEdgeToSupreVertex");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireEdgeToSupreVertex");
    }
}

void ag::DLLAlignmentStorage::fireReplacePath(const PathReplacement &rep) {
    if (rep.newPath().size() == 0)
        mergeIntoDominant(rep);
    else
        relocate(rep);
}

//    Q = [D]. A run of linked fragments over the edges into P's start, P's edges and the edges out of P's finish is
//    rewritten run-wise, so that it stays linked: P's part joins the fragment on the edge entering D before it (that
//    edge spells all of D), else the fragment on the edge leaving D after it (extended back to where the run entered D),
//    else becomes an inner fragment of D. Positions in D are positions in P (D spells P); an edge out of P's finish
//    starts at the finish's position in D.
void ag::DLLAlignmentStorage::mergeIntoDominant(const PathReplacement &rep) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    const RAGraphPath &old_path = rep.oldPath();
    Vertex &dominant = rep.getDominant();
    Vertex &start = old_path.getStart();
    Vertex &finish = old_path.getFinish();
    const size_t finish_position = rep.oldPosition(old_path.size());
    std::unordered_map<ConstEdgeId, size_t> index;
    for (size_t i = 0; i < old_path.size(); i++)
        index[old_path.getEdge(i).getId()] = i;
    auto on_path = [&](const Edge &e) {return index.find(e.getId()) != index.end();};
    auto is_outer = [&](const Edge &e) {return !on_path(e) && (e.getFinish() == start || e.getStart() == finish);};
    auto involved = [&](const AlignmentFragment &f) {return f.isEdgeSegment() && (on_path(*f.edge) || is_outer(*f.edge));};
    auto target = [&](Edge &e) -> Edge & {
        return (e.getFinish() == start && rep.startReplaced()) || (e.getStart() == finish && rep.endReplaced()) ?
               rep.replacement(e) : e;
    };
    std::vector<DLLPathPosition> candidates;
    auto collect = [&](const Edge &e) {
        if (hasRecords(e))
            for (DLLPathPosition &pos : getRecords(e.getId()))
                candidates.push_back(pos);
    };
    for (Edge &e : old_path.edges())
        collect(e);
    for (Edge &e : start.incoming())
        collect(e);
    for (Edge &e : finish)
        collect(e);
    for (DLLPathPosition pos : candidates) {
        if (!pos.valid() || pos.extracted() || !involved(*pos))
            continue;
//        Neighbouring fragments of the contig may be changed by threads merging other paths.
        std::lock_guard<std::mutex> guard(pathLock(pos->seg.contig()));
        DLLPathPosition head = pos;
        while (head.prev()->isLinked(*head) && involved(*head.prev()))
            head = head.prev();
        std::vector<DLLPathPosition> run = {head};
        while (run.back()->isLinked(*run.back().next()) && involved(*run.back().next()))
            run.push_back(run.back().next());
        std::vector<AlignmentFragment> result;
//        A visit of D: absorbed into result.back() if the run entered D along an edge, else collected in visit_seg from
//        visit_from in D.
        bool in_dominant = false;
        bool absorbed = false;
        size_t visit_from = 0;
        Segment<Contig> visit_seg;
        for (DLLPathPosition &p : run) {
            recordTouch(touched, *p);
            const AlignmentFragment &f = *p;
            Edge &e = *f.edge;
            if (on_path(e)) {
                if (!in_dominant) {
                    in_dominant = true;
                    absorbed = !result.empty();
                    visit_from = rep.oldPosition(index[e.getId()]) + f.cut_left;
                    visit_seg = f.seg;
                }
                if (absorbed)
                    result.back() = AlignmentFragment::EdgeSegment(result.back().seg.unite(f.seg), *result.back().edge,
                                                                   result.back().cut_left);
                else
                    visit_seg = visit_seg.unite(f.seg);
            } else if (in_dominant) {
//                Leaving D along an edge out of the finish, which now starts at D's start.
                in_dominant = false;
                if (absorbed)
                    result.push_back(AlignmentFragment::EdgeSegment(
                            Segment<Contig>(f.seg.contig(), f.seg.left - finish_position, f.seg.right), target(e), 0));
                else
                    result.push_back(AlignmentFragment::EdgeSegment(visit_seg.unite(f.seg), target(e), visit_from));
            } else {
                const size_t shift = e.getStart() == finish ? finish_position : 0;
                result.push_back(AlignmentFragment::EdgeSegment(f.seg, target(e), shift + f.cut_left));
            }
        }
        if (in_dominant && !absorbed)
            result.push_back(AlignmentFragment::InnerFragment(visit_seg, dominant, visit_from));
        for (const AlignmentFragment &f : result) {
            if (f.isEdgeSegment())
                insertBefore(*f.edge, head, f);
            else
                insertBefore(dominant, head, f);
        }
        for (DLLPathPosition &p : run)
            p.extract();
    }
    for (size_t i = 0; i <= old_path.size(); i++) {
        Vertex &v = old_path.getVertex(i);
        if (!rep.survives(i) && hasRecords(v))
            for (DLLPathPosition &pos : getRecords(v.getId()))
                if (pos.valid() && !pos.extracted()) {
                    std::lock_guard<std::mutex> guard(pathLock(pos->seg.contig()));
                    recordTouch(touched, *pos);
                    insertBefore(dominant, pos, AlignmentFragment::InnerFragment(pos->seg, dominant, rep.oldPosition(i) + pos->cut_left));
                    pos.extract();
                }
    }
    for (auto &p : touched) {
        logChange(*p.first, p.second, "fireReplacePath");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireReplacePath");
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
        logChange(*p.first, p.second, "fireMergePathToEdge");
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
        logChange(*p.first, p.second, "fireMergeTipsToEdge");
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
        logChange(*p.first, p.second, "fireSplitEdge");
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
//                        newv = a + core + b at [c0 - 1, c1 + 1): u -> newv spells u + b, newv -> w spells a + w. A copy
//                        of core (d = 0) is at [c0, c1).
                        Vertex &newv = resolution.get(inc, out);
                        const size_t d = newv.size() == core.size() ? 0 : 1;
                        touched[&contig].insert(newv.getId());
                        size_t nl = std::max(L, c0 - d), nr = std::min(R, c1 + d);
//                        A side gets its edge if the contig reaches beyond newv there, or covers u from its start
//                        (w to its end): that keeps u (w) in the path and linked to its neighbour.
                        bool to_u = L + d < c0 || left->cut_left == 0;
                        bool to_w = R > c1 + d || right->cut_right == 0;
                        if (to_u)
                            insertBefore(newv.incFront(), left, AlignmentFragment::EdgeSegment(
                                    Segment<Contig>(contig, L, nr), newv.incFront(), left->cut_left));
                        if (to_w)
                            insertBefore(newv.front(), left, AlignmentFragment::EdgeSegment(
                                    Segment<Contig>(contig, nl, R), newv.front(), nl - (c0 - d)));
                        if (!to_u && !to_w)
                            insertBefore(newv, left, AlignmentFragment::InnerFragment(
                                    Segment<Contig>(contig, nl, nr), newv, nl - (c0 - d)));
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
        logChange(*p.first, p.second, "fireResolveVertex");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireResolveVertex");
    }
}

//    Fragments move by positions: a removed vertex's into D; an edge of P's onto the Q edge leaving its surviving start
//    or entering its surviving finish (into D if neither end survives); an edge around a replaced end onto its
//    replacement. A fragment's position is measured from its edge's (vertex's) start, so it shifts by the difference of
//    the old and the new start positions.
void ag::DLLAlignmentStorage::relocate(const PathReplacement &rep) {
    std::unordered_map<Contig *, std::unordered_set<VertexId>> touched;
    const RAGraphPath &old_path = rep.oldPath();
    const RAGraphPath &new_path = rep.newPath();
    Vertex &dominant = rep.getDominant();
    auto new_position = [&](const Vertex &v) {return v == dominant ? rep.dominantPosition() : 0;};
    auto move_edge = [&](Edge &from, Edge &to, size_t shift) {
        if (hasRecords(from))
            for (DLLPathPosition &pos : getRecords(from.getId()))
                if (pos.valid() && !pos.extracted()) {
                    recordTouch(touched, *pos);
                    insertBefore(to, pos, AlignmentFragment::EdgeSegment(pos->seg, to, shift + pos->cut_left));
                    pos.extract();
                }
    };
//    id is a vertex or an edge id.
    auto move_into_dominant = [&](auto id, size_t shift) {
        if (hasRecords(*id))
            for (DLLPathPosition &pos : getRecords(id))
                if (pos.valid() && !pos.extracted()) {
                    recordTouch(touched, *pos);
                    insertBefore(dominant, pos, AlignmentFragment::InnerFragment(pos->seg, dominant, shift + pos->cut_left));
                    pos.extract();
                }
    };
    for (size_t i = 0; i <= old_path.size(); i++)
        if (!rep.survives(i))
            move_into_dominant(ConstVertexId(old_path.getVertex(i).getId()), rep.oldPosition(i) - rep.dominantPosition());
    for (size_t i = 0; i < old_path.size(); i++) {
        Edge &e = old_path.getEdge(i);
        EdgeId target;
        for (Edge &q : new_path.edges())
            if ((rep.survives(i) && q.getStart() == e.getStart()) ||
                    (!rep.survives(i) && rep.survives(i + 1) && q.getFinish() == e.getFinish()))
                target = q.getId();
        if (target.valid())
            move_edge(e, *target, rep.oldPosition(i) - new_position(target->getStart()));
        else if (!rep.survives(i))
            move_into_dominant(ConstEdgeId(e.getId()), rep.oldPosition(i) - rep.dominantPosition());
    }
    if (rep.startReplaced())
        for (Edge &in : old_path.getStart().incoming())
            move_edge(in, rep.replacement(in), 0);
    if (rep.endReplaced())
        for (Edge &out : old_path.getFinish())
            move_edge(out, rep.replacement(out),
                      rep.oldPosition(old_path.size()) - new_position(new_path.getFinish()));
    for (auto &p : touched) {
        logChange(*p.first, p.second, "fireReplacePath");
        if (drawer)
            drawer->draw(*p.first, {p.second.begin(), p.second.end()}, "fireReplacePath");
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
//    Besides vertices with inner fragments, vertices the path passes through between edge fragments: an edge fragment
//    without a cut on a side spells the whole vertex on that side.
    std::function<std::string(const Vertex &)> res = [this](const Vertex &v)->std::string {
        if (hasRecords(v)) return "blue";
        for (const Edge &e : v.incoming())
            if (hasRecords(e))
                for (const DLLPathPosition &pos : getRecords(e.getId()))
                    if (pos.valid() && !pos.extracted() && pos->cut_right == 0)
                        return "blue";
        for (const Edge &e : v)
            if (hasRecords(e))
                for (const DLLPathPosition &pos : getRecords(e.getId()))
                    if (pos.valid() && !pos.extracted() && pos->cut_left == 0)
                        return "blue";
        return "";
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
//    Vertex index for 4: each vertex with sample segments under its first min_size nucleotides.
    std::unordered_map<std::string, std::vector<std::pair<Vertex *, std::string>>> by_prefix;
    if (min_size != size_t(-1))
        for (Vertex &v : graph.vertices())
            if (!v.getSeq().empty() && !v.coverage_info.empty()) {
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
