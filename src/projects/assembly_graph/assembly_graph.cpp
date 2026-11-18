#include "assembly_graph.hpp"

using namespace ag;

const hashing::htype VertexData::default_hash = 1651253415;

VertexData VertexData::SPGData(bool cyclic, bool inf_left, bool inf_right) {
    VertexData res;
    res.cyclic = cyclic;
    res.inf_left = inf_left;
    res.inf_right = inf_right;
    return res;
}

VertexData VertexData::DBGData(hashing::htype hash) {
    VertexData res;
    res.hash = hash;
    return res;
}

Edge &AssemblyGraph::addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeIdType eid, EdgeIdType rcid) {
    return addEdgeLockFree(start, end, full_sequence, EdgeData(), eid, rcid);
}

Edge &
AssemblyGraph::addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeIdType eid,
                               EdgeIdType rcid) {
    return addEdgeLockFree(start, end, tseq, rctseq, EdgeData(), eid, rcid);
}

Edge &AssemblyGraph::addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeIdType eid, EdgeIdType rcid) {
    return addEdge(start, end, full_seq, EdgeData(), eid, rcid);
}

//    Only the ends of the path are locked, although the edge lists of their outside neighbours change too: merging is
//    not run in parallel for now.
Vertex &ag::AssemblyGraph::mergePath(const GraphPath &path, Vertex::id_type id) {
    Locker<VertexId> locker = Locker<VertexId>::FromVector({path.getStart().getId(), path.getStart().rc().getId(), path.getRCStart().getId(), path.getFinish().rc().getId()});
    if (path.calculateSize() == 0)
        return path.getStart();
    return mergePathLockFree(path, id);
}

Vertex &ag::AssemblyGraph::flipVertex(Edge &inc) {
    Vertex &core = inc.getFinish();
    std::vector<VertexId> to_lock = {core.getId(), core.rc().getId()};
    //Here should be collection of all neighbours for one big lock. Then lock of these neighbours, then check that they
    //are still alive. For now it's ok as is since we run flip only in parallel to flip. If this changes and we have to
    //flip in parallel to merge and/or multiplex, we will have to do the advanced locking.
    Locker<VertexId> locker = Locker<VertexId>::FromVector(to_lock);
    return flipVertexLockFree(inc);
}

void ag::AssemblyGraph::replacePathLockFree(const RAGraphPath &old_path, const RAGraphPath &new_path,
                                            size_t dominant_index, const std::unordered_map<EdgeId, EdgeId> &outer) {
    auto set_corporeal = [&](bool value) {
        for (Edge &e : new_path.edges()) {
            e.setCorporeal(value);
            e.rc().setCorporeal(value);
        }
        for (const auto &p : outer) {
            p.second->setCorporeal(value);
            p.second->rc().setCorporeal(value);
        }
    };
    set_corporeal(false);
    fireReplacePath(old_path, new_path, dominant_index, outer);
//    Collected first: the paths reach their vertices through their edges, which isolation deletes.
    std::vector<VertexId> removed;
    for (Vertex &v : old_path.vertices()) {
        bool survives = false;
        for (Vertex &w : new_path.vertices())
            survives = survives || v == w;
        if (!survives)
            removed.push_back(v.getId());
    }
    for (VertexId vid : removed)
        if (!vid->marked())
            isolateAndMarkLockFree(*vid);
    set_corporeal(true);
}

Vertex &ag::AssemblyGraph::mergePathLockFree(const GraphPath &path, Vertex::id_type id) {
    VERIFY(path.startClosed() && path.endClosed());
//    The path is read before new edges make following its codes ambiguous. A cycle is merged without its last edge,
//    which then closes the loop D -> D.
    RAGraphPath old_path = path.asRAPath();
    VERIFY(old_path.size() > 0);
    if (old_path.getStart() == old_path.getFinish())
        old_path = old_path.subPath(0, old_path.size() - 1);
    const size_t n = old_path.size();
    if (n == 0)
        return old_path.getStart();
    for (Vertex &v : old_path.innerVertices()) {
        VERIFY(v.inDeg() == 1 && v.outDeg() == 1);
    }
    Vertex &start = old_path.getStart();
    Vertex &finish = old_path.getFinish();
    const Sequence seq = old_path.Seq();
//    D is the vertex of the path spelling all of it, if there is one (non-redundancy: at most one), else a new vertex.
    VertexId dominant_id;
    for (Vertex &v : old_path.vertices())
        if (!dominant_id.valid() && v.size() == seq.size() && v.getSeq() == seq)
            dominant_id = v.getId();
    Vertex &dominant = dominant_id.valid() ? *dominant_id : addVertex(seq, seq.isCanonical() ? id : -id);
    const bool start_replaced = start != dominant;
    const bool end_replaced = finish != dominant;
    std::vector<EdgeId> incoming;
    std::vector<EdgeId> outgoing;
    if (start_replaced)
        for (Edge &in : start.incoming())
            incoming.push_back(in.getId());
    if (end_replaced)
        for (Edge &out : finish)
            outgoing.push_back(out.getId());
//    y -> start becomes y -> D (D -> D if y is the replaced finish): it spells y -> start continued to the end of D, and
//    its rc overlaps y.rc() as start.rc() did. finish -> z becomes D -> z: D ends with finish, and the rc spells
//    z.rc() -> finish.rc() continued to the end of D.rc(). Types are inherited. An edge already replaced (the loop
//    finish -> start, or the rc of an edge into start for a self-rc path) is not replaced again.
    std::unordered_map<EdgeId, EdgeId> outer;
    auto add_replacement = [&](Edge &old_edge, Vertex &from, Vertex &to, const Sequence &tseq, const Sequence &rctseq) {
        Edge &new_edge = addEdgeLockFree(from, to, tseq, rctseq);
        if (old_edge.getVirtualType() != VirtualType::unset)
            new_edge.setVirtualType(old_edge.getVirtualType());
        outer[old_edge.getId()] = new_edge.getId();
        outer[old_edge.rc().getId()] = new_edge.rc().getId();
    };
    for (EdgeId eid : incoming) {
        const bool from_finish = end_replaced && eid->getStart() == finish;
        add_replacement(*eid, from_finish ? dominant : eid->getStart(), dominant,
                        eid->truncSeq() + seq.Subseq(start.size()),
                        eid->rc().truncSeq() + (from_finish ? dominant.rc().getSeq().Subseq(finish.size()) : Sequence()));
    }
    for (EdgeId eid : outgoing)
        if (outer.find(eid) == outer.end())
            add_replacement(*eid, dominant, eid->getFinish(), eid->truncSeq(),
                            eid->rc().truncSeq() + dominant.rc().getSeq().Subseq(finish.size()));
    replacePathLockFree(old_path, RAGraphPath(dominant), 0, outer);
    return dominant;
}

size_t SmallestShift(const Sequence &seq) {
    Sequence dseq = seq + seq;
    size_t best_shift = 0;
    for(size_t i = 0; i < seq.size(); i++) {
        if(dseq.Subseq(i) < dseq.Subseq(best_shift))
            best_shift = i;
        if (dseq.Subseq(i, i + seq.size()) == !dseq.Subseq(i, i + seq.size())) {
            VERIFY(seq.size() % 2 == 0);
            if (dseq.Subseq(i) < dseq.Subseq(i + seq.size()/2))
                return i;
            else
                return i + seq.size()/2;
        }
    }
    return best_shift;
}

Vertex &ag::AssemblyGraph::mergeLoop(const GraphPath &path) {
    VERIFY(false);
}

Vertex & AssemblyGraph::flipVertexLockFree(Edge &inc) {
    VERIFY(inc.getFinish().inDeg() == 1);
    VERIFY(inc.isSuffix());//Necessary to preserve the splitting property of the graph.
    VERIFY(!inc.isPrefix());
    VERIFY(inc.getStart() != inc.getFinish());
//    The merge of the single edge: inc's start spells it, so it is the merged vertex and core's outgoing edges move to it.
    return mergePathLockFree(GraphPath(inc));
}

Vertex &AssemblyGraph::extendVertexLockFree(Edge &suffix_edge, Vertex::id_type id) {
    Vertex &x = suffix_edge.getStart();
    Vertex &w = suffix_edge.getFinish();
//    Untyped on a plain DBG; the new edges then stay untyped too.
    VERIFY(suffix_edge.isVirtualSuffix() || suffix_edge.isMerge() || suffix_edge.getVirtualType() == VirtualType::unset);
    VERIFY(x.outDeg() == 1);
    VERIFY(!suffix_edge.isPrefix() && !suffix_edge.isSuffix());
    VERIFY(x != w && x != w.rc());
    if (x == x.rc()) {
        VERIFY(w != w.rc());
//        x's only incoming edge is suffix_edge.rc() = w.rc() -> x, which spells rc(t) + x for suffix_edge's tail t, so
//        v = rc(t) + x + t.
        const Sequence seq = suffix_edge.rc().fullSeq() + suffix_edge.truncSeq();
        VERIFY(seq == !seq);
        Vertex &v = addVertex(seq, seq.isCanonical() ? id : -id);
//        w is a suffix of v; the rc w.rc() -> v spells w.rc() -> x -> w, with w.rc() a prefix of v.
        Edge &new_suffix = addEdgeLockFree(v, w, Sequence(), suffix_edge.rc().truncSeq() + suffix_edge.truncSeq());
        if (suffix_edge.getVirtualType() != VirtualType::unset)
            new_suffix.setVirtualType(suffix_edge.getVirtualType());
//        P = w.rc() -> x -> w becomes Q = w.rc() -> v -> w with D = v. Both are their own rc and P's ends survive.
        Vertex &u = w.rc();
        replacePathLockFree(RAGraphPath(u, {suffix_edge.rc().getId(), suffix_edge.getId()}, 0, 0),
                            RAGraphPath(u, {new_suffix.rc().getId(), new_suffix.getId()}, 0, 0), 1, {});
        return v;
    }
    std::vector<EdgeId> incoming;
    for (Edge &in : x.incoming()) {
        VERIFY(in.getStart() != x && in.getStart() != x.rc());
        incoming.push_back(in.getId());
    }
    Sequence seq = suffix_edge.fullSeq();
    Vertex &v = addVertex(seq, seq.isCanonical() ? id : -id);
    std::unordered_map<EdgeId, EdgeId> extensions;
//    y -> v spells y -> x -> w; v.rc() -> y.rc() overlaps y.rc() exactly as x.rc() did.
    for (EdgeId eid : incoming) {
        Edge &new_edge = addEdgeLockFree(eid->getStart(), v, eid->truncSeq() + suffix_edge.truncSeq(),
                                         eid->rc().truncSeq());
        if (eid->getVirtualType() != VirtualType::unset)
            new_edge.setVirtualType(eid->getVirtualType());
        extensions[eid] = new_edge.getId();
        extensions[eid->rc().getId()] = new_edge.rc().getId();
    }
//    w is a suffix of v; w.rc() -> v.rc() spells w.rc() -> x.rc().
    Edge &new_suffix = addEdgeLockFree(v, w, Sequence(), suffix_edge.rc().truncSeq());
    if (suffix_edge.getVirtualType() != VirtualType::unset)
        new_suffix.setVirtualType(suffix_edge.getVirtualType());
//    P = x -> w becomes Q = v -> w; x's incoming edges are replaced.
    replacePathLockFree(RAGraphPath(suffix_edge), RAGraphPath(new_suffix), 0, extensions);
    return v;
}

//    A real prefix or suffix edge is virtually so.
static void SetSPVirtualType(Edge &edge) {
    edge.setVirtualType(edge.isPrefix() && edge.isSuffix() ? VirtualType::merge :
                        edge.isSuffix() ? VirtualType::suffix : VirtualType::prefix);
}

Edge &ag::AssemblyGraph::addSPEdgeLockFree(Vertex &start, Vertex &end,
                                           ag::Edge::id_type eid,
                                           ag::Edge::id_type rcid) {
    Sequence tseq = end.getSeq().Subseq(std::min(start.size(), end.size()));
    Sequence rctseq = start.rc().getSeq().Subseq(std::min(start.size(), end.size()));
    Edge &res = addEdgeLockFree(start, end, tseq, rctseq, EdgeData(), eid, rcid);
    SetSPVirtualType(res);
    return res;
}

Edge &ag::AssemblyGraph::addSPEdgeLockFree(Vertex &start, Vertex &end, VirtualType type) {
    Sequence tseq = end.getSeq().Subseq(std::min(start.size(), end.size()));
    Sequence rctseq = start.rc().getSeq().Subseq(std::min(start.size(), end.size()));
    Edge &res = addEdgeLockFree(start, end, tseq, rctseq);
    res.setVirtualType(type);
    return res;
}

Edge & ag::AssemblyGraph::addSPEdge(Vertex &start, Vertex &end, ag::Edge::id_type eid,
                                         ag::Edge::id_type rcid) {
    Sequence tseq = end.getSeq().Subseq(std::min(start.size(), end.size()));
    Sequence rctseq = start.rc().getSeq().Subseq(std::min(start.size(), end.size()));
    Edge &res = addEdge(start, end, tseq, rctseq, EdgeData(), eid, rcid);
    SetSPVirtualType(res);
    return res;
}

Edge &AssemblyGraph::addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeIdType eid,
                             EdgeIdType rcid) {
    return addEdge(start, end, tseq, rctseq, EdgeData(), eid, rcid);
}

Vertex &AssemblyGraph::edgeToSupreVertex(Edge &edge, Vertex::id_type id) {
    Locker<VertexId> locker = Locker<VertexId>::FromVector({edge.getStart().getId(), edge.getFinish().rc().getId()});
    return edgeToSupreVertexLockFree(edge, id);
}

Vertex &AssemblyGraph::edgeToSupreVertexLockFree(Edge &edge, Vertex::id_type id) {
    VERIFY(!edge.isSuffix());
    VERIFY(!edge.isPrefix());
    Sequence seq = edge.fullSeq();
    Vertex &res = addVertex(seq, seq.isCanonical() ? id : -id);
    VERIFY(res.size() == res.rc().size());
    addEdgeLockFree(res, edge.getFinish(), Sequence(), edge.rc().truncSeq());
    if(edge != edge.rc())
        addEdgeLockFree(res.rc(), edge.getStart().rc(), Sequence(), edge.truncSeq());
//    res -> finish and res.rc() -> start.rc() are suffix edges: virtually suffix if the edge they take over from was
//    (edge, and edge.rc() for the second; untyped edges, during the conversion of a DBG, count as such), merge otherwise.
    const bool untyped = edge.getVirtualType() == VirtualType::unset;
    const VirtualType out_type = untyped || edge.isVirtualSuffix() ? VirtualType::suffix : VirtualType::merge;
    const VirtualType rc_out_type = untyped || edge.rc().isVirtualSuffix() ? VirtualType::suffix : VirtualType::merge;
    res.front().setVirtualType(out_type);
    if (res.rc().front().getVirtualType() == VirtualType::unset)
        res.rc().front().setVirtualType(rc_out_type);
    res.front().setCorporeal(false);
    res.front().rc().setCorporeal(false);
    res.rc().front().setCorporeal(false);
    res.incFront().setCorporeal(false);
    fireEdgeToSupreVertex(res, edge);
    removeEdgeLockFree(edge);
    res.front().setCorporeal(true);
    res.front().rc().setCorporeal(true);
    res.rc().front().setCorporeal(true);
    res.incFront().setCorporeal(true);
    return res;
}

//    connections are the plan's unique connections, taken before ExtendAdjacent: selecting them dereferences the edges,
//    which ExtendAdjacent may free. Only their stable ids are used here; the mapped extension edges are always live.
VertexResolutionPlan FixPlan(Vertex &core, const std::vector<InOutEdgePair> &connections,
                             const std::unordered_map<EdgeId, EdgeId> &extensions) {
    VertexResolutionPlan new_plan(core);
    for (const InOutEdgePair &p : connections)
        new_plan.add({*extensions.at(p.incomingId()), *extensions.at(p.outgoingId()), p.getSupport()});
    return new_plan;
}

std::unordered_map<EdgeId, EdgeId> AssemblyGraph::ExtendAdjacent(Vertex &core) {
    std::unordered_map<EdgeId, EdgeId> extensions;
    for (Vertex &v : ThisAndRC(core)) {
        for (Edge &edge : v.incoming()) {
            extensions[edge.getId()] = edge.getId();
            extensions[edge.rc().getId()] = edge.rc().getId();
        }
    }
    // Every corporeal edge entering `core` must become a suffix edge into core, so that after the
    // resolution core's whole neighbourhood is suffix/prefix structured. We must (a) snapshot the edges to
    // extend before mutating the graph - the incoming set itself changes as new suffix edges appear - and
    // (b) capture each edge's rc id before the operation frees it, since the old handles are dangling
    // afterwards. The extensions map is keyed by edge id (comparison does not dereference), so stale keys
    // stay valid as keys; only the mapped value is updated, always to the live edge.
    std::vector<std::pair<EdgeId, EdgeId>> to_extend;
    for (Vertex &v : ThisAndRC(core)) {
        for (Edge &edge : v.incoming()) {
            if (!edge.isSuffix())
                to_extend.emplace_back(edge.getId(), edge.rc().getId());
        }
    }
//    core has several incoming edges, so each x -> core is x's only outgoing edge (virtually suffix): extendVertex
//    replaces x by v, x extended to the end of core, which keeps x's segments and edges and enters core by a suffix
//    edge. A self-rc x (a hairpin core.rc() -> x -> core) is extended both ways at once. Any other loop body x
//    (core -> x -> core, or a hairpin through core.rc()) is extended twice: forward along x -> core, then backward
//    along the rc of the edge from core, which the first extension replaced by core -> v. So the edges of the map that
//    end at x (start at x.rc()) follow the extension to v (v.rc()).
    for (const std::pair<EdgeId, EdgeId> &key : to_extend) {
        Edge &edge = *extensions.at(key.first);
        if (edge.isSuffix())
            continue;
        Vertex &x = edge.getStart();
        std::vector<std::pair<EdgeId, VertexId>> into_x;
        std::vector<std::pair<EdgeId, VertexId>> from_x_rc;
        for (const auto &p : extensions) {
            if (p.second != edge.getId() && p.second->getFinish() == x)
                into_x.emplace_back(p.first, p.second->getStart().getId());
            if (p.second != edge.rc().getId() && p.second->getStart() == x.rc())
                from_x_rc.emplace_back(p.first, p.second->getFinish().getId());
        }
        Edge &suffix_edge = extendVertexLockFree(edge).front();
        extensions.at(key.first) = suffix_edge.getId();
        extensions.at(key.second) = suffix_edge.rc().getId();
        if (!x.marked())
            continue;
        Vertex &v = suffix_edge.getStart();
        for (const auto &p : into_x)
            for (Edge &e : *p.second)
                if (e.getFinish() == v)
                    extensions.at(p.first) = e.getId();
        for (const auto &p : from_x_rc)
            for (Edge &e : v.rc())
                if (e.getFinish() == *p.second)
                    extensions.at(p.first) = e.getId();
    }
    return std::move(extensions);
}

VertexResolutionResult
AssemblyGraph::resolveCoreVertex(Vertex &core, const VertexResolutionPlan &resolution) {
    //TODO. make it work in parallel. Probably lock vertex and neighbours.
    VERIFY(core.isCore());
    VERIFY(core.inDeg() > 1);
    VERIFY(core.outDeg() > 1);
    for (Edge &edge : core.incoming()) {VERIFY(edge.isSuffix());}
    for (Edge &edge : core) {VERIFY(edge.isPrefix());}
//    Every neighbour of core has sample segments, so read starts moved off core land on vertices with segments.
    for (Edge &edge : core.incoming()) {VERIFY(!edge.getStart().coverage_info.empty());}
    for (Edge &edge : core) {VERIFY(!edge.getFinish().coverage_info.empty());}
    VertexResolutionResult result(core);
    std::vector<EdgeId> new_edges;
    for(const InOutEdgePair &p : resolution.connectionsUnique()) {
        VERIFY(p.incoming().getFinish() == core);
        VERIFY(p.outgoing().getStart() == core);
//        If exactly one of the edges has a single partner, a new vertex would be equivalent to its neighbour on that side
//        (and would get no sample segments), so the pair is resolved by a direct edge instead. A 1-to-1 pair gets a copy
//        of core: u -> copy is a suffix edge, copy -> w a prefix edge. A copy of a self-rc core is self-rc, so it can only
//        resolve a pair that is its own rc; any other 1-to-1 pair at a self-rc core gets a direct edge.
        const bool single_in = resolution.incEdgeConnections(p.incoming()) == 1;
        const bool single_out = resolution.outEdgeConnections(p.outgoing()) == 1;
        const bool one_to_one = single_in && single_out;
        if (single_in != single_out || (one_to_one && core == core.rc() && p.incoming() != p.outgoing().rc())) {
            Edge &direct = addEdgeLockFree(p.incoming().getStart(), p.outgoing().getFinish(), p.outgoing().truncSeq(),
                                           p.incoming().rc().truncSeq());
//            u keeps a single outgoing edge (virtually suffix) if its edge has a single partner, w a single incoming one
//            (virtually prefix) if its edge has; both: the edge is to be merged.
            direct.setVirtualType(one_to_one ? VirtualType::merge :
                                  single_in ? VirtualType::suffix : VirtualType::prefix);
            result.add(direct, p);
            new_edges.push_back(direct.getId());
            direct.setCorporeal(false);
            direct.rc().setCorporeal(false);
        } else if (one_to_one) {
            Vertex &copy = addVertex(core.getSeq());
            result.add(copy, p);
            Edge &inc = addEdgeLockFree(p.incoming().getStart(), copy, p.incoming().truncSeq(), p.incoming().rc().truncSeq());
            inc.setVirtualType(VirtualType::suffix);
            new_edges.push_back(inc.getId());
            inc.setCorporeal(false);
            inc.rc().setCorporeal(false);
//            For a self-rc copy, copy -> w is the rc of u -> copy.
            if (copy != copy.rc()) {
                Edge &out = addEdgeLockFree(copy, p.outgoing().getFinish(), p.outgoing().truncSeq(), p.outgoing().rc().truncSeq());
                out.setVirtualType(VirtualType::prefix);
                new_edges.push_back(out.getId());
                out.setCorporeal(false);
                out.rc().setCorporeal(false);
            }
        } else {
            Vertex &newv = addVertex(p.getSeq());
            result.add(newv, p);
            Edge &out = addEdgeLockFree(p.incoming().getStart(), newv, p.outgoing().truncSeq().Subseq(0, 1), p.incoming().rc().truncSeq().Subseq(1));
            out.setVirtualType(VirtualType::prefix);
            new_edges.push_back(out.getId());
            out.setCorporeal(false);
            out.rc().setCorporeal(false);
            if(newv != newv.rc()) {
                Edge &inc = addEdgeLockFree(p.outgoing().getFinish().rc(), newv.rc(), p.incoming().rc().truncSeq().Subseq(0, 1), p.outgoing().truncSeq().Subseq(1));
                inc.setVirtualType(VirtualType::prefix);
                new_edges.push_back(inc.getId());
                inc.setCorporeal(false);
                inc.rc().setCorporeal(false);
            }
        }
    }
    fireResolveVertex(core, result);
    isolateAndMarkLockFree(core);
    for (EdgeId eid : new_edges) {
        eid->setCorporeal(true);
        eid->rc().setCorporeal(true);
    }
    return std::move(result);
}

VertexResolutionResult
AssemblyGraph::resolveVertex(Vertex &core, const VertexResolutionPlan &resolution) {
    //TODO. make it work in parallel. Probably lock vertex and neighbours.
    VERIFY(core.isCore());
    VERIFY(core.inDeg() > 1);
    VERIFY(core.outDeg() > 1);
    std::vector<InOutEdgePair> connections;
    for (const InOutEdgePair &p : resolution.connectionsUnique())
        connections.push_back(p);
    std::unordered_map<EdgeId, EdgeId> extensions = ExtendAdjacent(core);
    VertexResolutionPlan plan = FixPlan(core, connections, extensions);
    return resolveCoreVertex(core, plan);
}

IterableStorage<TransformingIterator<typename std::list<Edge>::const_iterator, const Edge>> Vertex::incoming() const {
    std::function<const Edge &(const Edge &)> transform = [](const Edge &edge) -> const Edge & {
        return edge.rc();
    };
    return {{rc().begin(), rc().end(), transform}, {rc().end(), rc().end(), transform}};
}

void Vertex::setRC(Vertex &other) {
    rc_ = &other;
    rc_->rc_ = this;
}

Vertex::Vertex(Vertex::id_type id, Sequence seq, VertexData data) : VertexData(std::move(data)), id(id), seq(std::move(seq)), rc_(nullptr) {
    canonical = this->seq <= !this->seq;
    omp_init_lock(&writelock);
}

Vertex::Vertex(Vertex::id_type id, bool canonical, VertexData data) : VertexData(std::move(data)), id(id), seq(), canonical(canonical), rc_(nullptr) {
    omp_init_lock(&writelock);
}

IterableStorage<TransformingIterator<typename std::list<Edge>::iterator, Edge>> Vertex::incoming() {
    std::function<Edge &(Edge &)> transform = [](Edge &edge) -> Edge & {
        return edge.rc();
    };
    return {{rc().begin(), rc().end(), transform}, {rc().end(), rc().end(), transform}};
}

void Vertex::updateMaxOutId(const std::array<int, 5> other) {
    for(size_t i = 0; i < max_out_id.size(); i++)
        updateMaxOutId(i + other[i] * 10);
}

bool Vertex::innerRemoveEdge(Edge &edge) {
    VERIFY(edge.getStart() == *this);
    if(edge._rc != nullptr) {
        edge._rc->_rc = nullptr;
    }
    auto it = std::find(outgoing_.begin(), outgoing_.end(), edge);
    if(it == outgoing_.end())
        return false;
//    Non-corporeal edges are never removed (see AssemblyGraph::removeEdgeLockFree).
    VERIFY(it->corporeal);
    _outDeg--;
    outgoing_.erase(it);
    return true;
}

bool Vertex::operator<(const Vertex &other) const {
    return (std::abs(id) << 1) + (id > 0) < (std::abs(other.id) << 1) + (other.id > 0);
}

bool Vertex::operator<=(const Vertex &other) const {
    return *this < other || *this == other;
}

bool Vertex::operator>(const Vertex &other) const {
    return (std::abs(id) << 1) + (id > 0) > (std::abs(other.id) << 1) + (other.id > 0);
}

bool Vertex::operator>=(const Vertex &other) const {
    return *this > other || *this == other;
}

void Vertex::clear() {
    outgoing_.clear();
    _outDeg = 0;
}

bool Vertex::isCore() const {
    return (inDeg() == 0 || !rc().front().isSuffix()) && (outDeg() == 0 || !front().isSuffix());
}

bool Vertex::isOuter() const {
    return (inDeg() == 0 || rc().front().isSuffix()) && (outDeg() == 0 || front().isSuffix());
}


Vertex &AssemblyGraph::addSelfRCVertex(VertexData data) {
    Vertex *res;
#pragma omp critical
    {
        typename Vertex::id_type id = maxVId + 1;
        Vertex &new_res = innerAddVertex(id, true, std::move(data));
        new_res.setRC(new_res);
        res = &new_res;
    }
    this->fireAddVertex(*res);
    return *res;
}

IterableStorage<SkippingIterator<typename AssemblyGraph::vertex_iterator_type>> AssemblyGraph::vertices(bool unique) & {
    std::function<bool(Vertex &)> use =
            [unique](Vertex &vertex) -> bool {
                return !unique || vertex.isCanonical();
            };
    SkippingIterator<vertex_iterator_type> begin(vertex_list.begin(), vertex_list.end(), use);
    SkippingIterator<vertex_iterator_type> end(vertex_list.end(), vertex_list.end(), use);
    return {begin, end};
}

IterableStorage<SkippingIterator<typename AssemblyGraph::const_vertex_iterator_type>> AssemblyGraph::vertices(bool unique) const & {
    std::function<bool(const Vertex &)> use =
            [unique](const Vertex &vertex) -> bool {
                return !unique || vertex.isCanonical();
            };
    SkippingIterator<const_vertex_iterator_type> begin(vertex_list.begin(), vertex_list.end(), use);
    SkippingIterator<const_vertex_iterator_type> end(vertex_list.end(), vertex_list.end(), use);
    return {begin, end};
}

IterableStorage<SkippingIterator<typename AssemblyGraph::vertex_iterator_type>> AssemblyGraph::verticesUnique() &{
    return vertices(true);
}

IterableStorage<SkippingIterator<typename AssemblyGraph::const_vertex_iterator_type>> AssemblyGraph::verticesUnique() const &{
    return vertices(true);
}

IterableStorage<ApplyingIterator<typename AssemblyGraph::vertex_iterator_type, Edge, 4>> AssemblyGraph::edges(bool unique) & {
    std::function<std::array<Edge*, 4>(Vertex &)> apply = [unique](Vertex &vertex) {
        if(vertex.outDeg() > 4) {
            std::cout << vertex.getSeq() << std::endl;
            for(Edge &e : vertex) {
                std::cout << e.truncSeq() << std::endl;
            }
        }
        VERIFY(vertex.outDeg() <= 4);
        std::array<Edge*, 4> res = {};
        size_t cur = 0;
        for(Edge &edge : vertex) {
            if(!unique || edge <= edge.rc()) {
                res[cur] = &edge;
                cur++;
            }
        }
        return res;
    };
    ApplyingIterator<vertex_iterator_type, Edge, 4> begin(vertex_list.begin(), vertex_list.end(), apply);
    ApplyingIterator<vertex_iterator_type, Edge, 4> end(vertex_list.end(), vertex_list.end(), apply);
    return {begin, end};
}

IterableStorage<ApplyingIterator<typename AssemblyGraph::const_vertex_iterator_type, const Edge, 4>> AssemblyGraph::edges(bool unique) const & {
    std::function<std::array<const Edge*, 4>(const Vertex &)> apply = [unique](const Vertex &vertex) {
        std::array<const Edge*, 4> res = {};
        size_t cur = 0;
        for(const Edge &edge : vertex) {
            if(!unique || edge <= edge.rc()) {
                res[cur] = &edge;
                cur++;
            }
        }
        return res;
    };
    ApplyingIterator<const_vertex_iterator_type, const Edge, 4> begin(vertex_list.begin(), vertex_list.end(), apply);
    ApplyingIterator<const_vertex_iterator_type, const Edge, 4> end(vertex_list.end(), vertex_list.end(), apply);
    return {begin, end};
}

IterableStorage<ApplyingIterator<typename AssemblyGraph::vertex_iterator_type, Edge, 4>> AssemblyGraph::edgesUnique() &{
    return edges(true);
}

IterableStorage<ApplyingIterator<typename AssemblyGraph::const_vertex_iterator_type, const Edge, 4>> AssemblyGraph::edgesUnique() const &{
    return edges(true);
}

void AssemblyGraph::removeIsolated() {
    vertex_storage_type newv;
    for(Vertex &v: verticesUnique())
        if(v.inDeg() == 0 && v.outDeg() == 0)
            this->isolateAndMarkLockFree(v);
    removeMarked();
}

//    TODO: move listening to marking. Rename marking into deleting, make it private.
void AssemblyGraph::removeMarked() {
    for (auto it = vertex_list.begin(); it != vertex_list.end();) {
        if (it->marked()) {
            this->fireDeleteVertex(*it);
            if(*it == it->rc()) {
                it = vertex_list.erase(it);
            } else {
                Vertex &rc = it->rc();
                it = vertex_list.erase(it);
                VERIFY(rc == *it);//We rely on the fact that vertex and its rc are adjacent in the list.
                it = vertex_list.erase(it);
            }
        } else {
            ++it;
        }
    }
}

void AssemblyGraph::resetMarkers() {
    for(Vertex &vertex : vertices()) {
        for(Edge &edge : vertex) {
            edge.mark(EdgeMarker::common);
        }
        vertex.unmark();
    }
}

Vertex::id_type AssemblyGraph::reserveVertexIds(size_t num) {
    maxVId += num;
    return maxVId - num + 1;
}

//    Vertex id selection and vertex_list insertion have no per-vertex id to lock on yet, so this critical
//    section is what makes concurrent addVertex/addVertexPair/addSelfRCVertex calls safe. Contract: no one
//    iterates over the graph's vertex list while vertices are being added.
Vertex &AssemblyGraph::addVertex(const Sequence &seq, const VertexData &data, typename Vertex::id_type id) {
    Vertex *res;
#pragma omp critical
    {
        if (id == 0) {
            if(seq <= !seq)
                id = maxVId + 1;
            else
                id = -maxVId - 1;
        }
        maxVId = std::max(std::abs(id), maxVId);
        Vertex &new_res = innerAddVertex(id, seq, data);
        Vertex &new_rc = seq == !seq ? new_res : innerAddVertex(-id, !seq, data.RC());
        new_res.setRC(new_rc);
        // new_res.setSeq(seq);
        res = &new_res;
    }
    this->fireAddVertex(*res);
    return *res;
}

size_t AssemblyGraph::edgeCount() const {
    size_t res = 0;
    for(auto &v : vertices())
        res += v.outDeg();
    return res;
}

AssemblyGraph::~AssemblyGraph() {
    for(Vertex &v : verticesUnique())
        if(!v.marked())
            isolateAndMarkLockFree(v);
    removeMarked();
}

Vertex &AssemblyGraph::addVertex(const Vertex &other_graph_vertex) {
    Vertex &res = addVertex(other_graph_vertex.getSeq(), other_graph_vertex, other_graph_vertex.getInnerId());
    res.updateMaxOutId(other_graph_vertex.getMaxOutId());
    res.rc().updateMaxOutId(other_graph_vertex.rc().getMaxOutId());
    return res;
}

Vertex &AssemblyGraph::innerAddVertex(typename Vertex::id_type id, bool canonical, VertexData data) {
    VERIFY(canonical == (id > 0));
    maxVId = std::max(std::abs(id), maxVId);
    vertex_list.emplace_back(id, canonical, std::move(data));
    return vertex_list.back();
}

bool id_cmp(Vertex::id_type a, Vertex::id_type b) {
    if (std::abs(a) != std::abs(b))
        return std::abs(a) < std::abs(b);
    return a > b;
}
Vertex &AssemblyGraph::innerAddVertex(typename Vertex::id_type id, Sequence seq, VertexData data) {
    VERIFY(seq.isCanonical() == (id > 0));
    maxVId = std::max(std::abs(id), maxVId);
    if (vertex_list.empty() || id_cmp(vertex_list.back().getInnerId(), id)) {
        vertex_list.emplace_back(id, std::move(seq), std::move(data));
        return vertex_list.back();
    } else if (id_cmp(id, vertex_list.front().getInnerId())) {
        vertex_list.emplace_front(id, std::move(seq), std::move(data));
        return vertex_list.front();
    } else {
        auto it = vertex_list.end();
        --it;
        while (id_cmp(id, it->getInnerId())) {
            VERIFY(it != vertex_list.begin());
            --it;
        }
        ++it;
        return *vertex_list.emplace(it, id, std::move(seq), std::move(data));
    }
}

Edge &AssemblyGraph::addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeData data,
                                     EdgeIdType eid, EdgeIdType rcid) {
    return addEdgeLockFree(start, end, full_sequence.Subseq(start.size()), full_sequence.rc().Subseq(end.size()), data, eid, rcid);
}

Edge &AssemblyGraph::addEdgeLockFree(Vertex &start, Vertex &end,
                                     const Sequence &tseq, const Sequence &rctseq,
                                     EdgeData data, EdgeIdType eid, EdgeIdType rcid) {
    for(Edge &edge: start) {
        if(edge.getFinish() == end && edge.truncSeq() == tseq) {
            return edge;
        }
    }
    Edge &res = start.innerAddEdge(end, tseq, data, eid);
    if(end.rc() != start || (tseq.size() > start.size() && tseq.rc().Subseq(end.size()) != rctseq.rc().Subseq(start.size()))) {
        Edge &rc_edge = end.rc().innerAddEdge(start.rc(), rctseq, data.RC(), rcid);
        res._rc = &rc_edge;
        rc_edge._rc = &res;
    } else {
        res._rc = &res;
    }
    VERIFY(res.fullSize() == res.rc().fullSize());
    this->fireAddEdge(res);
    return res;
}

void AssemblyGraph::removeEdgeLockFree(Edge &edge) {
    // Only corporeal (already-existing) edges are ever removed: graph operations add new edges as
    // non-corporeal, run listeners, remove old (corporeal) edges, then mark the new ones corporeal.
    VERIFY(edge.corporeal);
    this->fireDeleteEdge(edge);
    if(edge != edge.rc())
        edge.getFinish().rc().innerRemoveEdge(edge.rc());
    edge.getStart().innerRemoveEdge(edge);
}

void AssemblyGraph::removeEdge(Edge &edge) {
    Locker<VertexId> locker = Locker<VertexId>::FromVector({edge.getStart().getId(), edge.getFinish().rc().getId()});
    removeEdgeLockFree(edge);
}

Edge &AssemblyGraph::addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeData data, EdgeIdType eid, EdgeIdType rcid) {
    return addEdge(start, end, full_seq.Subseq(start.size()), full_seq.rc().Subseq(end.size()), data, eid, rcid);
}

Edge &AssemblyGraph::addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq,
                                              EdgeData data, EdgeIdType eid, EdgeIdType rcid) {
    Locker<VertexId> locker = Locker<VertexId>::FromVector({start.getId(), end.rc().getId()});
    return addEdgeLockFree(start, end, tseq, rctseq, std::move(data), eid, rcid);
}

void AssemblyGraph::isolateAndMarkLockFree(Vertex &vertex) {
    VERIFY(!vertex.marked());
    for(Vertex &v : ThisAndRC(vertex)) {
        while (v.outDeg() != 0) {
            removeEdgeLockFree(v.front());
        }
        v.mark();
    }
}

namespace {
//    The vertices whose edge lists isolating vertex changes: vertex, vertex.rc() and, for every edge e leaving either,
//    e.getFinish().rc() (which holds e.rc()). Sorted and unique.
    std::vector<VertexId> IsolationLockSet(Vertex &vertex) {
        std::vector<VertexId> res = {vertex.getId(), vertex.rc().getId()};
        for (Vertex &v : ThisAndRC(vertex))
            for (Edge &e : v)
                res.emplace_back(e.getFinish().rc().getId());
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        return res;
    }
}

void AssemblyGraph::isolateAndMark(Vertex &vertex) {
//    The lock set is read while vertex and vertex.rc() are locked, then taken in full. A concurrent isolation of a
//    neighbour may have changed it in between; then it is read again.
    while (true) {
        std::vector<VertexId> lock_set;
        {
            Locker<VertexId> own = Locker<VertexId>::FromVector({vertex.getId(), vertex.rc().getId()});
            lock_set = IsolationLockSet(vertex);
        }
        Locker<VertexId> locker = Locker<VertexId>::FromVector(lock_set);
        if (IsolationLockSet(vertex) == lock_set) {
            isolateAndMarkLockFree(vertex);
            return;
        }
    }
}

Edge &AssemblyGraph::mergePathToEdge(const GraphPath &path) {
    Locker<VertexId> locker = Locker<VertexId>::FromVector({path.getStart().getId(), path.getStart().rc().getId(), path.getRCStart().getId(), path.getFinish().rc().getId()});
    return mergePathToEdgeLockFree(path);
}
//    TODO: when TAGraphPath is thoroughly damned, put it back here as a parameter.
Edge &AssemblyGraph::mergePathToEdgeLockFree(const GraphPath &path) {
    VERIFY(!path.empty());
    VERIFY(path.endClosed() && path.startClosed());
    for(Vertex &v : path.innerVertices()) {
        VERIFY(!v.marked());
    }
    VERIFY(!path.isSingleton());
//    Only used in DBG now.
    VERIFY(path.truncLen() != 0 && path.RC().truncLen() != 0);
    // VERIFY(path.getStart() == path.getFinish() || path.getStart() == path.getFinish().rc() || (path.getStart().isJunction() && path.getFinish().isJunction()));
    SequenceBuilder sb;
    Sequence new_seq = path.Seq();
    Edge &new_edge = addEdgeLockFree(path.getStart(), path.getFinish(), new_seq);
    new_edge.setCorporeal(false);
    new_edge.rc().setCorporeal(false);
    VERIFY((path == path.RC()) == (new_edge == new_edge.rc()));
    this->fireMergePathToEdge(path.asRAPath(), new_edge);
    std::vector<VertexId> inner_vertices;
    for(Vertex & v: path.innerVertices()){
        inner_vertices.emplace_back(v.getId());
    }
    for(VertexId &vid : inner_vertices) {
        if(!vid->marked()||vid->outDeg() > 0 || vid->inDeg() > 0)
            isolateAndMarkLockFree(*vid);
    }
    new_edge.setCorporeal(true);
    new_edge.rc().setCorporeal(true);
    return new_edge;
}


GraphPath AssemblyGraph::splitEdge(Edge &edge, const std::vector<EdgePosition> &split_positions) {
    VERIFY(!split_positions.empty());
    Locker<VertexId> locker = Locker<VertexId>::FromVector({edge.getStart().getId(), edge.getFinish().getId()});
    VERIFY(split_positions.front().pos > 0);
    VERIFY(split_positions.back().pos < edge.truncSize());
    for(size_t i = 0; i + 1 < split_positions.size(); i++)
        VERIFY(split_positions[i].pos < split_positions[i+1].pos);
    RAGraphPath res;
    VertexId last_vertex = edge.getStart().getId();
    EdgePosition last_pos = EdgePosition(edge, 0);
    bool self_rc = edge == edge.rc();
    if(self_rc) {
        for(size_t i = 0; i < split_positions.size(); i++) {
            VERIFY(split_positions[i].pos + split_positions[split_positions.size() - 1 - i].pos == edge.truncSize())
        }
    }
    for(EdgePosition pos: split_positions) {
        if(!self_rc || pos.pos * 2 <= edge.truncSize()) {
            Vertex &new_vertex = addVertex(pos.kmerSeq());
            Edge & new_edge = addEdgeLockFree(*last_vertex, new_vertex,
                                              edge.truncSeq().Subseq(last_pos.pos, pos.pos),
                                              edge.rc().truncSeq().Subseq(pos.RC().pos, last_pos.RC().pos));
            new_edge.setCorporeal(false);
            new_edge.rc().setCorporeal(false);
            if(new_vertex == new_vertex.rc() || new_edge == new_edge.rc()) {
                self_rc = true;
            }
            last_vertex = new_vertex.getId();
            res += new_edge;
        } else {
            if(last_vertex->outDeg() == 0) {
                VERIFY(last_pos.pos * 2 < edge.truncSize());
                addEdgeLockFree(*last_vertex, last_vertex->rc(), edge.truncSeq().Subseq(last_pos.pos, pos.pos),
                                edge.rc().truncSeq().Subseq(pos.RC().pos, last_pos.RC().pos));
            }
            Edge &new_edge = last_vertex->front();
            res += new_edge;
            last_vertex = new_edge.getFinish().getId();
        }
        last_pos = pos;
    }
    if(!self_rc) {
        res += addEdgeLockFree(*last_vertex, edge.getFinish(),
                               edge.truncSeq().Subseq(last_pos.pos),
                               edge.rc().truncSeq().Subseq(0, last_pos.RC().pos));
        res.backEdge().setCorporeal(false);
        res.backEdge().rc().setCorporeal(false);
    } else {
        res += last_vertex->front();
    }
    this->fireSplitEdge(edge, res);
    removeEdgeLockFree(edge);
    for(Edge &e : res.edges()) {
        e.setCorporeal(true);
        e.rc().setCorporeal(true);
    }
    return GraphPath(res);
}

AlignmentForm::ConstAlignmentColumnIterator ChooseSplitColumn(Edge &leftEdge, Edge &rightEdge, const AlignmentForm &alignment) {
    Sequence left_from = leftEdge.fullSubseq(leftEdge.fullSize() - alignment.queryLength(), leftEdge.fullSize());
    Sequence right_to = rightEdge.fullSubseq(0, alignment.targetLength());
    std::vector<AlignmentForm::ConstAlignmentColumnIterator> columns;
    for(AlignmentForm::ConstAlignmentColumnIterator it = alignment.columns().begin(); it != alignment.columns().end(); ++it) {
        if(it.getQpos() + leftEdge.truncSize() >= alignment.queryLength() && it.getTpos() <= rightEdge.rc().truncSize()) {
            columns.emplace_back(it);
        }
    }
    AlignmentForm::ConstAlignmentColumnIterator r_col = columns[columns.size() / 2];
    AlignmentForm::ConstAlignmentColumnIterator l_col = columns[columns.size() / 2];
    while (true) {
        VERIFY(r_col.getQpos() + 1 < left_from.size());
        VERIFY(r_col != alignment.columns().end());
        AlignmentForm::AlignmentColumn col = *r_col;
        if (col.event == CigarEvent::M && left_from[col.qpos] == right_to[col.tpos] && left_from[col.qpos - 1] != right_to[col.tpos + 1]) {
            return r_col;
        }
        ++r_col;;
        VERIFY(l_col.getTpos() > 0);
        VERIFY(l_col != alignment.columns().begin());
        col = *l_col;
        if (col.event == CigarEvent::M && left_from[col.qpos] == right_to[col.tpos] && left_from[col.qpos - 1] != right_to[col.tpos + 1]) {
            return l_col;
        }
        --l_col;
    }
    return {};
}

Edge &AssemblyGraph::mergeTipsToEdge(Edge &leftEdge, Edge &rightEdge, AlignmentForm alignment) {
    VERIFY(leftEdge != rightEdge);
    VERIFY(alignment.queryLength() <= leftEdge.fullSize());
    VERIFY(alignment.targetLength() <= rightEdge.rc().fullSize());
    VERIFY(leftEdge.getFinish().outDeg() == 0);
    VERIFY(leftEdge.getFinish().inDeg() == 1);
    VERIFY(rightEdge.getStart().outDeg() == 1);
    VERIFY(rightEdge.getStart().inDeg() == 0);
//    DBG only: the tip vertices removed are k-mers (and so have no sample segments).
    VERIFY(leftEdge.getFinish().size() == k && rightEdge.getStart().size() == k);
    auto split_column = ChooseSplitColumn(leftEdge, rightEdge, alignment);
    AlignmentForm right_sub_alignment(split_column, alignment.columns().end());
    AlignmentForm left_sub_alignment(alignment.columns().begin(), split_column);
    Sequence left_part = leftEdge.fullSubseq(0, leftEdge.fullSize() - alignment.queryLength() + split_column.getQpos());
    Sequence right_part = rightEdge.fullSubseq(split_column.getTpos(), rightEdge.fullSize());
//        Making sure that the sequence is properly collapsed
    VERIFY(left_part[left_part.size() - 1] != right_part[0]);
    VERIFY(left_part[left_part.size() - 1] != right_part[1]);
    Sequence new_seq = left_part + right_part;
    AlignmentForm left_alignment = AlignmentForm::Equal(left_sub_alignment.queryLength()) + right_sub_alignment;
    AlignmentForm right_alignment = left_sub_alignment.Reverse() + AlignmentForm::Equal(right_sub_alignment.targetLength());
    Edge &new_edge = addEdge(leftEdge.getStart(), rightEdge.getFinish(), new_seq);
    new_edge.setCorporeal(false);
    new_edge.rc().setCorporeal(false);
    this->fireMergeTipsToEdge(new_edge, leftEdge, rightEdge, left_alignment, right_alignment);
    if(leftEdge != rightEdge.rc())
        isolateAndMarkLockFree(rightEdge.getStart());
    isolateAndMarkLockFree(leftEdge.getFinish());
    new_edge.setCorporeal(true);
    new_edge.rc().setCorporeal(true);
    return new_edge;
}

//    TODO: Make it run in parallel
void AssemblyGraph::resetEdgeCodes(logging::Logger &logger, size_t threads) {
    logger.info() << "Resetting edge codes" << std::endl;
    this->fireResetEdgeCodes(logger, threads, *this);
    size_t code_size = 0;
//    Virtually suffix and merge edges get empty codes too (their start has no other outgoing edge to tell them from), as
//    GraphPath::resetEdgeCodes assumes.
    for(Edge &edge : edges()) {
        code_size += edge.edge_code.size();
        if(edge.isSuffix() || edge.isVirtualSuffix() || edge.isMerge())
            edge.edge_code = Sequence();
        else
            edge.edge_code = edge.truncSeq().Subseq(0, 1);
    }
    logger.info() << "Finished resetting edge codes. Reduced total code size from " << code_size << " to " << edgeCount() << std::endl;
}

Vertex &AssemblyGraph::addVertexPair(VertexData data, typename Vertex::id_type id) {
    VERIFY(id >= 0);
    Vertex *res;
#pragma omp critical
    {
        if(id == 0)
            id = maxVId + 1;
//        Vertices are appended, so to keep vertex_list sorted by id_cmp (|id| ascending, positive right before
//        negative) id must exceed every existing |id|, and the positive vertex must go in first.
        VERIFY(vertex_list.empty() || std::abs(vertex_list.back().getInnerId()) < id);
        VertexData rc_data = data.RC();
        Vertex &new_res = innerAddVertex(id, true, std::move(data));
        Vertex &rc = innerAddVertex(-id, false, std::move(rc_data));
        new_res.setRC(rc);
        res = &new_res;
    }
    this->fireAddVertex(*res);
    return *res;
}

bool AssemblyGraph::checkConsistency(std::ostream &os) {
    size_t errors = 0;
    auto report = [&]() -> std::ostream & {
        errors++;
        return os << "Graph inconsistency: ";
    };
    std::unordered_set<const Vertex *> present;
    for (const Vertex &v : vertex_list)
        present.insert(&v);
    const Vertex *prev = nullptr;
    for (const Vertex &v : vertex_list) {
        if (prev != nullptr && !id_cmp(prev->getInnerId(), v.getInnerId()))
            report() << "vertex list not sorted: " << prev->getInnerId() << " before " << v.getInnerId() << "\n";
        prev = &v;
    }
    for (const Vertex &v : vertex_list) {
        const Vertex::id_type vid = v.getInnerId();
        if (vid == 0 || std::abs(vid) > maxVId)
            report() << "vertex " << vid << " has id outside of (0, maxVId = " << maxVId << "]\n";
        if (v.rc_ == nullptr || present.count(v.rc_) == 0) {
            report() << "vertex " << vid << " has no rc in the graph\n";
            continue;
        }
        const Vertex &rc = *v.rc_;
        if (rc.rc_ != &v)
            report() << "vertex " << vid << ": rc().rc() is not the vertex itself\n";
        if (v.isCanonical() != (vid > 0))
            report() << "vertex " << vid << ": canonical flag disagrees with id sign\n";
        if (v != rc && rc.getInnerId() != -vid)
            report() << "vertex " << vid << ": rc id is " << rc.getInnerId() << " instead of " << -vid << "\n";
        if (v.mark_ != rc.mark_)
            report() << "vertex " << vid << ": marked flag differs from its rc\n";
        if (v.mark_ && !v.outgoing_.empty())
            report() << "vertex " << vid << " is marked but still has outgoing edges\n";
        if (v._outDeg != v.outgoing_.size())
            report() << "vertex " << vid << ": outDeg counter " << v._outDeg << " != " << v.outgoing_.size() << " outgoing edges\n";
        if (v.seq.empty()) {
            report() << "vertex " << vid << " has empty sequence\n";
        } else {
            if (rc.seq != !v.seq)
                report() << "vertex " << vid << ": rc sequence is not the reverse complement\n";
            if (v.isCanonical() != (v.seq <= !v.seq))
                report() << "vertex " << vid << ": canonical flag disagrees with sequence\n";
            if (v.seq == !v.seq && v != rc)
                report() << "vertex " << vid << " has palindromic sequence but is not self-rc\n";
        }

        bool has_suffix = false;
        std::array<bool, 4> first_nucls = {false, false, false, false};
        std::unordered_set<EdgeIdType> ids;
        for (const Edge &e : v.outgoing_) {
            const std::string name = "edge " + e.getInnerId().str();
            if (!ids.insert(e.getInnerId()).second)
                report() << name << ": duplicate edge id at vertex " << vid << "\n";
            if (e.start != &v)
                report() << name << ": start is not the vertex owning it\n";
            if (e.getInnerId().vid != vid)
                report() << name << ": id does not start with start vertex id " << vid << "\n";
            const int nucl = e.truncSeq().empty() ? 4 : int(e.truncSeq()[0]);
            if (e.getInnerId().eid <= 0 || e.getInnerId().eid % 10 != nucl)
                report() << name << ": last id digit is not the first nucleotide (" << nucl << ")\n";
            else if (e.getInnerId().eid / 10 > v.max_out_id[nucl])
                report() << name << ": id exceeds start vertex max_out_id\n";
            if (!e.corporeal)
                report() << name << " is not corporeal\n";
            if (e.isSuffix()) {
                has_suffix = true;
            } else {
                if (first_nucls[e.truncSeq()[0]])
                    report() << name << ": another outgoing edge of vertex " << vid << " has the same first nucleotide\n";
                first_nucls[e.truncSeq()[0]] = true;
            }
            if (e.finish == nullptr || present.count(e.finish) == 0) {
                report() << name << ": finish vertex is not in the graph\n";
                continue;
            }
            if (e._rc == nullptr) {
                report() << name << " has no rc\n";
                continue;
            }
            const Edge &erc = *e._rc;
            if (erc._rc != &e)
                report() << name << ": rc().rc() is not the edge itself\n";
            if (erc.start != e.finish->rc_ || erc.finish != v.rc_) {
                report() << name << ": rc edge does not go finish.rc() -> start.rc()\n";
            } else if (std::none_of(erc.start->outgoing_.begin(), erc.start->outgoing_.end(),
                                    [&erc](const Edge &x) { return &x == &erc; })) {
                report() << name << ": rc edge is not among outgoing edges of its start\n";
            }
            if (e.intCov() != erc.intCov())
                report() << name << ": coverage " << e.intCov() << " differs from rc coverage " << erc.intCov() << "\n";
            if (v.seq.empty() || e.finish->seq.empty())
                continue;
            const Sequence full = e.fullSeq();
            const size_t fsize = e.finish->size();
            if (full.size() < fsize || full.Subseq(full.size() - fsize) != e.finish->seq)
                report() << name << ": full sequence does not end with finish vertex sequence\n";
            if (&e <= &erc && erc.start != nullptr && !erc.start->seq.empty() && erc.fullSeq() != !full)
                report() << name << ": rc edge sequence is not the reverse complement\n";
        }
        if (has_suffix && v.outgoing_.size() > 1)
            report() << "vertex " << vid << " has an outgoing suffix edge and other outgoing edges\n";
    }
    if (errors > 0) {
        os << "Graph inconsistency: " << errors << " violation(s) found" << std::endl;
        return false;
    }
    return this->fireCheckConsistency(os);
}