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

Vertex &ag::AssemblyGraph::addSPGVertex(Sequence seq, bool cyclic, bool inf_left, bool inf_right, int id) {
    Vertex &res =  addVertex(std::move(seq), ag::VertexData::SPGData(cyclic, inf_left, inf_right), id);
    return res;
}

Vertex &ag::AssemblyGraph::mergePath(const GraphPath &path) {
    VERIFY(path.startClosed() && path.endClosed());
    PathPosition p = path.firstPosition();
    while (p != path.lastPosition() && p.nextEdge().isSuffix())
        ++p;
    bool updown = true;
    for (PathPosition p1 = p; p1 != path.lastPosition(); ++p1) {
        if (!p1.nextEdge().isSuffix()) {
            updown = false;
            break;
        }
    }
    if (updown) {
        GraphPath left_path = path.subPath(path.firstPosition(), p);
        GraphPath right_path = path.subPath(p, path.lastPosition());
        Vertex & res = p.getVertex();
        if (left_path.calculateSize() > 1) {
            mergePathToEdge(left_path);
        }
        if (right_path.calculateSize() > 1) {
            mergePathToEdge(right_path);
        }
        return res;
    } else {
        Sequence seq = path.Seq();
        Vertex &res = addSPGVertex(seq, false, false, false);
        Edge &inc_edge = addSPEdgeLockFree(path.getStart(), res);
        inc_edge.setCorporeal(false);
        inc_edge.rc().setCorporeal(false);
        Edge &out_edge = res == res.rc() ? inc_edge.rc() : addSPEdgeLockFree(res, path.getFinish());
        out_edge.setCorporeal(false);
        out_edge.rc().setCorporeal(false);
        fireMergePath(path.asRAPath(), res);
        isolateAndMark(path.innerVertices().begin(), path.innerVertices().end());
        inc_edge.setCorporeal(true);
        inc_edge.rc().setCorporeal(true);
        out_edge.setCorporeal(true);
        out_edge.rc().setCorporeal(true);
        return res;
    }
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
    VERIFY(path.calculateSize() > 1 && path.getStart() == path.getFinish());
    Sequence new_seq = path.truncSeq();
    size_t shift = SmallestShift(new_seq);;
    Sequence new_seq_shifted = new_seq.Subseq(shift) + new_seq.Subseq(0, shift);
    Vertex &res = addSPGVertex(new_seq_shifted, true, false, false);
//    Edge &edge = addEdgeLockFree(res, res, seq, seq);
    fireMergeLoop(path, res);
    isolateAndMark(path.innerVertices().begin(), path.innerVertices().end());
    isolateAndMark(path.getStart());
    return res;
}

Edge &ag::AssemblyGraph::addSPEdgeLockFree(Vertex &start, Vertex &end,
                                                ag::Edge::id_type eid,
                                                ag::Edge::id_type rcid) {
    Sequence tseq = end.getSeq().Subseq(std::min(start.size(), end.size()));
    Sequence rctseq = start.rc().getSeq().Subseq(std::min(start.size(), end.size()));
    return addEdgeLockFree(start, end, tseq, rctseq, EdgeData(), eid, rcid);
}

Edge & ag::AssemblyGraph::addSPEdge(Vertex &start, Vertex &end, ag::Edge::id_type eid,
                                         ag::Edge::id_type rcid) {
    Sequence tseq = end.getSeq().Subseq(std::min(start.size(), end.size()));
    Sequence rctseq = start.rc().getSeq().Subseq(std::min(start.size(), end.size()));
    return addEdge(start, end, tseq, rctseq, EdgeData(), eid, rcid);
}

Edge &AssemblyGraph::addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeIdType eid,
                             EdgeIdType rcid) {
    return addEdge(start, end, tseq, rctseq, EdgeData(), eid, rcid);
}

Vertex &AssemblyGraph::edgeToSupreVertex(Edge &edge) {
    Locker<VertexId> locker = Locker<VertexId>::FromVector({edge.getStart().getId(), edge.getFinish().rc().getId()});
    VERIFY(!edge.isSuffix());
    VERIFY(!edge.isPrefix());
    Sequence seq = edge.fullSeq();
    Vertex &res = addVertex(seq);
    addEdgeLockFree(res, edge.getFinish(), Sequence(), edge.rc().truncSeq());
    if(edge != edge.rc())
        addEdgeLockFree(res.rc(), edge.getStart().rc(), Sequence(), edge.truncSeq());
    res.front().setCorporeal(false);
    res.front().rc().setCorporeal(false);
    res.rc().front().setCorporeal(false);
    res.rc().front().rc().setCorporeal(false);
    fireEdgeToSupreVertex(res, edge);
    removeEdgeLockFree(edge);
    res.front().setCorporeal(true);
    res.front().rc().setCorporeal(true);
    res.rc().front().setCorporeal(true);
    res.rc().front().rc().setCorporeal(true);
    return res;
}

VertexResolutionResult
AssemblyGraph::resolveVertex(Vertex &core, const VertexResolutionPlan &resolution) {
    VERIFY(core.isCore() && core.inDeg() > 0 && core.outDeg() > 0);
    VertexResolutionResult result(core);
    for(const InOutEdgePair &p : resolution.connectionsUnique()) {
        VERIFY(p.incoming().getFinish() == core);
        if(p.middle() != core)
            continue;
        VERIFY(p.outgoing().getStart() == core);
        VERIFY(p.incoming().isSuffix());
        VERIFY(p.outgoing().isPrefix());
        Sequence seq = p.getSeq();
        Vertex &newv = addVertex(seq);
        result.add(newv, p);
        Edge &out = addEdgeLockFree(p.incoming().getStart(), newv, newv.getSeq());
        out.setCorporeal(false);
        out.rc().setCorporeal(false);
        if(newv != newv.rc()) {
            Edge &inc = addEdgeLockFree(newv, p.outgoing().getFinish(), newv.getSeq());
            inc.setCorporeal(false);
            inc.rc().setCorporeal(false);
        }
    }
    fireResolveVertex(core, result);
    isolateAndMark(core);
    for(Vertex &v : result.newVertices()) {
        v.front().setCorporeal(true);
        v.front().rc().setCorporeal(true);
        v.rc().front().setCorporeal(true);
        v.rc().front().rc().setCorporeal(true);
    }
    return std::move(result);
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
    if(it->corporeal)
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
    typename Vertex::id_type id = maxVId + 1;
    Vertex &res = innerAddVertex(id, true, std::move(data));
    res.setRC(res);
    this->fireAddVertex(res);
    return res;
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
            this->isolateAndMark(v);
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

Vertex &AssemblyGraph::addVertex(const Sequence &seq, const VertexData &data, typename Vertex::id_type id) {
    if (id == 0) {
        if(seq <= !seq)
            id = maxVId + 1;
        else
            id = -maxVId - 1;
    }
    maxVId = std::max(std::abs(id), maxVId);
    Vertex & res = innerAddVertex(id, seq, data);
    Vertex & rc = seq == !seq ? res : innerAddVertex(-id, !seq, data.RC());
    res.setRC(rc);
    res.setSeq(seq);
    this->fireAddVertex(res);
    return res;
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
            isolateAndMark(v);
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

Vertex &AssemblyGraph::innerAddVertex(typename Vertex::id_type id, Sequence seq, VertexData data) {
    VERIFY(seq.isCanonical() == (id > 0));
    maxVId = std::max(std::abs(id), maxVId);
    vertex_list.emplace_back(id, std::move(seq), std::move(data));
    return vertex_list.back();
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

    return addEdgeLockFree(start, end, tseq, rctseq, data, eid, rcid);
}

void AssemblyGraph::removeEdgeLockFree(Edge &edge) {
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

void AssemblyGraph::isolateAndMark(Vertex &vertex) {
    VERIFY(!vertex.marked());
    for(Vertex &v : ThisAndRC(vertex)) {
        while (v.outDeg() != 0) {
            removeEdgeLockFree(v.front());
        }
        v.mark();
    }
}

//    TODO: when TAGraphPath is thoroughly damned, put it back here as a parameter.
Edge &AssemblyGraph::mergePathToEdge(const GraphPath &path) {
    VERIFY(!path.empty());
    VERIFY(path.endClosed() && path.startClosed());
    Locker<VertexId> locker = Locker<VertexId>::FromVector({path.getStart().getId(), path.getRCStart().getId()});
    for(Vertex &v : path.innerVertices()) {
        VERIFY(!v.marked());
    }
    VERIFY(!path.isSingleton());
    VERIFY(path.getStart() == path.getFinish() || path.getStart() == path.getFinish().rc() || (path.getStart().isJunction() && path.getFinish().isJunction()));
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
            isolateAndMark(*vid);
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
        isolateAndMark(rightEdge.getStart());
    isolateAndMark(leftEdge.getFinish());
    new_edge.setCorporeal(true);
    new_edge.rc().setCorporeal(true);
    return new_edge;
}

//    TODO: Make it run in parallel
void AssemblyGraph::resetEdgeCodes(logging::Logger &logger, size_t threads) {
    logger.info() << "Resetting edge codes" << std::endl;
    this->fireResetEdgeCodes(logger, threads, *this);
    size_t code_size = 0;
    for(Edge &edge : edges()) {
        code_size += edge.edge_code.size();
        if(!edge.isSuffix())
            edge.edge_code = edge.truncSeq().Subseq(0, 1);
    }
    logger.info() << "Finished resetting edge codes. Reduced total code size from " << code_size << " to " << edgeCount() << std::endl;
}

Vertex &AssemblyGraph::addVertexPair(VertexData data, typename Vertex::id_type id) {
    VERIFY(id >= 0);
    if(id == 0)
        id = maxVId + 1;
    Vertex &rc = innerAddVertex(-id, false, data.RC());
    Vertex &res = innerAddVertex(id, true, std::move(data));
    res.setRC(rc);
    this->fireAddVertex(res);
    return res;
}