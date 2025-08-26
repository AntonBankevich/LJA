#include "assembly_graph_base.hpp"

#include <utility>

using namespace ag;
Edge::Edge() : id(0, 0), start(nullptr), finish(nullptr), seq(), _rc(nullptr) {
//            TODO: Remove this!!! It exists only for Andreys code compilation but that code should be purged
}

Edge::Edge(Edge::id_type id, Vertex &_start, Vertex &_end, Sequence _seq, EdgeData data) :
        EdgeData(std::move(data)), id(id), start(&_start), finish(&_end), seq(std::move(_seq)), _rc(nullptr) {
}

bool Edge::operator<(const Edge &other) const {
    if(this == &other)
        return false;
    if(start != other.start)
        return *start < *other.start;
    return this->truncSeq() < other.truncSeq();
}

bool Edge::operator>(const Edge &other) const {
    if(this == &other)
        return false;
    if(start != other.start)
        return *start > *other.start;
    return other.truncSeq() < truncSeq();
}

Sequence Edge::fullSubseq(size_t from, size_t to) const {
    VERIFY(start->size() > 0);
    VERIFY(from <= start->size() + truncSize());
    VERIFY(to <= start->size() + truncSize());
    if (from >= start->size()) {
        return truncSeq().Subseq(from - start->size(), to - start ->size());
    } else if(to <= start->size()) {
        return start->getSeq().Subseq(from, to);
    } else {
        return getStart().getSeq().Subseq(from) + truncSeq().Subseq(0, to - start->size());
    }
}

Sequence Edge::suffix(size_t pos) const {
    VERIFY(pos <= truncSeq().size());
    size_t k = start->getSeq().size();
    if (pos >= k)
        return truncSeq().Subseq(pos - k, truncSeq().size());
    else {
        return getStart().getSeq().Subseq(pos) + truncSeq().Subseq(0, truncSeq().size());
    }
}

Sequence Edge::getSeq() const {
    return start->getSeq() + seq;
}

Sequence Edge::fullSeq() const {return getSeq();}

size_t Edge::getStartSize() const {return start->size();}

void Edge::DeleteEdge(Edge &edge) {
    Locker<VertexId> locker = Locker<VertexId>::FromVector({edge.start->getId(), edge.finish->getId()});
    DeleteEdgeLockFree(edge);
}

void Edge::DeleteEdgeLockFree(Edge &edge) {
    Vertex &start = *edge.start;
    Vertex &rcstart = *edge.finish;
    Edge &rcedge = edge.rc();
    if(edge != rcedge) {
        rcstart.innerRemoveEdge(rcedge);
    }
    start.innerRemoveEdge(edge);
}

size_t Edge::innerSize() const {
    size_t full = fullSize();
    size_t vsize = getStart().size() + getFinish().size();
    return full >= vsize ? full - vsize : 0;
}

size_t Edge::overlapSize() const {
    size_t full = fullSize();
    size_t vsize = getStart().size() + getFinish().size();
    return full >= vsize ? 0 : vsize - full;
}

bool Edge::isOuter() const { return getStart().outDeg() > 1 && getFinish().inDeg() > 1; }

bool Edge::isInner() const { return getStart().outDeg() == 1 && getFinish().inDeg() == 1; }

void Edge::incCov(int64_t delta) {
#pragma omp atomic
    cov += delta;
    VERIFY(cov < size_t(-1) >> 2)
}

void Vertex::checkConsistency() const {
    VERIFY(!seq.empty());
    VERIFY(isCanonical() || rc().isCanonical());
    VERIFY(*this == rc() || isCanonical() != rc().isCanonical());
    VERIFY(isCanonical() == seq <= !seq);
    for (const Edge &edge : outgoing_) {
        VERIFY(edge.isCanonical() || edge.rc().isCanonical());
//            VERIFY(edge.isCanonical() == edge.getSeq() <= !edge.getSeq());
        VERIFY(edge == edge.rc() || (edge.isCanonical() != edge.rc().isCanonical()));
        VERIFY(edge.rc().getFinish() == this->rc());
        VERIFY(std::find(edge.getFinish().rc().begin(), edge.getFinish().rc().end(), edge.rc()) != edge.getFinish().rc().end());
        VERIFY(edge.intCov() == edge.rc().intCov());
    }
}

void Vertex::setSeq(Sequence _seq) {
    lock();
    VERIFY(isCanonical() == _seq.isCanonical());
    if (getSeq().empty()) {
        seq = std::move(_seq);
        Sequence rc_seq = seq.rc();
        unlock();
        if(rc_ != nullptr) {
            rc_->lock();
            if(rc_->getSeq().empty()) {
                rc_->seq = std::move(rc_seq);
            }
            rc_->unlock();
        }
    } else {
        unlock();
    }
}

bool Vertex::hasOutgoing(unsigned char c) const {
    for (const Edge &edge : outgoing_) {
        if (edge.truncSeq()[0] == c) {
            return true;
        }
    }
    return false;
}

void Vertex::sortOutgoing() {
    outgoing_.sort();
//    std::sort(outgoing_.begin(), outgoing_.end());
}

bool Vertex::isJunction() const {
    return outDeg() != 1 || inDeg() != 1;
}

bool Vertex::hasOutgoingSuffix() const {
    for(Edge &edge : *this) {
        if(edge.truncSize() == 0)
            return true;
    }
    return false;
}

bool Vertex::isPalindrome() const {return *this == rc();}

Edge &Vertex::getOutgoing(unsigned char c) const {
    size_t cnt = 0;
    for (Edge &edge : outgoing_) {
        if (edge.corporeal && (edge.isSuffix() || edge.getCode()[0] == c)) {
            cnt++;
        }
    }
    VERIFY(cnt <= 1);
    for (Edge &edge : outgoing_) {
        if (edge.corporeal && (edge.isSuffix() || edge.getCode()[0] == c)) {
            return edge;
        }
    }
    std::cout << "Outgoing edge not found" << std::endl;
    std::cout << "Vertex: " << *this << " " << mark_ << " Missing outgoing nucleotide: " << nucl(c) << std::endl;
    std::cout << "Vertex seq: " << getSeq() << std::endl;
    for (const Edge &edge : outgoing_) {
        std::cout << "Outgoing code: " << edge.getCode() << std::endl;
    }
    VERIFY(false);
    return outgoing_.front();
}

Edge &Vertex::innerAddEdge(Vertex &end, const Sequence &tseq, EdgeData data, EdgeIdType eid) {
    if (!eid.valid()) {
        int code = tseq.empty() ? 4 : tseq[0];
        eid = {id, (max_out_id[code] + 1) * 10 + code};
    }
    VERIFY(eid.vid == id);
    updateMaxOutId(eid.eid);
    outgoing_.emplace_back(eid, *this, end, tseq, std::move(data));
    Edge &edge = outgoing_.back();
    _outDeg++;
    return edge;
}

std::vector<EdgePosition> EdgePosition::step() const {
    if (pos == edge->truncSize()) {
        std::vector<EdgePosition> res;
        auto &v = edge->getFinish();
        for (Edge &next : v) {
            res.emplace_back(next, 1);
        }
        return std::move(res);
    } else {
        return {{*edge, pos + 1}};
    }
}

std::ostream &std::operator<<(std::ostream &os, ag::EdgeSaveLabel val) {
    return os << val.fId << "_" << val.rcId;
}

std::ostream &ag::operator<<(std::ostream &os, const Vertex &vertex) {
    return os << vertex.getId();
}
std::ostream &ag::operator<<(std::ostream &os, const Edge &edge) {
    return os << edge.getId();
}

size_t std::hash<ag::EdgeIdType>::operator()(const ag::EdgeIdType &x) const {
    return std::hash<int>()(x.vid * 12343251) ^ std::hash<int>()(x.eid * 1294835);
}
