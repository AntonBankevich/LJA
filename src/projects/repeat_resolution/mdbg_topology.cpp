//
// Created by Andrey Bzikadze on 11/10/21.
//

#include "mdbg_topology.hpp"
using namespace repeat_resolution;

// ---------- RRVertexProperty ----------

[[nodiscard]] bool RRVertexProperty::IsCanonical() const {
    return seq.IsCanonical();
}

void RRVertexProperty::IncLeft(MDBGSeq prefix) {
    seq.Prepend(std::move(prefix));
}

void RRVertexProperty::IncRight(MDBGSeq suffix) {
    seq.Append(std::move(suffix));
}

void RRVertexProperty::TrimLeft(const uint64_t size) { seq.TrimLeft(size); }

void RRVertexProperty::TrimRight(const uint64_t size) { seq.TrimRight(size); }

MDBGSeq RRVertexProperty::GetSeqPrefix(size_t len, int64_t shift) const {
    return seq.Substr(shift > 0 ? shift : 0, len);
}

MDBGSeq RRVertexProperty::GetSeqSuffix(size_t len, int64_t shift) const {
    uint64_t pos = seq.Size() - len;
    if (shift > 0) {
        pos -= shift;
    }
    return seq.Substr(pos, len);
}

[[nodiscard]] double RRVertexProperty::Cov() const {
    return seq.Cov();
}

std::ostream &repeat_resolution::operator<<(std::ostream &os,
                                            const RRVertexProperty &vertex) {
    os << vertex.size();
    // TODO show coverage
    // os << "size=" << vertex.size() << "\\n"
    //     << "cov=" << vertex.Cov();
    return os;
}

bool RRVertexProperty::operator==(const RRVertexProperty &rhs) const {
    return seq==rhs.seq and frozen==rhs.frozen;
}

// ---------- RREdgeProperty ----------

[[nodiscard]] int64_t RREdgeProperty::Size() const {
    if (not seq.Empty()) {
        VERIFY(size_==seq.Size())
    }
    return size_;
}

[[nodiscard]] bool RREdgeProperty::IsCanonical() const {
    return seq.IsCanonical();
}

void RREdgeProperty::Merge(RRVertexProperty vertex, RREdgeProperty rhs) {
    // in case current edge has negative length, there is an overlap b/w vertices
    int64_t vertex_size = vertex.size();
    if (size_ < 0) {
        vertex.TrimLeft((uint64_t) std::min((int64_t) vertex.size(), -size_));
    }
    size_ += (int64_t) vertex_size + rhs.size_;
    if (rhs.size_ < 0) {
        vertex.TrimRight((uint64_t) std::min((int64_t) vertex.size(),
                                             -rhs.size_));
    }
    seq.Append(std::move(vertex.seq));
    seq.Append(std::move(rhs.seq));
    if (size_ > 0) {
        VERIFY(seq.Size()==size_);
    }
    if (rhs.unique) {
        unique = true;
    }
}

MDBGSeq RREdgeProperty::ExtractSeqPrefix(const size_t len) {
    VERIFY(seq.Size() >= len);
    MDBGSeq prefix = seq.Substr(0, len);
    seq.TrimLeft(len);
    size_ -= len;
    return prefix;
}

MDBGSeq RREdgeProperty::ExtractSeqSuffix(const size_t len) {
    VERIFY(seq.Size() >= len);
    MDBGSeq suffix = seq.Substr(seq.Size() - len, len);
    seq.TrimRight(len);
    size_ -= len;
    return suffix;
}

void RREdgeProperty::ShortenWithEmptySeq(size_t len) {
    VERIFY(seq.Empty());
    size_ -= len;
}

[[nodiscard]] double RREdgeProperty::Cov() const {
    return seq.Cov();
}

bool repeat_resolution::operator==(const RREdgeProperty &lhs,
                                   const RREdgeProperty &rhs) {
    return lhs.Index()==rhs.Index();
}

bool repeat_resolution::operator!=(const RREdgeProperty &lhs,
                                   const RREdgeProperty &rhs) {
    return not(lhs==rhs);
}

RREdgeProperty repeat_resolution::Add(const RRVertexProperty &lhs,
                                      const RRVertexProperty &rhs,
                                      const uint64_t index) {
    VERIFY(lhs.size()==rhs.size());
    // can assign uniqueness more carefully if we pass left&right edges
    return RREdgeProperty(
        /*index=*/index,
        /*(inner)seq=*/{},
        /*size=*/-((int64_t) lhs.size()) + 1,
        /*unique=*/false);
}

// -------- SuccinctEdgeInfo -------

std::ostream &
repeat_resolution::operator<<(std::ostream &os,
                              const RREdgeProperty &edge_property) {
    os << "index=" << edge_property.Index() << "\\n"
       << "size=" << edge_property.Size() << "\\n"
       << "unique=" << edge_property.IsUnique();
       // TODO allow coverage
       // << "unique=" << edge_property.IsUnique() << "\\n"
       // << "cov=" << edge_property.Cov();
    return os;
}

bool repeat_resolution::operator==(const SuccinctEdgeInfo &lhs,
                                   const SuccinctEdgeInfo &rhs) {
    return lhs.start_ind==rhs.start_ind and lhs.end_ind==rhs.end_ind and
        lhs.edge==rhs.edge and lhs.unique==rhs.unique;
}

bool repeat_resolution::operator!=(const SuccinctEdgeInfo &lhs,
                                   const SuccinctEdgeInfo &rhs) {
    return !repeat_resolution::operator==(lhs, rhs);
}
