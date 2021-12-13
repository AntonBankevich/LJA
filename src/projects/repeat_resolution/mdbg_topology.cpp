//
// Created by Andrey Bzikadze on 11/10/21.
//

#include "mdbg_topology.hpp"
using namespace repeat_resolution;

// ---------- MDBGSeq ----------

/*
constexpr char MDBGSeq::CharCompl(const char c) {
  constexpr std::string_view bases = "ATCG";
  return bases[bases.find(c) ^ 1];
}

MDBGSeq::MDBGSeq(std::list<char> seq) : seq{std::move(seq)} {}
MDBGSeq::MDBGSeq(std::string str) {
  std::move(str.begin(), str.end(), std::back_inserter(seq));
}

MDBGSeq::MDBGSeq(const Sequence &seq) : MDBGSeq(std::move(seq.str())) {}

std::string MDBGSeq::ToSequence() const {
  std::string str;
  std::move(seq.begin(), seq.end(), std::back_inserter(str));
  return str;
}

size_t MDBGSeq::size() const { return seq.size(); }

MDBGSeq MDBGSeq::RC() const {
  std::list<char> rc;
  std::transform(seq.crbegin(), seq.crend(), std::back_inserter(rc), CharCompl);
  return MDBGSeq(std::move(rc));
}

bool MDBGSeq::IsCanonical() const {
  for (auto [it, rit] = std::make_pair(seq.cbegin(), seq.crbegin());
       it != seq.cend(); ++it, ++rit) {
    char f = *it;
    char r = CharCompl(*rit);
    if (f < r) {
      return true;
    }
    if (r < f) {
      return false;
    }
  }
  return true;
}

bool MDBGSeq::Empty() const { return seq.empty(); }

void MDBGSeq::Append(MDBGSeq mdbg_seq) {
  seq.splice(seq.end(), std::move(mdbg_seq.seq));
}

void MDBGSeq::Prepend(MDBGSeq mdbg_seq) {
  seq.splice(seq.begin(), std::move(mdbg_seq.seq));
}

void MDBGSeq::TrimLeft(uint64_t size) {
  for (uint64_t i = 0; i < size; ++i) {
    seq.pop_front();
  }
}

void MDBGSeq::TrimRight(uint64_t size) {
  for (uint64_t i = 0; i < size; ++i) {
    seq.pop_back();
  }
}

MDBGSeq MDBGSeq::Substr(const uint64_t pos, const uint64_t len) const {
  VERIFY(pos + len <= seq.size());
  auto it = [this, &pos, &len] {
    if (pos < seq.size() - pos) {
      auto it = seq.cbegin();
      std::advance(it, pos);
      return it;
    }
    auto it = seq.crbegin();
    std::advance(it, seq.size() - pos);
    return it.base();
  }();
  std::list<char> substr_list;
  for (uint64_t cnt = 0; cnt < len; ++cnt, ++it) {
    substr_list.emplace_back(*it);
  }
  return MDBGSeq(substr_list);
}

bool MDBGSeq::operator==(const MDBGSeq &rhs) const { return seq == rhs.seq; }
*/

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

std::ostream &repeat_resolution::operator<<(std::ostream &os,
                                            const RRVertexProperty &vertex) {
  os << vertex.size();
  return os;
}

bool RRVertexProperty::operator==(const RRVertexProperty &rhs) const {
  return seq == rhs.seq and frozen == rhs.frozen;
}

// ---------- RREdgeProperty ----------

[[nodiscard]] int64_t RREdgeProperty::size() const {
  if (not seq.Empty()) {
    VERIFY(size_ == seq.Size())
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
    vertex.TrimLeft((uint64_t)std::min((int64_t)vertex.size(), -size_));
  }
  size_ += (int64_t)vertex_size + rhs.size_;
  if (rhs.size_ < 0) {
    vertex.TrimRight((uint64_t)std::min((int64_t)vertex.size(), -rhs.size_));
  }
  seq.Append(std::move(vertex.seq));
  seq.Append(std::move(rhs.seq));
  if (size_ > 0) {
    VERIFY(seq.Size() == size_);
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

bool repeat_resolution::operator==(const RREdgeProperty &lhs,
                                   const RREdgeProperty &rhs) {
  return lhs.Index() == rhs.Index();
}

bool repeat_resolution::operator!=(const RREdgeProperty &lhs,
                                   const RREdgeProperty &rhs) {
  return not(lhs == rhs);
}

RREdgeProperty repeat_resolution::Add(const RRVertexProperty &lhs,
                                      const RRVertexProperty &rhs,
                                      const uint64_t index) {
  VERIFY(lhs.size() == rhs.size());
  // can assign uniqueness more carefully if we pass left&right edges
  return RREdgeProperty(/*index=*/index, /*(inner)seq=*/{},
                        /*size=*/-((int64_t)lhs.size()) + 1, /*unique=*/false);
}

// -------- SuccinctEdgeInfo -------

std::ostream &
repeat_resolution::operator<<(std::ostream &os,
                              const RREdgeProperty &edge_property) {
  os << "index=" << edge_property.Index() << "\\n"
     << "size=" << edge_property.size() << "\\n"
     << "unique=" << edge_property.IsUnique();
  return os;
}

bool repeat_resolution::operator==(const SuccinctEdgeInfo &lhs,
                                   const SuccinctEdgeInfo &rhs) {
  return lhs.start_ind == rhs.start_ind and lhs.end_ind == rhs.end_ind and
         lhs.edge == rhs.edge and lhs.unique == rhs.unique;
}

bool repeat_resolution::operator!=(const SuccinctEdgeInfo &lhs,
                                   const SuccinctEdgeInfo &rhs) {
  return !repeat_resolution::operator==(lhs, rhs);
}
