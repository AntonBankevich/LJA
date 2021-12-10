//
// Created by Andrey Bzikadze on 12/10/21.
//

#include "mdbg_seq.hpp"

using namespace repeat_resolution;

// ---------- EdgeSegment ----------

EdgeSegment::EdgeSegment(const dbg::Edge *edge, const uint64_t start,
                         const uint64_t end)
    : edge{edge}, start{start}, end{end} {
  VERIFY(start <= end);
}

EdgeSegment::EdgeSegment(const dbg::Edge *edge)
    : EdgeSegment(edge, 0, edge->size()) {}

bool EdgeSegment::RightFull() const {
  VERIFY(edge != nullptr);
  return end == edge->size();
}

void EdgeSegment::TrimLeft(const uint64_t size) {
  VERIFY(Size() >= size);
  start += size;
}

void EdgeSegment::TrimRight(const uint64_t size) {
  VERIFY(Size() >= size);
  end -= size;
}

uint64_t EdgeSegment::Size() const {
  VERIFY(start <= end);
  return end - start;
}

size_t EdgeSegment::EdgeSize() const {
  VERIFY(edge != nullptr);
  return edge->size();
}

void EdgeSegment::ExtendRight(const EdgeSegment &segment) {
  VERIFY(edge == segment.edge);
  VERIFY(end == segment.start);
  end = segment.end;
}

EdgeSegment EdgeSegment::RC() const {
  VERIFY(edge != nullptr);
  // TODO check that this is safe
  return EdgeSegment(&edge->rc(), edge->size() - end, edge->size() - start);
}

Sequence EdgeSegment::ToSequence() const {
  VERIFY(edge != nullptr);
  return edge->seq.Subseq(start, end);
}

bool EdgeSegment::operator==(const EdgeSegment &rhs) const {
  return edge == rhs.edge and start == rhs.start and end == rhs.end;
}

// ---------- MDBGSeq2 ----------

[[nodiscard]] Sequence MDBGSeq2::ToSequence() const {
  std::vector<Sequence> sec_vec;
  for (const EdgeSegment segm : segms) {
    sec_vec.emplace_back(segm.ToSequence());
  }
  return Sequence::Concat(sec_vec);
}

[[nodiscard]] size_t MDBGSeq2::Size() const {
  size_t size{0};
  for (const EdgeSegment segm : segms) {
    size += segm.Size();
  }
  return size;
}

[[nodiscard]] size_t MDBGSeq2::ContainerSize() const {
  return segms.size();
}

[[nodiscard]] MDBGSeq2 MDBGSeq2::RC() const {
  std::list<EdgeSegment> segms_rc;
  for (const EdgeSegment segm : segms) {
    segms_rc.emplace_front(segm.RC());
  }
  return MDBGSeq2(std::move(segms_rc));
}

[[nodiscard]] bool MDBGSeq2::IsCanonical() const {
  // TODO maybe there is a faster way?
  return ToSequence() <= RC().ToSequence();
}

[[nodiscard]] bool MDBGSeq2::Empty() const {
  return segms.empty();
}

void MDBGSeq2::Append(MDBGSeq2 mdbg_seq) {
  if (mdbg_seq.Empty()) {
    return;
  }
  if (Empty()) {
    segms = std::move(mdbg_seq.segms);
  }
  EdgeSegment &back = segms.back();
  EdgeSegment &front = mdbg_seq.segms.front();

  if (back.edge == front.edge) {
    back.ExtendRight(front);
    mdbg_seq.segms.pop_front();
  }

  segms.splice(segms.end(), std::move(mdbg_seq.segms));
}

void MDBGSeq2::Prepend(MDBGSeq2 mdbg_seq) {
  MDBGSeq2 temp;
  temp.segms = std::move(segms);
  mdbg_seq.Append(std::move(temp));
  segms = std::move(mdbg_seq.segms);
}

void MDBGSeq2::TrimLeft(uint64_t size) {
  VERIFY(size <= Size());
  while (size > 0) {
    EdgeSegment &front = segms.front();
    if (front.Size() <= size) {
      size -= front.Size();
      segms.pop_front();
    } else {
      front.TrimLeft(size);
      size = 0;
    }
  }
}

void MDBGSeq2::TrimRight(uint64_t size) {
  VERIFY(size <= Size());
  while (size > 0) {
    EdgeSegment &back = segms.back();
    if (back.Size() <= size) {
      size -= back.Size();
      segms.pop_back();
    } else {
      back.TrimRight(size);
      size = 0;
    }
  }
}

[[nodiscard]] MDBGSeq2 MDBGSeq2::Substr(uint64_t pos, uint64_t len) const {
  VERIFY(pos + len <= Size());
  auto left = segms.begin();
  while(left->Size() <= pos) {
    pos -= left->Size();
    ++left;
  }
  auto right = left;
  while(right->Size() <= len) {
    len -= right->Size();
    ++right;
  }
  if (left == right) {
    return MDBGSeq2({EdgeSegment(left->edge, pos, pos + len)});
  }

  std::list<EdgeSegment> res;
  res.emplace_back(left->edge, pos, left->edge->size());
  ++left;
  for(; left != right; ++left) {
    res.emplace_back(left->edge);
  }
  res.emplace_back(right->edge, 0, pos + len);
  return MDBGSeq2(res);
}

[[nodiscard]] bool MDBGSeq2::operator==(const MDBGSeq2 &rhs) const {
  return segms == rhs.segms;
}
