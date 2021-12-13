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
  VERIFY(edge != nullptr and segment.edge != nullptr)
  VERIFY(*edge == *segment.edge);
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

// ---------- MDBGSeq ----------

[[nodiscard]] Sequence MDBGSeq::ToSequence() const {
  std::vector<Sequence> sec_vec;
  for (const EdgeSegment segm : segms) {
    sec_vec.emplace_back(segm.ToSequence());
  }
  return Sequence::Concat(sec_vec);
}

[[nodiscard]] size_t MDBGSeq::Size() const {
  size_t size{0};
  for (const EdgeSegment segm : segms) {
    size += segm.Size();
  }
  return size;
}

[[nodiscard]] size_t MDBGSeq::ContainerSize() const { return segms.size(); }

[[nodiscard]] MDBGSeq MDBGSeq::RC() const {
  std::list<EdgeSegment> segms_rc;
  for (const EdgeSegment segm : segms) {
    segms_rc.emplace_front(segm.RC());
  }
  return MDBGSeq(std::move(segms_rc));
}

[[nodiscard]] bool MDBGSeq::IsCanonical() const {
  // TODO maybe there is a faster way?
  return ToSequence() <= RC().ToSequence();
}

[[nodiscard]] bool MDBGSeq::Empty() const { return segms.empty(); }

void MDBGSeq::Append(MDBGSeq mdbg_seq) {
  if (mdbg_seq.Empty()) {
    return;
  }
  if (Empty()) {
    segms = std::move(mdbg_seq.segms);
  }
  EdgeSegment &back = segms.back();
  EdgeSegment &front = mdbg_seq.segms.front();

  if (back.edge == front.edge and back.end == front.start) {
    back.ExtendRight(front);
    mdbg_seq.segms.pop_front();
  }

  segms.splice(segms.end(), std::move(mdbg_seq.segms));
}

void MDBGSeq::Prepend(MDBGSeq mdbg_seq) {
  MDBGSeq temp;
  temp.segms = std::move(segms);
  mdbg_seq.Append(std::move(temp));
  segms = std::move(mdbg_seq.segms);
}

void MDBGSeq::TrimLeft(uint64_t size) {
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

void MDBGSeq::TrimRight(uint64_t size) {
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

[[nodiscard]] MDBGSeq MDBGSeq::Substr(uint64_t pos, const uint64_t len) const {
  VERIFY(pos + len <= Size());
  if (len == 0 or Empty()) {
    return MDBGSeq();
  }

  auto left = segms.begin();
  while (left->Size() <= pos) {
    pos -= left->Size();
    ++left;
    VERIFY(left != segms.end());
  }
  auto right = left;
  uint64_t end = pos + len;
  while (right->Size() < end) {
    end -= right->Size();
    ++right;
    VERIFY(right != segms.end());
  }
  if (left == right) {
    return MDBGSeq(
        {EdgeSegment(left->edge, left->start + pos, left->start + end)});
  }

  std::list<EdgeSegment> res;
  res.emplace_back(left->edge, left->start + pos, left->end);
  for (++left; left != right; ++left) {
    res.emplace_back(*left);
  }
  res.emplace_back(right->edge, right->start, right->start + end);
  MDBGSeq subseq(std::move(res));
  VERIFY(subseq.Size() == len);
  return subseq;
}

[[nodiscard]] bool MDBGSeq::operator==(const MDBGSeq &rhs) const {
  return segms == rhs.segms;
}
