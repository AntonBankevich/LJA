//
// Created by Andrey Bzikadze on 12/10/21.
//


#pragma once

#include <list>
#include "dbg/sparse_dbg.hpp"
#include "sequences/sequence.hpp"

namespace repeat_resolution {

// ---------- EdgeSegment ----------

struct EdgeSegment {
  const dbg::Edge *edge{nullptr};
  uint64_t start{0};
  uint64_t end{0};

  EdgeSegment(const dbg::Edge *edge, uint64_t start, uint64_t end);
  explicit EdgeSegment(const dbg::Edge *edge);

  EdgeSegment(const EdgeSegment &) = default;
  EdgeSegment(EdgeSegment &&) = default;

  [[nodiscard]] bool Empty() const { return start == end; }
  [[nodiscard]] bool LeftFull() const { return start == 0; }
  [[nodiscard]] bool RightFull() const;
  void TrimLeft(uint64_t size);
  void TrimRight(uint64_t size);
  [[nodiscard]] uint64_t Size() const;
  [[nodiscard]] size_t EdgeSize() const;
  void ExtendRight(const EdgeSegment &segment);
  [[nodiscard]] EdgeSegment RC() const;
  [[nodiscard]] Sequence ToSequence() const;

  [[nodiscard]] bool operator==(const EdgeSegment &rhs) const;
};

// ---------- MDBGSeq2 ----------

class MDBGSeq2 {
  std::list<EdgeSegment> segms;

public:
  MDBGSeq2(const dbg::Edge *edge, const uint64_t start, const uint64_t end) {
    segms.emplace_back(edge, start, end);
  }

  explicit MDBGSeq2(const dbg::Edge *edge) : MDBGSeq2(edge, 0, edge->size()) {}
  explicit MDBGSeq2(std::list<EdgeSegment> segms) : segms{std::move(segms)} {}
  MDBGSeq2() = default;

  [[nodiscard]] Sequence ToSequence() const;
  [[nodiscard]] size_t Size() const;
  [[nodiscard]] size_t ContainerSize() const;
  [[nodiscard]] MDBGSeq2 RC() const;
  [[nodiscard]] bool IsCanonical() const;
  [[nodiscard]] bool Empty() const;
  void Append(MDBGSeq2 mdbg_seq);
  void Prepend(MDBGSeq2 mdbg_seq);
  void TrimLeft(uint64_t size);
  void TrimRight(uint64_t size);
  [[nodiscard]] MDBGSeq2 Substr(uint64_t pos, uint64_t len) const;
  [[nodiscard]] bool operator==(const MDBGSeq2 &rhs) const;
};

} // namespace repeat_resolution