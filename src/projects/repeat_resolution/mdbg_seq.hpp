//
// Created by Andrey Bzikadze on 12/10/21.
//

#pragma once

#include "dbg/sparse_dbg.hpp"
#include "sequences/sequence.hpp"
#include <list>

namespace repeat_resolution {

// ---------- EdgeSegment ----------

struct EdgeSegment {
    const dbg::Edge *edge{nullptr};
    uint64_t start{0};
    uint64_t end{0};

    EdgeSegment(const dbg::Edge *edge, uint64_t start, uint64_t end);

    EdgeSegment(const EdgeSegment &) = default;
    EdgeSegment(EdgeSegment &&) = default;

    [[nodiscard]] uint64_t GetStK() const { return edge->start()->seq.size(); }
    [[nodiscard]] bool Empty() const { return start==end; }
    [[nodiscard]] bool LeftFull() const { return start==0; }
    [[nodiscard]] bool RightFull() const;
    void TrimLeft(uint64_t size);
    void TrimRight(uint64_t size);
    [[nodiscard]] uint64_t Size() const;
    [[nodiscard]] size_t EdgeSize() const;
    void ExtendRight(const EdgeSegment &segment);
    [[nodiscard]] EdgeSegment RC() const;
    [[nodiscard]] Sequence ToSequence() const;
    [[nodiscard]] double Cov() const;
    [[nodiscard]] bool operator==(const EdgeSegment &rhs) const;
};

// ---------- MDBGSeq ----------

class MDBGSeq {
    std::list<EdgeSegment> segms{};
    double cov{0};
    uint64_t size{0};

    static double CovListEdgeSegm(const std::list<EdgeSegment> &segms);
    static void Swap(MDBGSeq &lhs, MDBGSeq &rhs);

 public:
    MDBGSeq(const dbg::Edge *edge, uint64_t start, uint64_t end);
    explicit MDBGSeq(std::list<EdgeSegment> segms);
    MDBGSeq() = default;

    [[nodiscard]] Sequence ToSequence() const;
    [[nodiscard]] size_t Size() const;
    [[nodiscard]] size_t ContainerSize() const;
    [[nodiscard]] MDBGSeq RC() const;
    [[nodiscard]] bool IsCanonical() const;
    [[nodiscard]] bool Empty() const;
    void Append(MDBGSeq mdbg_seq);
    void Prepend(MDBGSeq mdbg_seq);
    void TrimLeft(uint64_t size_);
    void TrimRight(uint64_t size_);
    [[nodiscard]] MDBGSeq Substr(uint64_t pos, uint64_t len) const;
    [[nodiscard]] double Cov() const;
    [[nodiscard]] bool operator==(const MDBGSeq &rhs) const;
};

} // namespace repeat_resolution