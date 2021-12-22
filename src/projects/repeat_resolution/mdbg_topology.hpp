//
// Created by Andrey Bzikadze on 10/14/21.
//

#pragma once

#include "common/verify.hpp"
#include "dbg/sparse_dbg.hpp"
#include "sequences/sequence.hpp"
#include <functional>
#include <graphlite/graphlite.hpp>
#include <graphlite/serialize.hpp>
#include "mdbg_seq.hpp"

namespace repeat_resolution {

// ---------- RRVertexProperty ----------

using RRVertexType = uint64_t;

class RRVertexProperty {
    MDBGSeq seq;
    bool frozen{false};

 public:
    friend class RREdgeProperty;
    RRVertexProperty(MDBGSeq seq, const bool frozen)
        : seq{std::move(seq)}, frozen{frozen} {}

    RRVertexProperty(const RRVertexProperty &) = delete;
    RRVertexProperty(RRVertexProperty &&) = default;
    RRVertexProperty &operator=(const RRVertexProperty &) = delete;
    RRVertexProperty &operator=(RRVertexProperty &&) = default;

    [[nodiscard]] bool IsCanonical() const;
    [[nodiscard]] uint64_t size() const { return seq.Size(); }
    [[nodiscard]] const MDBGSeq &Seq() const { return seq; }
    [[nodiscard]] bool IsFrozen() const { return frozen; }

    void Freeze() { frozen = true; }

    void IncLeft(MDBGSeq prefix);
    void IncRight(MDBGSeq suffix);
    void TrimLeft(uint64_t inc = 1);
    void TrimRight(uint64_t inc = 1);

    [[nodiscard]] MDBGSeq GetSeqPrefix(size_t len, int64_t shift = 0) const;
    [[nodiscard]] MDBGSeq GetSeqSuffix(size_t len, int64_t shift = 0) const;

    [[nodiscard]] double Cov() const;

    [[nodiscard]] bool operator==(const RRVertexProperty &rhs) const;
};

std::ostream &operator<<(std::ostream &os, const RRVertexProperty &vertex);

// ---------- RREdgeProperty ----------

using RREdgeIndexType = uint64_t;

class RREdgeProperty {
    RREdgeIndexType index{0};
    MDBGSeq seq;
    int64_t size_{0};
    bool unique{false};

 public:
    RREdgeProperty(const RREdgeIndexType index,
                   MDBGSeq seq,
                   const int64_t size_,
                   bool unique)
        : index{index},
          seq{std::move(seq)},
          size_{size_},
          unique{unique} {}

    RREdgeProperty(const RREdgeProperty &) = delete;
    RREdgeProperty(RREdgeProperty &&) = default;
    RREdgeProperty &operator=(const RREdgeProperty &) = delete;
    RREdgeProperty &operator=(RREdgeProperty &&) = default;

    [[nodiscard]] int64_t Size() const;

    [[nodiscard]] bool IsCanonical() const;
    [[nodiscard]] bool IsUnique() const { return unique; }
    [[nodiscard]] const MDBGSeq &Seq() const { return seq; }

    [[nodiscard]] RREdgeIndexType Index() const { return index; }

    void Merge(RRVertexProperty vertex, RREdgeProperty rhs);

    MDBGSeq ExtractSeqPrefix(size_t len);
    MDBGSeq ExtractSeqSuffix(size_t len);

    void ShortenWithEmptySeq(size_t len);

    [[nodiscard]] double Cov() const;
};

bool operator==(const RREdgeProperty &lhs, const RREdgeProperty &rhs);
bool operator!=(const RREdgeProperty &lhs, const RREdgeProperty &rhs);

RREdgeProperty Add(const RRVertexProperty &lhs, const RRVertexProperty &rhs,
                   RREdgeIndexType index);

std::ostream &operator<<(std::ostream &os, const RREdgeProperty &edge_property);

// ---------- SuccinctEdgeInfo ----------

struct SuccinctEdgeInfo {
    RRVertexType start_ind{0};
    RRVertexType end_ind{0};
    const dbg::Edge *edge{nullptr};
    bool unique{false};
};

bool operator==(const SuccinctEdgeInfo &lhs, const SuccinctEdgeInfo &rhs);

bool operator!=(const SuccinctEdgeInfo &lhs, const SuccinctEdgeInfo &rhs);
} // End namespace repeat_resolution