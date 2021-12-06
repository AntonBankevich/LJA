//
// Created by Andrey Bzikadze on 10/14/21.
//

#pragma once

#include "common/verify.hpp"
#include <functional>
#include <graphlite/graphlite.hpp>
#include <graphlite/serialize.hpp>

namespace repeat_resolution {
using RRVertexType = uint64_t;

class RRVertexProperty {
  std::list<char> seq;
  bool frozen{false};

public:
  friend class RREdgeProperty;
  RRVertexProperty(std::list<char> seq, const bool frozen)
      : seq{std::move(seq)}, frozen{frozen} {}

  RRVertexProperty(const RRVertexProperty &) = delete;
  RRVertexProperty(RRVertexProperty &&) = default;
  RRVertexProperty &operator=(const RRVertexProperty &) = delete;
  RRVertexProperty &operator=(RRVertexProperty &&) = default;

  [[nodiscard]] uint64_t size() const { return seq.size(); }
  [[nodiscard]] bool IsFrozen() const { return frozen; }
  [[nodiscard]] const std::list<char> &Seq() const { return seq; }

  void freeze() { frozen = true; }

  void IncLeft(std::list<char> prefix);
  void IncRight(std::list<char> suffix);
  void DecLeft(uint64_t inc = 1);
  void DecRight(uint64_t inc = 1);

  [[nodiscard]] std::list<char> GetSeqPrefix(size_t len,
                                             int64_t shift = 0) const;
  [[nodiscard]] std::list<char> GetSeqSuffix(size_t len,
                                             int64_t shift = 0) const;
};

std::ostream &operator<<(std::ostream &os, const RRVertexProperty &vertex);

bool operator==(const RRVertexProperty &lhs, const RRVertexProperty &rhs);

using EdgeIndexType = uint64_t;

class RREdgeProperty {
  EdgeIndexType index{0};
  std::list<char> seq;
  int64_t size_{0};
  bool unique{false};

public:
  RREdgeProperty(const EdgeIndexType index, std::list<char> seq,
                 const int64_t size_, bool unique)
      : index{index}, seq{std::move(seq)}, size_{size_}, unique{unique} {}

  RREdgeProperty(const RREdgeProperty &) = delete;
  RREdgeProperty(RREdgeProperty &&) = default;
  RREdgeProperty &operator=(const RREdgeProperty &) = delete;
  RREdgeProperty &operator=(RREdgeProperty &&) = default;

  [[nodiscard]] int64_t size() const;

  [[nodiscard]] bool IsUnique() const { return unique; }

  [[nodiscard]] EdgeIndexType Index() const { return index; }
  [[nodiscard]] const std::list<char> &Seq() const { return seq; }

  void Merge(RRVertexProperty vertex, RREdgeProperty rhs);

  std::list<char> ExtractSeqPrefix(size_t len);
  std::list<char> ExtractSeqSuffix(size_t len);

  void ShortenWithEmptySeq(size_t len);
};

bool operator==(const RREdgeProperty &lhs, const RREdgeProperty &rhs);
bool operator!=(const RREdgeProperty &lhs, const RREdgeProperty &rhs);

RREdgeProperty Add(const RRVertexProperty &lhs, const RRVertexProperty &rhs,
                   EdgeIndexType index);

std::ostream &operator<<(std::ostream &os, const RREdgeProperty &edge_property);

struct SuccinctEdgeInfo {
  RRVertexType start_ind{0};
  RRVertexProperty start_prop;
  RRVertexType end_ind{0};
  RRVertexProperty end_prop;
  int64_t infix_size{0};
  std::list<char> seq;
  bool unique{false};
};

bool operator==(const SuccinctEdgeInfo &lhs, const SuccinctEdgeInfo &rhs);

bool operator!=(const SuccinctEdgeInfo &lhs, const SuccinctEdgeInfo &rhs);
} // End namespace repeat_resolution