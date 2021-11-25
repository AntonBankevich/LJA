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

struct RRVertexProperty {
  uint64_t len{0};
  bool frozen{false};

  void freeze() {
    frozen = true;
  }
};

std::ostream &operator<<(std::ostream &os, const RRVertexProperty &vertex);

bool operator==(const RRVertexProperty &lhs, const RRVertexProperty &rhs);

using EdgeIndexType = uint64_t;

class RREdgeProperty {
  EdgeIndexType index{0};
  std::list<char> seq;
  bool unique{false};

public:
  RREdgeProperty(const EdgeIndexType index, std::list<char> seq, bool unique)
      : index{index}, seq{std::move(seq)}, unique{unique} {}

  [[nodiscard]] uint64_t size() const { return seq.size(); }

  [[nodiscard]] bool is_unique() const { return unique; }

  [[nodiscard]] EdgeIndexType get_index() const { return index; }
  [[nodiscard]] const std::list<char> &get_seq() const { return seq; }

  void assert_incidence(const RREdgeProperty &rhs,
                        const uint64_t overlap_len) const {
    VERIFY(size() > overlap_len and rhs.size() > overlap_len);
    auto it = [this, overlap_len]() {
      auto rit = seq.rbegin();
      for (uint64_t i = 0; i < overlap_len; ++i) {
        ++rit;
      }
      return rit.base();
    }();
    auto it_rhs = rhs.seq.begin();
    while (it != seq.end()) {
      VERIFY(*it == *it_rhs);
      ++it, ++it_rhs;
    }
  }

  void append(const RREdgeProperty &rhs, uint64_t overlap_len) {
    VERIFY(rhs.size() > overlap_len);
    auto it = rhs.get_seq().begin();
    std::advance(it, overlap_len);
    seq.push_back(*it);
  }

  void prepend(const RREdgeProperty &lhs, uint64_t overlap_len) {
    VERIFY(lhs.size() > overlap_len);
    auto it = lhs.get_seq().rbegin();
    std::advance(it, overlap_len);
    seq.push_front(*it);
  }

  void merge(RREdgeProperty rhs, uint64_t overlap_len) {
    VERIFY(size() > overlap_len and rhs.size() > overlap_len);
    auto lhs_it = [this, &overlap_len]() {
      auto it = seq.rbegin();
      for (size_t i = 0; i < overlap_len; ++i) {
        it++;
      }
      return it.base();
    }();

    for (auto rhs_it = rhs.seq.begin(); lhs_it != seq.end();
         ++lhs_it, ++rhs_it) {
      VERIFY(*rhs_it == *lhs_it);
    }
    for (size_t i = 0; i < overlap_len; ++i) {
      rhs.seq.pop_front();
    }
    seq.splice(seq.end(), std::move(rhs.seq));
    if (rhs.unique) {
      unique = true;
    }
  }
};

bool operator==(const RREdgeProperty &lhs, const RREdgeProperty &rhs);
bool operator!=(const RREdgeProperty &lhs, const RREdgeProperty &rhs);

RREdgeProperty add(const RREdgeProperty &lhs, const RREdgeProperty &rhs,
                   uint64_t overlap_len, EdgeIndexType index);

std::ostream &operator<<(std::ostream &os, const RREdgeProperty &edge_property);

struct SuccinctEdgeInfo {
  RRVertexType start_ind{0};
  RRVertexProperty start_prop;
  RRVertexType end_ind{0};
  RRVertexProperty end_prop;
  std::list<char> seq;
  bool unique{false};
};

inline bool operator==(const SuccinctEdgeInfo &lhs,
                       const SuccinctEdgeInfo &rhs) {
  return lhs.start_ind == rhs.start_ind and lhs.start_prop == rhs.start_prop and
         lhs.end_ind == rhs.end_ind and lhs.end_prop == rhs.end_prop and
         lhs.seq == rhs.seq and lhs.unique == rhs.unique;
}

inline bool operator!=(const SuccinctEdgeInfo &lhs,
                       const SuccinctEdgeInfo &rhs) {
  return !operator==(lhs, rhs);
}
} // End namespace repeat_resolution