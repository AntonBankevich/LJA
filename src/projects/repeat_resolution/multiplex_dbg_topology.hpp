//
// Created by Andrey Bzikadze on 10/14/21.
//

#pragma once

#include "common/verify.hpp"
#include <functional>
#include <graphlite/graphlite.hpp>
#include <graphlite/serialize.hpp>

namespace repeat_resolution {
struct RRVertexType {
  uint64_t index{0};
  uint64_t len{0};
  bool frozen{false};
};

inline bool operator==(const RRVertexType &lhs, const RRVertexType &rhs) {
  return lhs.index == rhs.index;
}
inline bool operator!=(const RRVertexType &lhs, const RRVertexType &rhs) {
  return !operator==(lhs, rhs);
}
inline bool operator<(const RRVertexType &lhs, const RRVertexType &rhs) {
  return lhs.index < rhs.index;
}
inline bool operator>(const RRVertexType &lhs, const RRVertexType &rhs) {
  return operator<(rhs, lhs);
}
inline bool operator<=(const RRVertexType &lhs, const RRVertexType &rhs) {
  return !operator>(lhs, rhs);
}
inline bool operator>=(const RRVertexType &lhs, const RRVertexType &rhs) {
  return !operator<(lhs, rhs);
}

std::ostream &operator<<(std::ostream &os, const RRVertexType &vertex);
} // namespace repeat_resolution

namespace std {
template <> struct hash<repeat_resolution::RRVertexType> {
  size_t operator()(const repeat_resolution::RRVertexType &s) const noexcept {
    return hash<uint64_t>()(s.index);
  }
};
} // End namespace std

namespace repeat_resolution {

using EdgeIndexType = uint64_t;

class RREdgeProperty {
  uint64_t index{0};
  std::list<char> seq;
  bool unique{false};

public:
  RREdgeProperty(const uint64_t index, std::list<char> seq, bool unique)
      : index{index}, seq{std::move(seq)}, unique{unique} {}

  [[nodiscard]] uint64_t size() const { return seq.size(); }

  [[nodiscard]] bool is_unique() const { return unique; }

  [[nodiscard]] uint64_t get_index() const { return index; }
  [[nodiscard]] std::list<char> get_seq() const { return seq; }

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

  void append(char c) { seq.push_back(c); }
  void prepend(char c) { seq.push_front(c); }

  void merge(RREdgeProperty rhs, uint64_t overlap_len) {
    VERIFY(size() > overlap_len and rhs.size() > overlap_len);
    auto lhs_it = [this, &overlap_len]() {
      auto it = seq.rbegin();
      for (size_t i = 0; i < overlap_len; ++i) {
        it++;
      }
      return it.base();
    }();
    auto rhs_it = rhs.seq.begin();
    for (/*empty*/; lhs_it != seq.end(); ++lhs_it, ++rhs_it) {
      VERIFY(*rhs_it == *lhs_it);
    }
    for (size_t i = 0; i < overlap_len; ++i) {
      rhs.seq.pop_front();
    }
    seq.splice(seq.end(), rhs.seq);
    if (rhs.unique) {
      unique = true;
    }
  }
};

std::ostream &operator<<(std::ostream &os, const RREdgeProperty &edge_property);

struct SuccinctEdgeInfo {
  RRVertexType start;
  RRVertexType end;
  std::list<char> seq;
  bool unique{false};
};

inline bool operator==(const SuccinctEdgeInfo &lhs,
                       const SuccinctEdgeInfo &rhs) {
  return lhs.start == rhs.start and lhs.end == rhs.end and
         lhs.seq == rhs.seq and lhs.unique == rhs.unique;
}

inline bool operator!=(const SuccinctEdgeInfo &lhs,
                       const SuccinctEdgeInfo &rhs) {
  return !operator==(lhs, rhs);
}
} // End namespace repeat_resolution