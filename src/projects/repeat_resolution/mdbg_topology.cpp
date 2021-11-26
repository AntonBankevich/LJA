//
// Created by Andrey Bzikadze on 11/10/21.
//

#include "mdbg_topology.hpp"
using namespace repeat_resolution;

std::ostream &repeat_resolution::operator<<(std::ostream &os,
                                            const RRVertexProperty &vertex) {
  os << vertex.len;
  return os;
}

std::ostream &
repeat_resolution::operator<<(std::ostream &os,
                              const RREdgeProperty &edge_property) {
  os << edge_property.size() << "\\n" << edge_property.is_unique();
  return os;
}

bool repeat_resolution::operator==(const RRVertexProperty &lhs,
                                   const RRVertexProperty &rhs) {
  return lhs.len == rhs.len and lhs.frozen == rhs.frozen;
}

bool repeat_resolution::operator==(const RREdgeProperty &lhs,
                                   const RREdgeProperty &rhs) {
  return lhs.get_index() == rhs.get_index();
}

bool repeat_resolution::operator!=(const RREdgeProperty &lhs,
                                   const RREdgeProperty &rhs) {
  return not (lhs == rhs);
}

RREdgeProperty repeat_resolution::add(const RREdgeProperty &lhs,
                                      const RREdgeProperty &rhs,
                                      const uint64_t vertex_size,
                                      const uint64_t index) {
  VERIFY(lhs.size() > vertex_size);
  VERIFY(rhs.size() > vertex_size);
  const bool unique = lhs.is_unique() or rhs.is_unique();
  std::list<char> new_seq = [&lhs, &rhs, &vertex_size]() {
    const std::list<char> &seq = lhs.get_seq();
    auto lhs_it = lhs.get_seq().rbegin();
    for (size_t i = 0; i < vertex_size - 1; ++i) {
      lhs_it++;
    }
    char lhs_char = *lhs_it;
    std::list<char> new_seq{lhs_char};
    auto rhs_it = rhs.get_seq().begin();
    for (uint64_t i = 0; i < vertex_size; ++i) {
      new_seq.emplace_back(*rhs_it);
      ++rhs_it;
    }
    return new_seq;
  }();
  return {index, std::move(new_seq), unique};
}
