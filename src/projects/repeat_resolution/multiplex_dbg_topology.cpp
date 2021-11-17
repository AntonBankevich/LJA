//
// Created by Andrey Bzikadze on 11/10/21.
//

#include "multiplex_dbg_topology.hpp"
using namespace repeat_resolution;

std::ostream &repeat_resolution::operator<<(std::ostream &os,
                                            const RRVertexType &vertex) {
  os << '"' << vertex.index << "\\n" << vertex.len << '"';
  return os;
}

std::ostream &
repeat_resolution::operator<<(std::ostream &os,
                              const RREdgeProperty &edge_property) {
  os << edge_property.size() << "\\n" << edge_property.is_unique();
  return os;
}
