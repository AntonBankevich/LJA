//
// Created by Andrey Bzikadze on 11/25/21.
//

#pragma once

#include "mdbg.hpp"

namespace repeat_resolution {

class MDBGSimpleVertexProcessor {

  void process_0in_1pout(MultiplexDBG &graph, const RRVertexType &vertex);
  void process_1pin_0out(MultiplexDBG &graph, const RRVertexType &vertex);

public:
  void process(MultiplexDBG &graph, const RRVertexType &vertex,
               uint64_t n_iter);
};

class MDBGComplexVertexProcessor {

  std::unordered_map<EdgeIndexType, RRVertexType>
  split_vertex(MultiplexDBG &graph, const RRVertexType &vertex);

  void process_11(
      MultiplexDBG &graph, const RRVertexType &vertex,
      const RRVertexType &left_vertex, const RRVertexType &right_vertex,
      MultiplexDBG::NeighborsIterator e1_it,
      MultiplexDBG::NeighborsIterator e2_it, const EdgeIndexType &edge1,
      const EdgeIndexType &edge2,
      std::unordered_map<EdgeIndexType, EdgeIndexType> &where_edge_merged);

  void process_not11(
      MultiplexDBG &graph, const RRVertexType &vertex,
      const RRVertexType &left_vertex, const RRVertexType &right_vertex,
      MultiplexDBG::NeighborsIterator e1_it,
      MultiplexDBG::NeighborsIterator e2_it, const EdgeIndexType &edge1,
      const EdgeIndexType &edge2,
      const std::unordered_set<EdgeIndexType> &edge1_neighbors,
      const std::unordered_set<EdgeIndexType> &edge2_neighbors,
      std::unordered_map<EdgeIndexType, EdgeIndexType> &where_edge_merged);

public:
  void process(MultiplexDBG &graph, const RRVertexType &vertex);
};

} // End namespace repeat_resolution