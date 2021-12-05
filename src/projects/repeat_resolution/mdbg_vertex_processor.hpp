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

  std::pair<std::unordered_map<EdgeIndexType, RRVertexType>,
            std::vector<RRVertexType>>
  split_vertex(MultiplexDBG &graph, const RRVertexType &vertex);

public:
  void process(MultiplexDBG &graph, const RRVertexType &vertex);
};

} // End namespace repeat_resolution