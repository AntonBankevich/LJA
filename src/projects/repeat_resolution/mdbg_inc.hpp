//
// Created by Andrey Bzikadze on 11/25/21.
//

#pragma once

#include "mdbg.hpp"
#include "mdbg_vertex_processor.hpp"

namespace repeat_resolution {

class MultiplexDBGIncreaser {
  uint64_t start_k{1};
  uint64_t saturating_k{1};
  logging::Logger &logger;
  bool debug{true};
  MDBGSimpleVertexProcessor simple_vertex_processor;
  MDBGComplexVertexProcessor complex_vertex_processor;

private:
  void process_vertex(MultiplexDBG &graph, const RRVertexType &vertex);
  static void collapse_short_edges_into_vertices(MultiplexDBG & graph);
  static void collapse_edge(MultiplexDBG &graph, MultiplexDBG::ConstIterator s_it,
                     MultiplexDBG::NeighborsIterator e_it);

public:
  MultiplexDBGIncreaser(uint64_t start_k, uint64_t saturating_k,
                        logging::Logger &logger, bool debug);

  void Increment(MultiplexDBG &graph);

  void IncrementN(MultiplexDBG &graph, uint64_t N);
};

} // End namespace repeat_resolution