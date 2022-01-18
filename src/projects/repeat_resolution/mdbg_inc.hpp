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
    void ProcessVertex(MultiplexDBG &graph, const RRVertexType &vertex,
                       uint64_t max_iter,
                       std::set<Sequence> &merged_self_loops);
    static void CollapseShortEdgesIntoVertices(MultiplexDBG &graph);
    static void CollapseEdge(MultiplexDBG &graph,
                             MultiplexDBG::ConstIterator s_it,
                             MultiplexDBG::NeighborsIterator e_it);
    [[nodiscard]] uint64_t GetNiterWoComplex(const MultiplexDBG &graph) const;

 public:
    MultiplexDBGIncreaser(uint64_t start_k, uint64_t saturating_k,
                          logging::Logger &logger, bool debug);

    void Increase(MultiplexDBG &graph,
                  bool unite_simple,
                  uint64_t max_iter = 1);

    void IncreaseN(MultiplexDBG &graph, uint64_t N, bool unite_simple);

    void IncreaseUntilSaturation(MultiplexDBG &graph, bool unite_simple = true);
};

} // End namespace repeat_resolution