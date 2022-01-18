//
// Created by Andrey Bzikadze on 11/25/21.
//

#pragma once

#include "mdbg.hpp"

namespace repeat_resolution {

class MDBGSimpleVertexProcessor {

    void Process0In1Pout(MultiplexDBG &graph, const RRVertexType &vertex);
    void Process1Pin0Out(MultiplexDBG &graph, const RRVertexType &vertex);

 public:
    void Process(MultiplexDBG &graph, const RRVertexType &vertex,
                 uint64_t n_iter);
};

class MDBGComplexVertexProcessor {

    std::pair<std::unordered_map<RREdgeIndexType, RRVertexType>,
              std::vector<RRVertexType>>
    SplitVertex(MultiplexDBG &graph, const RRVertexType &vertex);

 public:
    void Process(MultiplexDBG &graph,
                 const RRVertexType &vertex,
                 std::set<Sequence> &merged_self_loops);
};

} // End namespace repeat_resolution