//
// Created by Andrey Bzikadze on 11/25/21.
//

#pragma once

#include "mdbg.hpp"

namespace repeat_resolution {

class MDBGSimpleVertexProcessor {

public:
  void process(MultiplexDBG &graph, const RRVertexType &vertex);
};

class MDBGComplexVertexProcessor {

public:
  void process(MultiplexDBG &graph, const RRVertexType &vertex);
};

} // End namespace repeat_resolution