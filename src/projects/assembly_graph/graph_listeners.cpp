#include "graph_listeners.hpp"
#include "assembly_graph.hpp"
using namespace ag;
ResolutionListener::ResolutionListener(ResolutionFire &fire, const std::string &name) : AbstractListener(fire, name) {}

AssemblyGraph &ResolutionListener::getGraph() {
    return static_cast<AssemblyGraph &>(getFire<ResolutionFire>());
}
