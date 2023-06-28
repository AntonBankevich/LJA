#include "dbg/multi_graph.hpp"
#include <common/cl_parser.hpp>

using namespace multigraph;
int main(int argc, char **argv) {
    MultiGraph mg = MultiGraphHelper::LoadGFA(argv[1], false);
    size_t k = std::stoull(argv[2]);
    mg = MultiGraphHelper::TransformToVertexGraph(mg, k);
    std::cout << "dbg " << mg.size() << " " << mg.edgeNumber() << std::endl;
    for(Vertex & v: mg.vertices()) {
        if(v.inDeg() == 1 && v.outDeg() == 1 && v[0] != v.rc()[0].rc()) {
            std::cout << "1-in-1-out vertex " << v.getId() << " " << v[0].getId() << v.rc()[0].rc().getId() << std::endl;
        }
    }
    return 0;
}
