#include "dbg/multi_graph.hpp"
#include <common/cl_parser.hpp>

using namespace multigraph;
int main(int argc, char **argv) {
    MultiGraph mg = MultiGraphHelper::LoadGFA(argv[1], false);
    size_t k = std::stoull(argv[2]);
    mg = MultiGraphHelper::TransformToEdgeGraph(mg, k);
    std::cout << "dbg " << mg.size() << " " << mg.edgeCount() << std::endl;
    for(MGVertex & v: mg.vertices()) {
        if(v.inDeg() == 1 && v.outDeg() == 1 && v.front() != v.rc().front().rc()) {
            std::cout << "1-in-1-out getVertex " << v.getId() << " " << v.front().getId() << v.rc().front().rc().getId() << std::endl;
        }
    }
    return 0;
}
