#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <common/dir_utils.hpp>
#include "lja/multi_graph.hpp"
using namespace multigraph;
int main(int argc, char **argv) {
    MultiGraph mg;
    size_t k = std::stoull(argv[2]);
    mg.LoadGFA(argv[1], false);
    mg = mg.DBG(k);
    std::cout << "dbg " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    for(auto & it: mg.vertices) {
        Vertex *v = &it.second;
        if(v->inDeg() == 1 && v->outDeg() == 1 && v->outgoing[0] != v->rc->outgoing[0]->rc) {
            std::cout << "1-in-1-out vertex " << v->id << " " << v->outgoing[0]->getId() << v->rc->outgoing[0]->rc->getId() << std::endl;
        }
    }
    return 0;
}
