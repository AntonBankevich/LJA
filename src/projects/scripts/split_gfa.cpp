#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <common/dir_utils.hpp>
#include "lja/multi_graph.hpp"
using namespace multigraph;
int main(int argc, char **argv) {
    MultiGraph mg;
    mg.LoadGFA(argv[1], true);
    size_t cnt = 1;
    std::experimental::filesystem::path dir = argv[2];
    ensure_dir_existance(dir);
    std::cout << "dbg " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    for(const std::vector<const Vertex *> &comp : mg.split()) {
        std::cout << comp.size() << std::endl;
        mg.printVertexGFA(dir / (itos(cnt) + ".gfa"), comp);
        cnt++;
    }
    mg = mg.DBG(501);
    std::cout << "bulge " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.BulgeSubgraph();
    std::cout << "merge " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.Merge();
    std::cout << "bulge " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.BulgeSubgraph();
    std::cout << "merge " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.Merge();
    std::cout << "final " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg.printEdgeGFA(dir / "final.gfa");
    mg.printEdges(dir / "consensus.fasta", true);

    return 0;
}
