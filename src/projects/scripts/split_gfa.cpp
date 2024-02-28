#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <common/dir_utils.hpp>
#include "dbg/multi_graph.hpp"
using namespace multigraph;
int main(int argc, char **argv) {
    MultiGraph mg =MultiGraphHelper::LoadGFA(argv[1], true);
    size_t cnt = 1;
    std::experimental::filesystem::path dir = argv[2];
    ensure_dir_existance(dir);
    std::cout << "dbg " << mg.size() << " " << mg.edgeCount() << std::endl;
    std::cout << "component\tsize" << std::endl;
    for(const std::vector<ConstVertexId> &comp : MultiGraphHelper::split(mg)) {
        std::cout << cnt << ".gfa\t" << comp.size() << std::endl;
        MultiGraphHelper::printVertexGFA(dir / (itos(cnt) + ".gfa"), comp);
        cnt++;
    }
    return 0;
}
