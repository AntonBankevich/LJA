#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <common/dir_utils.hpp>
#include <assembly_graph/data_structures/splitters.hpp>
#include <assembly_graph/visualization.hpp>
#include "dbg/multi_graph.hpp"
using namespace multigraph;
int main(int argc, char **argv) {
    MultiGraph mg =MultiGraphHelper::LoadGFA(argv[1], true);
    size_t cnt = 1;
    std::experimental::filesystem::path dir = argv[2];
    ensure_dir_existance(dir);
    std::cout << "dbg " << mg.size() << " " << mg.edgeCount() << std::endl;
    std::cout << "component\tsize" << std::endl;
    ag::Printer printer;
    for(const ag::Component &comp : ag::CCSplitter().splitGraph(mg)) {
        std::cout << cnt << ".gfa\t" << comp.size() << std::endl;
        printer.printGFA(dir / (itos(cnt) + ".gfa"), comp, false);
        cnt++;
    }
    return 0;
}
