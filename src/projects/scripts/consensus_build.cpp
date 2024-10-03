//
// Created by Anton Zamyatin on 6/18/24.
//
#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <common/dir_utils.hpp>
#include <error_correction/diploidy_analysis.hpp>

#include "assembly_graph/visualization.hpp"
#include "dbg/multi_graph.hpp"
#include "dbg/visualization.hpp"

using namespace multigraph;
typedef typename ag::BaseVertex<multigraph::MGTraits>::VertexId VID;
typedef typename ag::BaseEdge<multigraph::MGTraits>::EdgeId EID;

int main(int argc, char **argv) {
    MultiGraph mdbg = MultiGraphHelper::LoadGFA(argv[1], true);
    mdbg = MultiGraphHelper::TransformToEdgeGraph(mdbg, 5001);
    /*MultiGraphHelper::printDot2(mdbg, "test.dot");
    std::ofstream os("graph.txt");
    for (auto &v: mdbg.vertices()) {
        os << v.getInnerId();
        for (auto &e: v) {
            os << " " << e.getFinish().getInnerId();
        }
        os << std::endl;
    }
    os.close();*/
    BulgePathFinder<MGTraits> finder(mdbg, -1);

    std::unordered_set<ConstVertexId> vSet;
    for (BulgePath<MGTraits> &path : finder.paths) {
        std::cout << "bulge path len: " << path.size()  << std::endl;
        for (size_t i = 0; i < path.size() ; ++i) {
            std::cout << path[i].first->getStart().getInnerId() << ':'
                      << path[i].first->getFinish().getInnerId() << std::endl;
            vSet.emplace(path[i].first->getStart().getId());
            vSet.emplace(path[i].first->getFinish().getId());
        }
    }
    ag::Component<MGTraits> cmp(mdbg);
    ObjInfo<Vertex> vertexInfo = VertexPrintStyles<MGTraits>::defaultLabeler() +
                              VertexPrintStyles<MGTraits>::vertexSetColorer(vSet) +
                              VertexPrintStyles<MGTraits>::defaultTooltiper();
    ObjInfo<Edge> edgeInfo = ObjInfo<Edge>::Tooltiper([](const Edge &e){return std::string("black");});
    Printer<multigraph::MGTraits> printer(vertexInfo, edgeInfo);
    std::ofstream os("path.dot");
    printer.printDot(os, cmp);
    os.close();
    return 0;
}
