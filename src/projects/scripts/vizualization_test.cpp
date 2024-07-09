//
// Created by Anton Zamyatin on 7/8/24.
//

#include "assembly_graph/visualization.hpp"
#include "dbg/multi_graph.hpp"

int main(int argc, char **argv) {
    multigraph::MultiGraph mdbg = multigraph::MultiGraphHelper::LoadGFA(argv[1], true);
    mdbg = multigraph::MultiGraphHelper::TransformToEdgeGraph(mdbg, 5001);
    multigraph::MultiGraphHelper::printDot2(mdbg, "print2dot.dot");
    typedef typename ag::BaseVertex<multigraph::MGTraits>::VertexId VID;
    typedef typename ag::BaseEdge<multigraph::MGTraits>::EdgeId EID;
    ag::Component<multigraph::MGTraits> cmp(mdbg);
    ObjInfo<VID> vertexInfo = VertexPrintStyles<multigraph::MGTraits>::defaultDotInfo(cmp);
    ObjInfo<EID> edgeInfo = ObjInfo<EID>::Tooltiper([](const EID &eid){return std::string("black");});
    Printer<multigraph::MGTraits> printer(vertexInfo, edgeInfo);
    std::ofstream os("printer_test.dot");
    printer.printDot(os, cmp);
    os.close();
    return 0;
}