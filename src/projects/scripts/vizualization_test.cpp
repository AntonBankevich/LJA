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
    ObjInfo<EID> edgeInfo = EdgePrintStyles<multigraph::MGTraits>::defaultDotInfo();
    Printer<multigraph::MGTraits> printer(vertexInfo, edgeInfo);
    printer.printDot("printer_test.dot", cmp);
    printer.printGFA("printer_test.gfa", cmp);
    return 0;
}