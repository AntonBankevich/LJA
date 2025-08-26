//
// Created by Anton Zamyatin on 7/8/24.
//

#include "assembly_graph/visualization.hpp"
#include "dbg/multi_graph.hpp"

int main(int argc, char **argv) {
    multigraph::MultiGraph mdbg = multigraph::MultiGraphHelper::LoadGFA(argv[1], true);
    logging::Logger logger;
    mdbg = multigraph::MultiGraphHelper::TransformToEdgeGraph(logger, mdbg, 5001);
    multigraph::MultiGraphHelper::printDot2(mdbg, "print2dot.dot");
    ag::Component cmp(mdbg);
    ag::ObjInfo<multigraph::Vertex> vertexInfo = ag::VertexPrintStyles::defaultDotInfo();
    ag::ObjInfo<multigraph::Edge> edgeInfo = ag::EdgePrintStyles::defaultDotInfo();
    ag::Printer printer(vertexInfo, edgeInfo);
    printer.printDot("printer_test.dot", mdbg);
    printer.printGFA("printer_test.gfa", mdbg);
    printer.printExtendedGFA("printer_test_ext.gfa", mdbg);
    return 0;
}