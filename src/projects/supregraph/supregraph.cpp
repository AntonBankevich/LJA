#include "supregraph.hpp"

spg::SPGVertex &spg::SupreGraph::outerEdgeToVertex(spg::SPGEdge &edge) {
    VERIFY(edge.isOuter());
    Vertex &newv = addSPGVertex(edge.getSeq(), false, false, false);
    edge.getStart().addSPEdgeLockFree(newv);
    if(newv != newv.rc())
        newv.addSPEdgeLockFree(edge.getFinish());
    edge.getStart().removeEdgeLockFree(edge);
    return newv;
}

void spg::SupreGraph::IsolateAndMark(spg::SPGVertex &v) {
    while(v.outDeg() != 0) v.removeEdgeLockFree(v.front());
    while(v.rc().outDeg() != 0) v.rc().removeEdgeLockFree(v.rc().front());
    v.mark();
    v.rc().mark();
}
