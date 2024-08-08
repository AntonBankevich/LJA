#include "supregraph.hpp"
#include "assembly_graph/assembly_graph.hpp"

spg::SPGVertex &spg::SupreGraph::outerEdgeToVertex(spg::SPGEdge &edge) {
    VERIFY(edge.isOuter());
    Vertex &newv = addSPGVertex(edge.getSeq(), false, false, false);
    edge.getStart().addSPEdgeLockFree(newv);
    if(newv != newv.rc())
        newv.addSPEdgeLockFree(edge.getFinish());
    edge.getStart().removeEdgeLockFree(edge);
    return newv;
}

void spg::SupreGraph::IsolateAndMark(spg::SPGVertex &vertex) {
    VERIFY(!vertex.marked());
    for(Vertex &v : ag::ThisAndRC(vertex))
        while(v.outDeg() != 0) {
            fireDeleteEdge(v.front());
            fireDeleteEdge(v.front().rc());
            v.removeEdgeLockFree(v.front());
        }
    for(Vertex &v : ag::ThisAndRC(vertex)) {
        fireDeleteVertex(v);
        v.mark();
    }
}

spg::VertexResolutionResult
spg::SupreGraph::resolveVertex(spg::SPGVertex &core, const spg::VertexResolutionPlan &resolution) {
    VERIFY(core.isCore() && core.inDeg() > 0 && core.outDeg() > 0);
    VertexResolutionResult result(core);
    for(const InOutEdgePair &p : resolution.connectionsUnique()) {
        VERIFY(p.incoming().getFinish() == core);
        if(p.middle() != core)
            continue;
        VERIFY(p.outgoing().getStart() == core);
        VERIFY(p.incoming().isSuffix());
        VERIFY(p.outgoing().isPrefix());
        Sequence seq = p.getSeq();
        Vertex &newv = addSPGVertex(seq, false, false, false);
        result.add(newv, p);
//        for(Vertex &v : ag::ThisAndRC(newv))
//            fireAddVertex(v);
        Edge &out = p.incoming().getStart().addSPEdgeLockFree(newv);
        for(Edge &e : ag::ThisAndRC(out))
            fireAddEdge(e);
        if(newv != newv.rc()) {
            Edge &inc = newv.addSPEdgeLockFree(p.outgoing().getFinish());
            for(Edge &e : ag::ThisAndRC(inc))
                fireAddEdge(e);
        }
    }
    fireResolveVertex(core, result);
    if(core != core.rc())
        fireResolveVertex(core.rc(), result.RC());
    IsolateAndMark(core);
    return std::move(result);
}

spg::SPGVertex &spg::SupreGraph::addSPGVertex(Sequence seq, bool cyclic, bool inf_left, bool inf_right, int id) {
    Vertex &res =  addVertex(std::move(seq), SPGVertexData(cyclic, inf_left, inf_right), id);
    for(Vertex &v : ag::ThisAndRC(res))
        fireAddVertex(v);
    return res;
}

spg::SPGVertex &spg::SupreGraph::mergePath(const spg::GraphPath &path) {
    VERIFY(path.startClosed() && path.endClosed());
    Sequence seq = path.Seq();
    VertexId resId;
    if(path.start().getSeq() == seq)
        resId = path.start().getId();
    if(path.finish().getSeq() == seq)
        resId = path.finish().getId();
    Vertex &res = resId.valid() ? *resId : addSPGVertex(seq, false, false, false);
    if(res != path.start()) {
        Edge &new_edge = path.start().addSPEdgeLockFree(res);
        for(Edge &e : ag::ThisAndRC(new_edge))
            fireAddEdge(e);
    }
    if(res != path.finish() && res != res.rc()) {
        Edge &new_edge = res.addSPEdgeLockFree(path.finish());
        for(Edge &e : ag::ThisAndRC(new_edge))
            fireAddEdge(e);
    }
    fireMergePath(path, res);
    if(res != res.rc())
        fireMergePath(path.RC(), res.rc());
    if(res != res.rc())
        for(size_t i = path.size() - 1; i > 0; i--) {
            IsolateAndMark(path.getVertex(i));
        }
    else
        for(size_t i = path.size() - 1; 2 * i > path.size(); i--) {
            IsolateAndMark(path.getVertex(i));
        }
    return res;
}

spg::SPGVertex &spg::SupreGraph::mergeLoop(const spg::GraphPath &path) {
    VERIFY(path.start() == path.finish() || (path.size() > 0 && path.frontEdge() == path.backEdge() && path.cutRight() + path.cutLeft() == path.backEdge().truncSize()));
    Sequence seq = path.Seq();
    VertexId resId;
    Vertex &res = addSPGVertex(seq, true, false, false);
//    for(Vertex &v : ag::ThisAndRC(res))
//        fireAddVertex(v);
    fireMergeLoop(path, res);
    if(res != res.rc())
        fireMergeLoop(path.RC(), res.rc());
    if(res != res.rc())
        for(size_t i = path.size() - 1; i > 0; i--) {
            IsolateAndMark(path.getVertex(i));
        }
    else
        for(size_t i = path.size() - 1; 2 * i > path.size(); i--) {
            IsolateAndMark(path.getVertex(i));
        }
    return res;
}
