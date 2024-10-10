#include "supregraph.hpp"
#include "assembly_graph/assembly_graph.hpp"

spg::SPGVertex &spg::SupreGraph::outerEdgeToVertex(spg::SPGEdge &edge) {
    VERIFY(edge.isOuter());
    Vertex &newv = addSPGVertex(edge.getSeq(), false, false, false);
    addSPEdgeLockFree(edge.getStart(), newv);
    if(newv != newv.rc())
        addSPEdgeLockFree(newv, edge.getFinish());
    removeEdgeLockFree(edge);
    return newv;
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
        Edge &out = addSPEdgeLockFree(p.incoming().getStart(), newv);
        for(Edge &e : ag::ThisAndRC(out))
            fireAddEdge(e);
        if(newv != newv.rc()) {
            Edge &inc = addSPEdgeLockFree(newv, p.outgoing().getFinish());
            for(Edge &e : ag::ThisAndRC(inc))
                fireAddEdge(e);
        }
    }
    fireResolveVertex(core, result);
    if(core != core.rc())
        fireResolveVertex(core.rc(), result.RC());
    isolateAndMark(core);
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
        Edge &new_edge = addSPEdgeLockFree(path.start(), res);
    }
    if(res != path.finish() && res != res.rc()) {
        Edge &new_edge = addSPEdgeLockFree(res, path.finish());
    }
    fireMergePath(path, res);
    if(res != res.rc())
        fireMergePath(path.RC(), res.rc());
    if(res != res.rc())
        for(size_t i = path.size() - 1; i > 0; i--) {
            isolateAndMark(path.getVertex(i));
        }
    else
        for(size_t i = path.size() - 1; 2 * i > path.size(); i--) {
            isolateAndMark(path.getVertex(i));
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
            isolateAndMark(path.getVertex(i));
        }
    else
        for(size_t i = path.size() - 1; 2 * i > path.size(); i--) {
            isolateAndMark(path.getVertex(i));
        }
    return res;
}

spg::SPGEdge &spg::SupreGraph::addSPEdgeLockFree(spg::SPGVertex &start, spg::SPGVertex &end,
                                                 ag::BaseEdge<spg::SPGTraits>::id_type eid,
                                                 ag::BaseEdge<spg::SPGTraits>::id_type rcid) {
    Sequence tseq = end.getSeq().Subseq(std::min(start.size(), end.size()));
    Sequence rctseq = start.rc().getSeq().Subseq(std::min(start.size(), end.size()));
    return addEdgeLockFree(start, end, tseq, rctseq, EdgeData(), eid, rcid);
}

spg::SPGEdge & spg::SupreGraph::addSPEdge(spg::SPGVertex &start, spg::SPGVertex &end, ag::BaseEdge<spg::SPGTraits>::id_type eid,
                           ag::BaseEdge<spg::SPGTraits>::id_type rcid) {
    Sequence tseq = end.getSeq().Subseq(std::min(start.size(), end.size()));
    Sequence rctseq = start.rc().getSeq().Subseq(std::min(start.size(), end.size()));
    return addEdge(start, end, tseq, rctseq, EdgeData(), eid, rcid);
}
