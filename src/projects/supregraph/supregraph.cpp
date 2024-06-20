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
        Edge &out = addSPEdgeLockFree(p.incoming().getStart(), newv);
        out.setCorporeal(false);
        out.rc().setCorporeal(false);
        if(newv != newv.rc()) {
            Edge &inc = addSPEdgeLockFree(newv, p.outgoing().getFinish());
            inc.setCorporeal(false);
            inc.rc().setCorporeal(false);
        }
    }
    fireResolveVertex(core, result);
    if(core != core.rc())
        fireResolveVertex(core.rc(), result.RC());
    isolateAndMark(core);
    for(Vertex &v : result.newVertices()) {
        v.front().setCorporeal(true);
        v.front().rc().setCorporeal(true);
        v.rc().front().setCorporeal(true);
        v.rc().front().rc().setCorporeal(true);
    }
    return std::move(result);
}

spg::SPGVertex &spg::SupreGraph::addSPGVertex(Sequence seq, bool cyclic, bool inf_left, bool inf_right, int id) {
    Vertex &res =  addVertex(std::move(seq), SPGVertexData(cyclic, inf_left, inf_right), id);
    return res;
}

spg::SPGVertex &spg::SupreGraph::mergePath(const spg::GraphPath &path) {
    VERIFY(path.startClosed() && path.endClosed());
    Sequence seq = path.Seq();
    if(path.truncLen() == 0 || path.RC().truncLen() == 0) {
        spg::Edge &edge = mergePathToEdge(path);
        return edge.isSuffix() ? edge.getStart() : edge.getFinish();
    }
    SequenceBuilder sb;
    for(Edge &edge : path.edges()) {
        sb.append(edge.getCode());
    }
    Vertex &res = addSPGVertex(seq, false, false, false);
    Edge &inc_edge = addSPEdgeLockFree(path.getStart(), res);
    inc_edge.setCorporeal(false);
    inc_edge.rc().setCorporeal(false);
    Edge &out_edge = addSPEdgeLockFree(res, path.getFinish());
    out_edge.setCorporeal(false);
    out_edge.rc().setCorporeal(false);
    fireMergePath(path.asEdgeIds(), res);
    if(res != res.rc())
        fireMergePath(path.RC().asEdgeIds(), res.rc());
    isolateAndMark(path.innerVertices().begin(), path.innerVertices().end());
    inc_edge.setCorporeal(true);
    inc_edge.rc().setCorporeal(true);
    out_edge.setCorporeal(true);
    out_edge.rc().setCorporeal(true);
    return res;
}

spg::SPGVertex &spg::SupreGraph::mergeLoop(const spg::GraphPath &path) {
    VERIFY(path.getStart() == path.getFinish() || (!path.empty() && path.frontEdge() == path.backEdge() &&
            path.rightCut() +
            path.leftCut() == path.backEdge().truncSize()));
    Sequence seq = path.Seq();
    VertexId resId;
    Vertex &res = addSPGVertex(seq, true, false, false);
    ag::GraphPath<SPGTraits> new_path(path.getStart());
    for(Edge &edge : path.edges())
        new_path += edge;
    fireMergeLoop(new_path, res);
    if(res != res.rc())
        fireMergeLoop(new_path.RC(), res.rc());
    isolateAndMark(path.vertices().begin(), path.vertices().end());
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
