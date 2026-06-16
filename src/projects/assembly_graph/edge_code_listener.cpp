#include "edge_code_listener.hpp"

using namespace ag;

void EdgeCodeListener::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    VERIFY(new_vertex.inDeg() == 1);
    SequenceBuilder sb;
    for(Edge &e : path.edges()) sb.append(e.edge_code);
    new_vertex.rc().front().rc().edge_code = sb.BuildSequence();
}

void EdgeCodeListener::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    SequenceBuilder sb;
    for(Edge &e : path.edges()) sb.append(e.edge_code);
    new_edge.edge_code = sb.BuildSequence();
    VERIFY(new_edge.edge_code.startsWith(path.frontEdge().edge_code));
}

void EdgeCodeListener::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &,
                                           const AlignmentForm &) {
    new_edge.edge_code = left.edge_code + right.edge_code;
}

void EdgeCodeListener::fireMergeLoop(const GraphPath &path, Vertex &new_vertex) {}

void EdgeCodeListener::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
    for(Edge &e : split.edges()) {
        if(e != split.frontEdge())
            e.edge_code = e.getCode() + edge.edge_code;
        else
            e.edge_code = edge.edge_code;
    }
}

void EdgeCodeListener::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    for(Vertex &new_vertex : resolution.newVertices()) {
        new_vertex.rc().front().rc().edge_code = resolution.get(new_vertex).outgoing().edge_code;
    }
}

void EdgeCodeListener::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    v.rc().front().rc().edge_code = e.edge_code;
}

void EdgeCodeListener::fireAddEdge(Edge &new_edge) {new_edge.edge_code = new_edge.isSuffix() ? Sequence() : new_edge.truncSeq().Subseq(0, 1);}
