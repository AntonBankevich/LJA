#include "edge_code_listener.hpp"
#include "assembly_graph.hpp"
#include <array>

using namespace ag;

bool EdgeCodeListener::checkConsistency(std::ostream &os) {
    const AssemblyGraph &graph = getGraph();
    size_t errors = 0;
    auto report = [&]() -> std::ostream & {
        errors++;
        return os << "Edge code inconsistency: ";
    };
    for (const Vertex &v : graph.vertices()) {
        std::array<bool, 4> seen = {false, false, false, false};
        for (const Edge &e : v) {
            const std::string name = "edge " + e.getInnerId().str();
            if (e.isCodeEmpty()) {
                if (v.outDeg() != 1)
                    report() << name << ": edge with empty code (or virtual suffix) is not the only outgoing edge of its start\n";
                continue;
            }
            if (e.isSuffix()) {
                report() << name << ": suffix edge has non-empty code\n";
                continue;
            }
            if (e.edge_code[0] != e.truncSeq()[0])
                report() << name << ": code does not start with the edge's first nucleotide\n";
            if (seen[e.edge_code[0]])
                report() << name << ": another outgoing edge of vertex " << v.getInnerId() << " has the same first code nucleotide\n";
            seen[e.edge_code[0]] = true;
        }
    }
    if (errors > 0)
        os << "Edge code inconsistency: " << errors << " violation(s) found" << std::endl;
    return errors == 0;
}

void EdgeCodeListener::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    // The merged vertex has a single incoming edge spelling the whole path; its code is the
    // concatenation of the path edges' codes. The single outgoing edge is a suffix edge (empty code,
    // set by fireAddEdge). Start/end keep their external edges and codes unchanged.
    VERIFY(new_vertex.inDeg() == 1);
    SequenceBuilder sb;
    for(Edge &e : path.edges()) sb.append(e.edge_code);
    new_vertex.incFront().edge_code = sb.BuildSequence();
}

void EdgeCodeListener::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    SequenceBuilder sb;
    for(Edge &e : path.edges()) sb.append(e.edge_code);
    new_edge.edge_code = sb.BuildSequence();
    VERIFY(new_edge.edge_code.startsWith(path.frontEdge().edge_code));
}

void EdgeCodeListener::fireMergePrefixPath(const RAGraphPath &path, Edge &new_edge) {
    SequenceBuilder sb;
    for(Edge &e : path.edges()) sb.append(e.edge_code);
    new_edge.edge_code = sb.BuildSequence();
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
//    The pair in -> out becomes u -> new_vertex -> w: u -> new_vertex carries the whole code and new_vertex -> w, the only
//    outgoing edge of new_vertex, none, so paths through the pair keep their codes.
    for(Vertex &new_vertex : resolution.newVertices()) {
        new_vertex.incFront().edge_code = resolution.get(new_vertex).outgoing().edge_code;
        new_vertex.front().edge_code = Sequence();
    }
    for(const auto &it : resolution.newEdges())
        it.first->edge_code = it.second.incoming().edge_code + it.second.outgoing().edge_code;
}

void EdgeCodeListener::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    v.incFront().edge_code = e.edge_code;
}

// Every new edge replaces a two-edge path through core, so its code is the concatenation of that
// path's codes. Fired once per orientation:
// - forward (edge = inc, the unique incoming edge of core): new edge for out-edge e spells inc->e;
// - mirrored (edge = out, the unique outgoing edge of the rc side): new edge for in-edge e spells e->out.
void EdgeCodeListener::fireFlipVertex(Edge &edge, const std::unordered_map<EdgeId, EdgeId> &extensions) {
    if (edge.isSuffix()) {
        for (Edge &e : edge.getFinish())
            extensions.at(e.getId())->edge_code = edge.edge_code + e.edge_code;
    } else {
        for (Edge &e : edge.getStart().incoming())
            extensions.at(e.getId())->edge_code = e.edge_code + edge.edge_code;
    }
}

void EdgeCodeListener::fireAddEdge(Edge &new_edge) {new_edge.edge_code = new_edge.isSuffix() ? Sequence() : new_edge.truncSeq().Subseq(0, 1);}
