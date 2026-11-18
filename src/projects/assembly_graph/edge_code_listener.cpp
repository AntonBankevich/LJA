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

// P's codes go where PathReplacement::oldCodes says; every other new edge keeps the codes of the edge it replaces (Q's
// edge out of D has none of its own).
void EdgeCodeListener::fireReplacePath(const PathReplacement &rep) {
    const Sequence codes = rep.oldCodes();
    const bool to_vs = !rep.hasVP() && !rep.startReplaced();
    const bool to_end = to_vs && !rep.hasVS();
    VERIFY(!to_end || rep.endReplaced() || codes.empty());
    if (rep.hasVP())
        rep.vpEdge().edge_code = codes;
    if (rep.hasVS())
        rep.vsEdge().edge_code = to_vs ? codes : Sequence();
    if (rep.startReplaced())
        for (Edge &in : rep.oldPath().getStart().incoming())
            rep.replacement(in).edge_code = rep.codesToReplacedStart() ? in.edge_code + codes : in.edge_code;
    if (rep.endReplaced())
        for (Edge &out : rep.oldPath().getFinish())
            if (!rep.closesLoop(out))
                rep.replacement(out).edge_code = to_end ? codes + out.edge_code : out.edge_code;
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
//    The pair in -> out becomes u -> new_vertex -> w: u -> new_vertex carries the whole code and new_vertex -> w, the only
//    outgoing edge of new_vertex, none, so paths through the pair keep their codes.
//    A copy of core keeps core's edges' codes: u -> copy that of in, copy -> w that of out.
    for(Vertex &new_vertex : resolution.newVertices()) {
        const InOutEdgePair &pair = resolution.get(new_vertex);
        if (new_vertex.size() == core.size()) {
            new_vertex.incFront().edge_code = pair.incoming().edge_code;
            new_vertex.front().edge_code = pair.outgoing().edge_code;
        } else {
            new_vertex.incFront().edge_code = pair.outgoing().edge_code;
            new_vertex.front().edge_code = Sequence();
        }
    }
    for(const auto &it : resolution.newEdges())
        it.first->edge_code = it.second.incoming().edge_code + it.second.outgoing().edge_code;
}

void EdgeCodeListener::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    v.incFront().edge_code = e.edge_code;
}

void EdgeCodeListener::fireAddEdge(Edge &new_edge) {new_edge.edge_code = new_edge.isSuffix() ? Sequence() : new_edge.truncSeq().Subseq(0, 1);}
