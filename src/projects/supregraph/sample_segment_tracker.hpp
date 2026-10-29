#pragma once

#include "assembly_graph/assembly_graph.hpp"

namespace spg {
    using ag::AssemblyGraph;
    using ag::Vertex;
    using ag::Edge;
    using ag::RAGraphPath;
    using ag::VertexResolutionResult;

//    Owns the identity/structure of each vertex's sample (anchor) segments -- the minimal substrings that
//    uniquely identify a vertex -- and keeps their windows correctly shifted/combined as the graph is
//    edited. Deliberately a pure graph listener: it has no business knowing about reads (anchors are a
//    property of graph structure only), so every segment it creates here starts with zero support --
//    CoverageTracker is the only thing that ever writes a real support value (either by replaying reads,
//    or, for fireResolveVertex's ambiguous case where no read can be mid-reroute yet to replay, by patching
//    straight from VertexResolutionResult's own getSupport() -- see CoverageTracker::fireResolveVertex).
//    Supregraph-only: DBG already has its own implicit k+1-mer anchoring and never attaches this listener.
//    Meant to be constructed once, right before the DBG->SPG conversion loop, so its fireEdgeToSupreVertex
//    override catches every edge as it converts.
    class SampleSegmentTracker : public ag::ResolutionListener {
    private:
        size_t k;
    public:
        SampleSegmentTracker(AssemblyGraph &graph, size_t k) :
                ag::ResolutionListener(graph, "SampleSegmentTracker"), k(k) {}

        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
//        Supregraphs never merge a path into an edge or split an edge (that operation is DBG-only);
        virtual void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {VERIFY(false);}
        virtual void fireSplitEdge(Edge &edge, const RAGraphPath &split) {VERIFY(false);}

    };
}
