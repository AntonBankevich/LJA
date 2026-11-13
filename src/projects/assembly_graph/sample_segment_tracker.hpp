#pragma once

#include "graph_listeners.hpp"

namespace ag {
//    Owns the identity/structure of each vertex's sample (anchor) segments -- the minimal substrings that
//    uniquely identify a vertex -- and keeps their windows correctly shifted/combined as the graph is
//    edited. Deliberately a pure graph listener: it has no business knowing about reads (anchors are a
//    property of graph structure only), so every segment it creates here starts with zero support --
//    CoverageTracker is the only thing that ever writes a real support value (either by replaying reads,
//    or, for fireResolveVertex's ambiguous case where no read can be mid-reroute yet to replay, by patching
//    straight from VertexResolutionResult's own getSupport() -- see CoverageTracker::fireResolveVertex).
//    A member of every AssemblyGraph, declared right after EdgeCodeListener, so it is the second listener to
//    react to every operation: all other listeners can rely on sample segments being up to date. Kpomer
//    chunks are created when edgeToSupreVertex converts an edge between two k-mers (k = the graph's k); in a
//    DBG no vertex has sample segments.
    class SampleSegmentTracker : public ResolutionListener {
    public:
        explicit SampleSegmentTracker(ResolutionFire &graph) : ResolutionListener(graph, "SampleSegmentTracker") {}

        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
//        Used on DBGs, whose vertices have no sample segments, so the inner vertices it removes have none to lose.
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
//        No vertex is removed and the new split vertices start without sample segments.
        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override {}
        void fireMergePrefixPath(const RAGraphPath &path, Edge &new_edge) override;
        void fireFlipVertex(Vertex &core, const std::unordered_map<EdgeId, EdgeId> &extensions) override;
//        Samples lie within their vertex (kpomer chunks at least k + 1 long), are strictly ordered (starts and finishes
//        non-decreasing, never both equal), v.rc() holds exactly the mirror of v's samples, weight and raw_support are
//        the sums over the samples, and k-mers have none.
        bool checkConsistency(std::ostream &os) override;
    };
}
