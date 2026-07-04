#pragma once
#include "graph_listeners.hpp"

namespace ag {
//    Maintains Edge::getCode(): on any merge, a new edge's code is defined as the exact concatenation of
//    the merged edges' codes (not recomputed independently). This is deliberate, not incidental — it is
//    what lets path/read listeners elsewhere (e.g. AlignedReadStorageMaintenance, SuffixTracker) treat a
//    path that merely passes through a merged run as needing zero updates: matching an outgoing edge by
//    its code's first character, and advancing a path position by a code's length, both behave identically
//    whether the run is still N separate edges or one concatenated edge.
    class EdgeCodeListener : public ResolutionListener {
    public:
        EdgeCodeListener(ResolutionFire &graph) : ResolutionListener(graph, "EdgeCodeListener") {}
        void fireAddEdge(Edge &new_edge) override;
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &, const AlignmentForm &) override;
        //        TODO: implement properly
        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) override;
        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override;;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;;
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
    };
}