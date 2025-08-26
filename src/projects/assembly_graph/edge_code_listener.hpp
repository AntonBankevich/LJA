#pragma once
#include "graph_listeners.hpp"

namespace ag {
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
        void fireAddSupreVertex(Vertex &v, Edge &e) override;
    };
}