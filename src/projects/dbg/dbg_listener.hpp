#pragma once

#include "assembly_graph/graph_listeners.hpp"

namespace dbg {
//    Base for listeners that only work on a de Bruijn graph: every vertex is a k-mer. Operations that create longer
//    vertices (the supregraph ones) are prohibited: they abort and cannot be overridden. Such listeners must be
//    detached before the graph is converted to a supregraph.
    class DBGListener : public ag::ResolutionListener {
    public:
        using ag::ResolutionListener::ResolutionListener;

        void fireEdgeToSupreVertex(ag::Vertex &v, ag::Edge &e) final {VERIFY(false);}
        void fireReplacePath(const ag::PathReplacement &replacement) final {VERIFY(false);}
        void fireMergeLoop(const ag::GraphPath &path, ag::Vertex &new_vertex) final {VERIFY(false);}
        void fireResolveVertex(ag::Vertex &core, const ag::VertexResolutionResult &resolution) final {VERIFY(false);}
    };
}
