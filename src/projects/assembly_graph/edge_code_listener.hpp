#pragma once
#include "graph_listeners.hpp"

namespace ag {
    template<class Traits>
    class EdgeCodeListener : public ResolutionListener<Traits> {
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Vertex::VertexId VertexId;
        EdgeCodeListener(ResolutionFire<Traits> &graph) : ResolutionListener<Traits>(graph, "EdgeCodeListener") {}

        void fireMergePath(const std::vector<EdgeId> &path, Vertex &new_vertex) override {
            VERIFY(new_vertex.inDeg() == 1);
            SequenceBuilder sb;
            for(EdgeId eid : path) sb.append(eid->edge_code);
            new_vertex.rc().front().rc().edge_code = sb.BuildSequence();
        };

        void fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) override {
            SequenceBuilder sb;
            for(EdgeId eid : path) sb.append(eid->edge_code);
            new_edge.edge_code = sb.BuildSequence();
            VERIFY(new_edge.edge_code.startsWith(path.front()->edge_code));
        };

        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &, const AlignmentForm &) override {
            new_edge.edge_code = left.edge_code + right.edge_code;
        }

        //        TODO: implement properly
        void fireMergeLoop(const GraphPath<Traits> &path, Vertex &new_vertex) override {VERIFY(false)};

        void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) override {
            for(EdgeId eid : split) {
                if(eid != split.front())
                    eid->edge_code = eid->firstNucl() + edge.edge_code;
                else
                    eid->edge_code = edge.edge_code;
                if(eid != split.back())
                    eid->rc().edge_code = eid->rc().firstNucl() + edge.rc().edge_code;
                else
                    eid->rc().edge_code = edge.rc().edge_code;
            }
        };

        void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) override {
            for(Vertex &new_vertex : resolution.newVertices()) {
                new_vertex.rc().front().rc().edge_code = resolution.get(new_vertex).outgoing().edge_code;
            }
        };

        void fireAddSupreVertex(Vertex &v, Edge &e) override {
            v.front().rc().edge_code = e.rc().edge_code;
            v.rc().front().rc().edge_code = e.edge_code;
        }
    };
}