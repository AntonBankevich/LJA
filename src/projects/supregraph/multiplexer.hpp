#pragma once

#include "vertex_resolution.hpp"
#include "supregraph.hpp"
#include <unordered_set>
namespace spg {

    class Multiplexer {
    private:
        SupreGraph &graph;
        DecisionRule &rule;
        std::unordered_set<VertexId> core_queue;//Store only canonical vertices
    public:
        Multiplexer(SupreGraph &graph, DecisionRule &rule) : graph(graph), rule(rule) {
            for(Vertex &v : graph.verticesUnique()) {
                VERIFY(v.isCanonical());
                if(v.isCore() && v.outDeg() > 0 && v.inDeg() > 0) {
                    core_queue.insert(v.getId());
                }
            }
        }
        Multiplexer(Multiplexer &&) = delete;
        Multiplexer(Multiplexer &) = delete;

        std::vector<VertexId> multiplex(Vertex &vertex) {
//            VERIFY(vertex.inDeg() > 1 && vertex.outDeg() > 1);
            core_queue.erase(vertex.getId());
            VertexResolutionPlan rr = rule.judge(vertex);
            if(!rr.empty()) {
                std::vector<VertexId> res = graph.resolveVertex(vertex, rr);
                std::vector<VertexId> candidates;
                for(VertexId vid : res) {
                    candidates.emplace_back(vid);
                    for(Edge &edge : *vid) {
                        candidates.emplace_back(edge.getFinish().getId());
                        candidates.emplace_back(edge.getFinish().rc().getId());
                    }
                }
                for(VertexId vid : candidates) {
                    if(vid->isCanonical() && vid->isCore() && vid->outDeg() > 0 && vid->inDeg() > 0)
                        core_queue.insert(vid);
                }
                return res;
            } else
                return {};
        }

        std::vector<VertexId> multiplex() {
            Vertex &next = **core_queue.begin();
            core_queue.erase(core_queue.begin());
            return multiplex(next);
        }

        bool finished() const {
            return core_queue.empty();
        }

        //        Graph should be in proper form. No outer edges, no unbranching vertices.
        void fullMultiplex() {
            while(!finished()) {
                multiplex();
            }
        }
    };
}