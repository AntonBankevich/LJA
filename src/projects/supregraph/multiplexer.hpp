#pragma once
#include "supregraph.hpp"
namespace spg {
    class DecisionRule {
    public:
        virtual std::vector<std::pair<EdgeId, EdgeId>> judge(Vertex &v) = 0;
        virtual ~DecisionRule() = default;
    };

    class RandomDecisionRule : public DecisionRule {
    public:
        virtual std::vector<std::pair<EdgeId, EdgeId>> judge(Vertex &v) override {
            std::vector<std::pair<EdgeId, EdgeId>> res;
            auto out_it = v.begin();
            auto inc = v.incoming();
            auto in_it = inc.begin();
            while(out_it != v.end() || in_it != inc.end()) {
                if(out_it == v.end()) --out_it;
                if(in_it == inc.end()) --in_it;
                res.emplace_back(in_it->getId(), out_it->getId());
                res.emplace_back(out_it->rc().getId(), in_it->rc().getId());
                ++in_it;
                ++out_it;
            }
            std::sort(res.begin(), res.end());
            res.erase(std::unique(res.begin(), res.end()), res.end());
            return std::move(res);
        }
    };

    class Multiplexer {
    private:
        SupreGraph &graph;
        DecisionRule &rule;
        std::unordered_set<VertexId> core_queue;//Store only canonical vertices
    public:
        Multiplexer(SupreGraph &graph, DecisionRule &rule) : graph(graph), rule(rule) {
            for(Vertex &v : graph.verticesUnique()) {
                if(v.isCanonical() && v.isCore() && v.outDeg() > 0 && v.inDeg() > 0) {
                    core_queue.insert(v.getId());
                }
            }
        }
        Multiplexer(Multiplexer &&) = delete;
        Multiplexer(Multiplexer &) = delete;

        std::vector<VertexId> multiplex(Vertex &vertex) {
//            VERIFY(vertex.inDeg() > 1 && vertex.outDeg() > 1);
            core_queue.erase(vertex.getId());
            std::vector<std::pair<EdgeId, EdgeId>> rr = rule.judge(vertex);
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
            return multiplex(**core_queue.begin());
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