#pragma once

#include "vertex_resolution.hpp"
#include "supregraph.hpp"
#include "read_storage.hpp"
#include <unordered_set>
#include <common/logging.hpp>

namespace spg {

    class Multiplexer {
    private:
        SupreGraph &graph;
        DecisionRule &rule;
        size_t max_core_length;
        std::unordered_set<VertexId> core_queue;//Store only canonical vertices
    public:
        Multiplexer(SupreGraph &graph, DecisionRule &rule, size_t max_core_length) : graph(graph), rule(rule), max_core_length(max_core_length) {
            for(Vertex &v : graph.verticesUnique()) {
                VERIFY(v.isCanonical());
                if(v.isCore() && v.outDeg() > 0 && v.inDeg() > 0) {
                    core_queue.insert(v.getId());
                }
            }
        }
        Multiplexer(Multiplexer &&) = delete;
        Multiplexer(Multiplexer &) = delete;

        VertexResolutionResult multiplex(logging::Logger &logger, size_t threads, Vertex &vertex, PathIndex &index) {
//            VERIFY(vertex.inDeg() > 1 && vertex.outDeg() > 1);
            core_queue.erase(vertex.getId());
            std::cout << "Multiplexing step " << vertex.getId() << std::endl;
            VertexResolutionPlan rr = rule.judge(vertex);
            std::cout << "Judgement: " << rr << std::endl;
            if(!rr.empty()) {
                std::cout << "Starting to resolve" << std::endl;
                for(auto it : rr.connectionsUnique()) {
                    std::cout << it.first << " " << it.second << std::endl;
                }
                VERIFY(index.checkReadIndexConsistency());
                VertexResolutionResult res = graph.resolveVertex(vertex, rr);
                VERIFY(index.checkReadIndexConsistency());
                std::vector<VertexId> candidates;
                for(Vertex &new_vertex : res.newVertices()) {
                    for(Vertex &v : ag::ThisAndRC(new_vertex)) {
                        candidates.emplace_back(v.getId());
                        for (Edge &edge: v) {
                            candidates.emplace_back(edge.getFinish().getId());
                            candidates.emplace_back(edge.getFinish().rc().getId());
                        }
                    }
                }
                for(VertexId vid : candidates) {
                    if(vid->isCanonical() && vid->isCore() && vid->outDeg() > 0 && vid->inDeg() > 0 && vid->size() < max_core_length &&
                            (vid->isJunction() || GraphPath::WalkForward(vid->front()).finish().isJunction()))
                        core_queue.insert(vid);
                }
                std::cout << "Result: " << res << std::endl;
                return res;
            } else {
                std::cout << "Skipped" << std::endl;
                return {vertex};
            }
        }

        VertexResolutionResult multiplex(logging::Logger &logger, size_t threads, PathIndex &index) {
            Vertex &next = **core_queue.begin();
            core_queue.erase(core_queue.begin());
            return multiplex(logger, threads, next, index);
        }

        bool finished() const {
            return core_queue.empty();
        }

        //        Graph should be in proper form. No outer edges, no unbranching vertices.
        void fullMultiplex(logging::Logger &logger, size_t threads, PathIndex &index) {
            while(!finished()) {
                multiplex(logger, threads, index);
            }
            graph.removeMarked();
        }
    };
}