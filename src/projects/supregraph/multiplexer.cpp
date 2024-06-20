#include <assembly_graph/ag_algorithms.hpp>
#include "multiplexer.hpp"

using namespace spg;

Multiplexer::Multiplexer(SupreGraph &graph, ag::AlignedReadStorage<SPGTraits> &reads, DecisionRule &rule, size_t max_core_length) :
                            graph(graph), reads(reads), rule(rule), max_core_length(max_core_length) {
    for(Vertex &v : graph.verticesUnique()) {
        VERIFY(v.isCanonical());
        if(v.isCore() && v.outDeg() > 0 && v.inDeg() > 0) {
            core_queue.insert(v.getId());
        }
    }
}

VertexResolutionResult Multiplexer::multiplex(logging::Logger &logger, size_t threads, Vertex &vertex) {
//            VERIFY(vertex.inDeg() > 1 && vertex.outDeg() > 1);
    core_queue.erase(vertex.getId());
    if(vertex.marked())
        return {vertex};
    if(vertex.inDeg() == 1 && vertex.outDeg() == 1 &&
                (!vertex.front().getFinish().isJunction() || !vertex.rc().front().getFinish().isJunction())) {
        GraphPath path = ag::PathHelper<SPGTraits>::WalkForward(vertex.front());
        if(path.getFinish() != vertex) {
            path = ag::PathHelper<SPGTraits>::WalkForward(vertex.rc().front()).RC() + path;
        }
        Vertex &res = path.getStart().isJunction() ? graph.mergePath(path) : graph.mergeLoop(path);
        return {res};
    }
    logger.trace() << "Multiplexing step " << vertex.getId() << std::endl;
    VertexResolutionPlan rr = rule.judge(vertex);
    logger.trace() << "Judgement: " << rr << std::endl;
    if(!rr.empty()) {
        logger.trace() << "Starting to resolve" << std::endl;
        for(auto it : rr.connectionsUnique()) {
            logger.trace() << it.incoming().getId() << " " << it.outgoing().getId() << std::endl;
        }
        VertexResolutionResult res = graph.resolveVertex(vertex, rr);
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
            if(vid->marked())
                continue;
            if(vid->isCanonical() && vid->isCore() && vid->outDeg() > 0 && vid->inDeg() > 0 && vid->size() < max_core_length &&
               (vid->isJunction() || ag::PathHelper<SPGTraits>::WalkForward(vid->front()).getFinish().isJunction()))
                core_queue.insert(vid);
        }
        logger.trace() << "Result: " << res << std::endl;
        return res;
    } else {
        logger.trace() << "Skipped" << std::endl;
        return {vertex};
    }
}

VertexResolutionResult Multiplexer::multiplex(logging::Logger &logger, size_t threads) {
    Vertex &next = **core_queue.begin();
    core_queue.erase(core_queue.begin());
    return multiplex(logger, threads, next);
}

void Multiplexer::fullMultiplex(logging::Logger &logger, size_t threads) {
    while(!finished()) {
        multiplex(logger, threads);
    }
    graph.removeMarked();
}
