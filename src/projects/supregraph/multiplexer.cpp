#include <assembly_graph/ag_algorithms.hpp>
#include "multiplexer.hpp"

#include "spoa/include/spoa/graph.hpp"

using namespace spg;

Multiplexer::Multiplexer(ag::AssemblyGraph &graph, ag::AlignedReadStorage &reads, DecisionRule &rule, size_t max_core_length) :
                            graph(graph), reads(reads), rule(rule), max_core_length(max_core_length) {
    for(Vertex &v : graph.verticesUnique()) {
        VERIFY(v.isCanonical());
        if(v.isCore() && v.outDeg() > 0 && v.inDeg() > 0) {
            pushCore(v);
        }
    }
}

std::vector<VertexId> Multiplexer::multiplex(logging::Logger &logger, size_t threads, Vertex &vertex) {
//            VERIFY(vertex.inDeg() > 1 && vertex.outDeg() > 1);
            logger.trace() << "Premultiplexing " << vertex.getId() << " " << vertex.inDeg() << " " << vertex.outDeg() << " " << vertex.isCore() << std::endl;
    if(!vertex.isCore())
        return {vertex.getId()};
    logger.trace() << "Multiplexing vertex " << vertex.getId() << std::endl;
    VertexResolutionPlan rr = rule.judge(vertex);
    logger.trace() << "Judgement: " << rr << std::endl;
    if (!rr.empty()) {
        logger.trace() << "Starting to resolve" << std::endl;
        for (auto it: rr.connectionsUnique()) {
            logger.trace() << it.incoming().getId() << " " << it.outgoing().getId() << std::endl;
        }
        ag::VertexResolutionResult vrres = graph.resolveVertex(vertex, rr);
        std::vector<VertexId> candidates;
        std::vector<VertexId> res;
        for (Vertex &new_vertex: vrres.newVertices()) {
            res.emplace_back(new_vertex.getId());
	    logger.trace() << "Push merge " << new_vertex.getId() << std::endl;
            merge_queue.emplace_back(new_vertex.getId());
        }
        logger.trace() << "Result: " << res << std::endl;
        return std::move(res);
    } else {
        logger.trace() << "Could not resolve" << std::endl;
        return {};
    }
}

std::vector<VertexId> Multiplexer::merge(logging::Logger &logger, size_t threads, Vertex &vertex) {
//            VERIFY(vertex.inDeg() > 1 && vertex.outDeg() > 1);
    if(vertex.marked())
        return {};
    VERIFY(!vertex.isJunction());
    logger.trace() << "Processing vertex " << vertex.getId() << std::endl;
    Edge &start = ag::PathHelper::WalkForward(vertex.rc().front()).backEdge().rc();
    ag::GraphPath path = ag::PathHelper::WalkForward(start);
    if(!path.getStart().isJunction() && path.getStart() != path.getFinish()) {
        VERIFY(path.getStart() == path.getFinish().rc());
        path = path + ag::PathHelper::WalkForward(path.getFinish().front());
    } else {
        logger.trace() << "Push " << path.getStart().getId() << " " << path.getFinish().getId() << std::endl;
        pushCore(path.getStart());
        pushCore(path.getFinish());
    }
    if(path.calculateSize() > 2 || (path.calculateSize() == 2 && (!path.frontEdge().isPrefix() || !path.backEdge().isSuffix()))) {
        logger.trace() << "Merging path " << path.str() << std::endl;
        // TODO: switch to Supregraph and move this functionality to it!!!
        Vertex &res = ag::MergePathSPG(path, graph);
        // Vertex &res = path.getStart().isJunction() ? graph.mergePath(path) : graph.mergeLoop(path);
        return {res.getId()};
    } else {
        return {};
    }
}

std::vector<VertexId> Multiplexer::process(logging::Logger &logger, size_t threads) {
    if(merge_queue.empty()) {
        return multiplex(logger, threads, popCore());
    } else {
        Vertex &next = *merge_queue.back();
        merge_queue.pop_back();
        return merge(logger, threads, next);
    }
}

void Multiplexer::fullMultiplex(logging::Logger &logger, size_t threads) {
    while(!finished()) {
        process(logger, threads);
    }
    graph.removeMarked();
}

Vertex &Multiplexer::popCore() {
    Vertex &res = *core_queue.begin()->second;
    core_queue.erase(core_queue.begin());
    return res;
}

void Multiplexer::pushCore(Vertex &vertex) {core_queue.emplace(vertex.size(), vertex.getCanonical().getId());}
