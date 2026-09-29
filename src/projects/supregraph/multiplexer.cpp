#include <assembly_graph/ag_algorithms.hpp>
#include "multiplexer.hpp"

#include "spoa/include/spoa/graph.hpp"

using namespace spg;

Multiplexer::Multiplexer(ag::AssemblyGraph &graph, ag::AlignedReadStorage &reads, DecisionRule &rule, size_t max_core_length) :
                            graph(graph), reads(reads), rule(rule), max_core_length(max_core_length) {
    reset();
}
// std::vector<VertexId> Multiplexer::flipRepeat(logging::Logger &logger, size_t threads, Vertex &vertex) {
//     VERIFY(vertex.isOuter());
//     if (vertex.incFrontVertex().outDeg() == 1 && vertex.frontVertex().inDeg() == 1 &&
//         vertex.incFrontVertex().inDeg() != 1 && vertex.frontVertex().outDeg() != 1) {
//         if (rule.judgeFlip(vertex)) {
//             std::vector<VertexId> res;
//             ag::VertexResolutionResult res1 = graph.resolveVertex(vertex, VertexResolutionPlan::SimplePlan(vertex.incFrontVertex()));
//             for (Vertex &new_vertex: res1.newVertices()) {res.emplace_back(new_vertex.getId());}
//             if (vertex != vertex.rc()) {
//                 ag::VertexResolutionResult res2 = graph.resolveVertex(vertex, VertexResolutionPlan::SimplePlan(vertex.frontVertex()));
//                 for (Vertex &new_vertex: res2.newVertices()) {res.emplace_back(new_vertex.getId());}
//             }
//             pushCore(vertex);
//             return res;
//         }
//     } else {
//         pushCore(vertex.incFrontVertex());
//         pushCore(vertex.frontVertex());
//     }
//     return {};
// }

std::vector<VertexId> Multiplexer::multiplex(logging::Logger &logger, size_t threads, Vertex &vertex) {
//            VERIFY(vertex.inDeg() > 1 && vertex.outDeg() > 1);
    if (vertex.marked())
        return {};
    logger.trace() << "Premultiplexing " << vertex.getId() << " " << vertex.inDeg() << " " << vertex.outDeg() << " " << vertex.isCore() << std::endl;
    if(!vertex.isCore()) {
        if (vertex.isOuter()) {
            ag::GraphPath path = vertex;
            path.push_front(vertex.incFront());
            path += vertex.front();
            if (path.getStart().inDeg() > 1 && path.getStart().outDeg() == 1 && path.getFinish().outDeg() > 1 && path.getFinish().inDeg() == 1) {
                if (rule.judgeFlip(path)) {
                    graph.resolveVertex(path.getStart(), VertexResolutionPlan::SimplePlan(path.getStart()));
                    if (vertex != vertex.rc())
                        graph.resolveVertex(path.getFinish(), VertexResolutionPlan::SimplePlan(path.getFinish()));
                }
            }
        }
    }
    if (!vertex.isCore())
        return {};
    logger.trace() << "Multiplexing vertex " << vertex.getId() << std::endl;
    VertexResolutionPlan rr = rule.judge(vertex);
    logger.trace() << "Judgement: " << rr << std::endl;
    if (vertex.inDeg() + vertex.outDeg() == 1 || !rr.empty()) {
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
    ag::GraphPath path = ag::PathHelper::WalkForward(vertex.front());
    ag::GraphPath rpath = ag::PathHelper::WalkForward(vertex.rc().front()).RC();
    if (path.getFinish().isJunction()) {
        if (rpath.getStart().isJunction()) {
            path = rpath + path;
        } else {
            VERIFY(rpath.getStart() == rpath.getFinish().rc());
            path = path.RC() + rpath + path;
        }
    } else {
        if (path.getFinish() == path.getStart().rc()) {
            path = rpath + path;
        } else {
            VERIFY(path.getStart() == path.getFinish());
        }
    }
    // logger.trace() << "Push " << path.getStart().getId() << " " << path.getFinish().getId() << std::endl;
    // pushCore(path.getStart());
    // pushCore(path.getFinish());
    if(path.calculateSize() > 2 || (path.calculateSize() == 2 && (!path.frontEdge().isPrefix() || !path.backEdge().isSuffix()))) {
        logger.trace() << "Merging path " << path.str() << std::endl;
        VERIFY(!path.frontEdge().isSuffix());
        VERIFY(!path.backEdge().isPrefix());
        // TODO: switch to Supregraph and move this functionality to it!!!
        Vertex &res = ag::MergePathOrLoop(path, graph);
        logger.trace() << "Push outer" << res << std::endl;
        pushCore(res);
        // Vertex &res = path.getStart().isJunction() ? graph.mergePath(path) : graph.mergeLoop(path);
        return {res.getId()};
    } else {
        logger.trace() << "Push core " << vertex << std::endl;
        pushCore(vertex);
        return {};
    }
}

std::vector<VertexId> Multiplexer::process(logging::Logger &logger, size_t threads) {
    if (!merge_queue.empty()) {
        Vertex &next = *merge_queue.back();
        merge_queue.pop_back();
        return merge(logger, threads, next);
    }
    if (!core_queue.empty()) {
        return multiplex(logger, threads, popCore());
    }
    VERIFY(false);
    return {};
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

void Multiplexer::reset() {
    core_queue.clear();
    merge_queue.clear();
    for(Vertex &v : graph.verticesUnique()) {
        VERIFY(v.isCanonical());
        if(v.isCore()) {
            pushCore(v);
        }
    }
}

void Multiplexer::pushCore(Vertex &vertex) {core_queue.emplace(vertex.size(), vertex.getCanonical().getId());}
