//
// Created by Andrey Bzikadze on 11/25/21.
//

#include "mdbg_vertex_processor.hpp"

using namespace repeat_resolution;

void MDBGSimpleVertexProcessor::Process0In1Pout(MultiplexDBG &graph,
                                                const RRVertexType &vertex) {
    RRVertexProperty &v_prop = graph.node_prop(vertex);
    auto[out_nbr_begin, out_nbr_end] = graph.out_neighbors(vertex);
    std::vector<decltype(out_nbr_begin)> neighbors_its;
    for (auto it = out_nbr_begin; it!=out_nbr_end; ++it) {
        neighbors_its.emplace_back(it);
    }
    for (auto it : neighbors_its) {
    // for (auto it = out_nbr_begin; it!=out_nbr_end; ++it) {
        RRVertexType new_vertex = graph.GetNewVertex(v_prop.Seq());
        graph.MoveEdge(vertex, it, new_vertex, it->first);
        graph.IncreaseVertex(new_vertex, 1);
    }
    graph.remove_nodes(vertex); // careful: Iterator is invalidated
}

void MDBGSimpleVertexProcessor::Process1Pin0Out(MultiplexDBG &graph,
                                                const RRVertexType &vertex) {
    RRVertexProperty &v_prop = graph.node_prop(vertex);
    auto[in_nbr_begin, in_nbr_end] = graph.in_neighbors(vertex);
    std::vector<decltype(in_nbr_begin)> neighbors_its;
    for (auto it = in_nbr_begin; it!=in_nbr_end; ++it) {
        neighbors_its.emplace_back(it);
    }
    for (auto it : neighbors_its) {
    // for (auto it = in_nbr_begin; it!=in_nbr_end; ++it) {
        RREdgeIndexType edge_index = it->second.prop().Index();
        RRVertexType new_vertex = graph.GetNewVertex(v_prop.Seq());
        // need to construct a NeighborIterator pointing to vertex
        auto out_nbr = graph.out_neighbors(it->first).first;
        while (out_nbr->second.prop().Index()!=edge_index) {
            ++out_nbr;
        }
        graph.MoveEdge(it->first, out_nbr, it->first, new_vertex);
        graph.IncreaseVertex(new_vertex, 1);
    }
    graph.remove_nodes(vertex); // careful: Iterator is invalidated
}

void MDBGSimpleVertexProcessor::Process(MultiplexDBG &graph,
                                        const RRVertexType &vertex,
                                        uint64_t n_iter) {
    const int indegree = graph.count_in_neighbors(vertex);
    const int outdegree = graph.count_out_neighbors(vertex);
    VERIFY(indegree < 2 or outdegree < 2);

    VERIFY_MSG(indegree!=1 or outdegree!=1,
               "no vertexes on nonbranching paths allowed");
    RRVertexProperty &v_prop = graph.node_prop(vertex);
    if (indegree==0 and outdegree==0) {
        // Isolates should be skipped
    } else if (indegree==0 and outdegree==1) {
        // tip. Only increment length
        graph.IncreaseVertex(vertex, n_iter);
    } else if (indegree==1 and outdegree==0) {
        // tip. Only increment length
        graph.IncreaseVertex(vertex, n_iter);
    } else if (indegree==0 and outdegree > 1) {
        // "Starting" vertex
        VERIFY(n_iter==1);
        Process0In1Pout(graph, vertex);

    } else if (indegree > 1 and outdegree==0) {
        // "Finishing" vertex
        VERIFY(n_iter==1);
        Process1Pin0Out(graph, vertex);

    } else if (indegree==1 and outdegree > 1) {
        graph.IncreaseVertex(vertex, n_iter);

    } else if (indegree > 1 and outdegree==1) {
        graph.IncreaseVertex(vertex, n_iter);
    }
}

std::pair<std::unordered_map<RREdgeIndexType, RRVertexType>,
          std::vector<RRVertexType>>
MDBGComplexVertexProcessor::SplitVertex(MultiplexDBG &graph,
                                        const RRVertexType &vertex) {
    RRVertexProperty &v_prop = graph.node_prop(vertex);
    std::unordered_map<RREdgeIndexType, RRVertexType> edge2vertex;
    std::vector<RRVertexType> new_vertices;

    auto[in_nbr_begin, in_nbr_end] = graph.in_neighbors(vertex);
    std::vector<decltype(in_nbr_begin)> in_neighbors_its;
    for (auto it = in_nbr_begin; it!=in_nbr_end; ++it) {
        in_neighbors_its.emplace_back(it);
    }
    for (auto it : in_neighbors_its) {
        const RRVertexType &neighbor = it->first;
        const RREdgeIndexType edge_index = it->second.prop().Index();
        RRVertexType new_vertex = graph.GetNewVertex(v_prop.Seq());
        new_vertices.emplace_back(new_vertex);
        auto e_it = graph.out_neighbors(neighbor).first;
        while (e_it->second.prop().Index()!=edge_index) {
            ++e_it;
        }
        graph.MoveEdge(neighbor, e_it, neighbor, new_vertex);
        graph.IncreaseVertex(new_vertex, 1);
        edge2vertex.emplace(edge_index, neighbor);
    }

    auto[out_nbr_begin, out_nbr_end] = graph.out_neighbors(vertex);

    std::vector<decltype(out_nbr_begin)> out_neighbors_its;
    for (auto it = out_nbr_begin; it!=out_nbr_end; ++it) {
        out_neighbors_its.emplace_back(it);
    }
    for (auto it : out_neighbors_its) {
        const RREdgeIndexType edge_index = it->second.prop().Index();
        RRVertexType new_vertex = graph.GetNewVertex(v_prop.Seq());
        new_vertices.emplace_back(new_vertex);
        graph.MoveEdge(vertex, it, new_vertex, it->first);
        graph.IncreaseVertex(new_vertex, 1);
        // here we use map[key] = value instead of map.emplace(key, value)
        // because edge_index is already saved in the dict if the edge is a loop
        edge2vertex[edge_index] = new_vertex;
    }
    return {edge2vertex, new_vertices};
}

void MDBGComplexVertexProcessor::Process(MultiplexDBG &graph,
                                         const RRVertexType &vertex,
                                         std::set<Sequence> &merged_self_loops) {
    const RRVertexProperty &v_prop = graph.node_prop(vertex);

    auto[ac_s2e, ac_e2s] = graph.GetEdgepairsVertex(vertex);

    /*
    {
        auto[in_edges, out_edges] = graph.GetNeighborEdgesIndexes(vertex);
        for (const RREdgeIndexType &edge : in_edges) {
            VERIFY(ac_s2e.find(edge) != ac_s2e.end())
        }
        for (const RREdgeIndexType &edge : out_edges) {
            VERIFY(ac_e2s.find(edge) != ac_e2s.end())
        }
    }
     */

    auto[edge2vertex, new_vertices] = SplitVertex(graph, vertex);

    for (const auto &[edge1, edge1_neighbors] : ac_s2e) {
        for (const auto &edge2 : edge1_neighbors) {
            // const RREdgeIndexType edge1 = FindMergeEdgeId(edge1_);
            const RRVertexType left_vertex = edge2vertex.at(edge1);
            auto e1_it = graph.FindOutEdgeIterator(left_vertex, edge1);

            // const RREdgeIndexType edge2 = FindMergeEdgeId(edge2_);
            const RRVertexType right_vertex = edge2vertex.at(edge2);
            auto e2_it = graph.FindOutEdgeIterator(right_vertex, edge2);

            graph.AddConnectingEdge(e1_it, right_vertex, e2_it);
        }
    }

    for (const RRVertexType &new_vertex : new_vertices) {
        const uint64_t indegree = graph.count_in_neighbors(new_vertex);
        const uint64_t outdegree = graph.count_out_neighbors(new_vertex);
        if (indegree==1 and outdegree==1) {
            auto in_rev_it = graph.in_neighbors(new_vertex).first;
            const RRVertexType &left_vertex = in_rev_it->first;
            auto out_rev_it = graph.out_neighbors(new_vertex).first;
            const RRVertexType &right_vertex = out_rev_it->first;

            if (left_vertex==new_vertex) {
                // self-loop should be skipped
                Sequence seq = graph.node_prop(new_vertex).Seq().ToSequence();
                merged_self_loops.emplace(std::move(seq));
                continue;
            }
            {
                RRVertexType v = right_vertex;
                while (v!=new_vertex and graph.count_out_neighbors(v)==1
                    and graph.count_in_neighbors(v)==1) {
                    v = graph.out_neighbors(v).first->first;
                }
                if (v==new_vertex) {
                    // cycle
                    Sequence s = graph.node_prop(new_vertex).Seq().ToSequence();
                    if (merged_self_loops.count(!s)) {
                        // rev-compl of this vertex in graph — not merge
                        continue;
                    }
                }
            }

            auto in_it = graph.FindOutEdgeIterator(left_vertex,
                                                   in_rev_it->second.prop()
                                                       .Index());
            auto out_it = graph.out_neighbors(new_vertex).first;
            graph.MergeEdges(left_vertex, in_it, out_it);
        }
    }
    graph.remove_nodes(vertex);
}