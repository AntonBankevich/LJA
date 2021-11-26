//
// Created by Andrey Bzikadze on 11/25/21.
//

#include "mdbg_vertex_processor.hpp"

using namespace repeat_resolution;

void MDBGSimpleVertexProcessor::process_0in_1pout(MultiplexDBG &graph,
                                                  const RRVertexType &vertex) {
  RRVertexProperty &v_prop = graph.node_prop(vertex);
  auto [out_nbr_begin, out_nbr_end] = graph.out_neighbors(vertex);
  for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
    RRVertexType new_vertex = graph.get_new_vertex(v_prop.len + 1);
    graph.move_edge(vertex, it, new_vertex, it->first);
  }
  graph.remove_nodes(vertex); // careful: Iterator is invalidated
}
void MDBGSimpleVertexProcessor::process_1pin_0out(MultiplexDBG &graph,
                                                  const RRVertexType &vertex) {
  RRVertexProperty &v_prop = graph.node_prop(vertex);
  auto [in_nbr_begin, in_nbr_end] = graph.in_neighbors(vertex);
  for (auto it = in_nbr_begin; it != in_nbr_end; ++it) {
    RRVertexType new_vertex = graph.get_new_vertex(v_prop.len + 1);
    // need to construct a NeighborIterator pointing to vertex
    auto out_nbr = graph.out_neighbors(it->first).first;
    while (out_nbr->first != vertex) {
      ++out_nbr;
    }
    graph.move_edge(it->first, out_nbr, it->first, new_vertex);
  }
  graph.remove_nodes(vertex); // careful: Iterator is invalidated
}
void MDBGSimpleVertexProcessor::process_1in_1pout(MultiplexDBG &graph,
                                                  const RRVertexType &vertex) {
  RRVertexProperty &v_prop = graph.node_prop(vertex);
  auto in_nbr_begin = graph.in_neighbors(vertex).first;
  RREdgeProperty &in_edge = in_nbr_begin->second.prop();
  auto [out_nbr_begin, out_nbr_end] = graph.out_neighbors(vertex);
  for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
    RREdgeProperty &out_edge = it->second.prop();
    out_edge.prepend(in_edge, v_prop.len);
  }
  ++v_prop.len;
}
void MDBGSimpleVertexProcessor::process_1pin_1out(MultiplexDBG &graph,
                                                  const RRVertexType &vertex) {
  RRVertexProperty &v_prop = graph.node_prop(vertex);
  auto out_nbr_begin = graph.out_neighbors(vertex).first;
  RREdgeProperty &out_edge = out_nbr_begin->second.prop();
  auto [in_nbr_begin, in_nbr_end] = graph.in_neighbors(vertex);
  for (auto it = in_nbr_begin; it != in_nbr_end; ++it) {
    RREdgeProperty &in_edge = it->second.prop();
    in_edge.append(out_edge, v_prop.len);
  }
  ++v_prop.len;
}

void MDBGSimpleVertexProcessor::process(MultiplexDBG &graph,
                                        const RRVertexType &vertex) {
  const int indegree = graph.count_in_neighbors(vertex);
  const int outdegree = graph.count_out_neighbors(vertex);
  VERIFY(indegree < 2 or outdegree < 2);

  VERIFY_MSG(indegree != 1 or outdegree != 1,
             "no vertexes on nonbranching paths allowed");
  RRVertexProperty &v_prop = graph.node_prop(vertex);
  if (indegree == 0 and outdegree == 0) {
    // Isolates should be skipped
  } else if (indegree == 0 and outdegree == 1) {
    // tip. Only increment length
    ++v_prop.len;
  } else if (indegree == 1 and outdegree == 0) {
    // tip. Only increment length
    ++v_prop.len;

  } else if (indegree == 0 and outdegree > 1) {
    // "Starting" vertex
    process_0in_1pout(graph, vertex);

  } else if (indegree > 1 and outdegree == 0) {
    // "Finishing" vertex
    process_1pin_0out(graph, vertex);

  } else if (indegree == 1 and outdegree > 1) {
    process_1in_1pout(graph, vertex);

  } else if (indegree > 1 and outdegree == 1) {
    process_1pin_1out(graph, vertex);
  }
}

std::unordered_map<EdgeIndexType, RRVertexType>
MDBGComplexVertexProcessor::split_vertex(MultiplexDBG &graph,
                                         const RRVertexType &vertex) {
  RRVertexProperty &v_prop = graph.node_prop(vertex);
  std::unordered_map<EdgeIndexType, RRVertexType> edge2vertex;
  auto [in_nbr_begin, in_nbr_end] = graph.in_neighbors(vertex);
  for (auto it = in_nbr_begin; it != in_nbr_end; ++it) {
    const RRVertexType &neighbor = it->first;
    const EdgeIndexType edge_index = it->second.prop().get_index();
    RRVertexType new_vertex = graph.get_new_vertex(v_prop.len + 1);
    auto e_it = graph.out_neighbors(neighbor).first;
    while (e_it->second.prop().get_index() != edge_index) {
      ++e_it;
    }
    graph.move_edge(neighbor, e_it, neighbor, new_vertex);
    edge2vertex[edge_index] = neighbor;
  }

  auto [out_nbr_begin, out_nbr_end] = graph.out_neighbors(vertex);
  for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
    const EdgeIndexType edge_index = it->second.prop().get_index();
    RRVertexType new_vertex = graph.get_new_vertex(v_prop.len + 1);
    graph.move_edge(vertex, it, new_vertex, it->first);
    edge2vertex[edge_index] = new_vertex;
  }
  return edge2vertex;
}

void MDBGComplexVertexProcessor::process_11(
    MultiplexDBG &graph, const RRVertexType &vertex,
    const RRVertexType &left_vertex, const RRVertexType &right_vertex,
    MultiplexDBG::NeighborsIterator e1_it,
    MultiplexDBG::NeighborsIterator e2_it, const EdgeIndexType &edge1,
    const EdgeIndexType &edge2,
    std::unordered_map<EdgeIndexType, EdgeIndexType> &where_edge_merged) {
  if (edge1 != edge2) {
    graph.merge_edges(left_vertex, e1_it, right_vertex, e2_it,
                      graph.node_prop(vertex).len);
    where_edge_merged.emplace(edge2, edge1);
  } else {
    // isolated loop
    VERIFY(left_vertex == right_vertex);
    RRVertexType vertex2remove = e1_it->first;
    graph.move_edge(left_vertex, e1_it, left_vertex, left_vertex);
    graph.remove_nodes(vertex2remove);
    --graph.node_prop(left_vertex).len;
  }
}

void MDBGComplexVertexProcessor::process_not11(
    MultiplexDBG &graph, const RRVertexType &vertex,
    const RRVertexType &left_vertex, const RRVertexType &right_vertex,
    MultiplexDBG::NeighborsIterator e1_it,
    MultiplexDBG::NeighborsIterator e2_it, const EdgeIndexType &edge1,
    const EdgeIndexType &edge2,
    const std::unordered_set<EdgeIndexType> &edge1_neighbors,
    const std::unordered_set<EdgeIndexType> &edge2_neighbors,
    std::unordered_map<EdgeIndexType, EdgeIndexType> &where_edge_merged) {

  const EdgeIndexType new_index =
      graph.add_connecting_edge(e1_it, right_vertex, e2_it);

  if (edge1_neighbors.size() == 1 and edge2_neighbors.size() >= 2) {
    VERIFY(graph.count_out_neighbors(e1_it->first) == 1);
    auto new_edge_it = graph.out_neighbors(e1_it->first).first;
    graph.merge_edges(left_vertex, e1_it, e1_it->first, new_edge_it,
                      graph.node_prop(e1_it->first).len);

  } else if (edge1_neighbors.size() >= 2 and edge2_neighbors.size() == 1) {
    VERIFY(graph.count_in_neighbors(right_vertex) == 1);
    auto new_edge_it = graph.out_neighbors(e1_it->first).first;
    while (new_edge_it->second.prop().get_index() != new_index) {
      ++new_edge_it;
    }
    graph.merge_edges(e1_it->first, new_edge_it, right_vertex, e2_it,
                      graph.node_prop(right_vertex).len);
  }
}

void MDBGComplexVertexProcessor::process(MultiplexDBG &graph,
                                         const RRVertexType &vertex) {
  const RRVertexProperty &v_prop = graph.node_prop(vertex);

  auto [ac_s2e, ac_e2s] = graph.get_edgepairs_vertex(vertex);
  std::unordered_map<EdgeIndexType, RRVertexType> edge2vertex =
      split_vertex(graph, vertex);

  std::unordered_map<EdgeIndexType, EdgeIndexType> where_edge_merged;
  auto FindMergeEdgeId = [&where_edge_merged](const EdgeIndexType edge_ind_) {
    EdgeIndexType edge_ind{edge_ind_};
    while (where_edge_merged.find(edge_ind) != where_edge_merged.end()) {
      edge_ind = where_edge_merged.at(edge_ind);
    }
    return edge_ind;
  };

  for (const auto &[edge1_, edge1_neighbors] : ac_s2e) {
    for (const auto &edge2_ : edge1_neighbors) {
      const EdgeIndexType edge1 = FindMergeEdgeId(edge1_);
      const RRVertexType left_vertex = edge2vertex.at(edge1);
      auto e1_it = graph.find_edge_iterator(left_vertex, edge1);

      const EdgeIndexType edge2 = FindMergeEdgeId(edge2_);
      const RRVertexType right_vertex = edge2vertex.at(edge2);
      const std::unordered_set<EdgeIndexType> &edge2_neighbors = ac_e2s[edge2];
      auto e2_it = graph.find_edge_iterator(right_vertex, edge2);

      if (edge1_neighbors.size() == 1 and edge2_neighbors.size() == 1) {
        process_11(graph, vertex, left_vertex, right_vertex, e1_it, e2_it,
                   edge1, edge2, where_edge_merged);
      } else {
        process_not11(graph, vertex, left_vertex, right_vertex, e1_it, e2_it,
                      edge1, edge2, edge1_neighbors, edge2_neighbors,
                      where_edge_merged);
      }
    }
  }
  graph.remove_nodes(vertex);
}