//
// Created by Andrey Bzikadze on 11/25/21.
//

#include "mdbg_inc.hpp"

using namespace repeat_resolution;

void MultiplexDBGIncreaser::process_vertex(MultiplexDBG &graph,
                                           const RRVertexType &vertex) {
  if (graph.node_prop(vertex).frozen) {
    return;
  }
  const int indegree = graph.count_in_neighbors(vertex);
  const int outdegree = graph.count_out_neighbors(vertex);
  if (indegree >= 2 and outdegree >= 2) {
    complex_vertex_processor.process(graph, vertex);
  } else {
    simple_vertex_processor.process(graph, vertex);
  }
}

void MultiplexDBGIncreaser::collapse_edge(
    MultiplexDBG &graph, MultiplexDBG::ConstIterator s_it,
    MultiplexDBG::NeighborsIterator e_it) {
  RRVertexType s = *s_it;
  RRVertexType e = e_it->first;
  VERIFY(s != e);
  VERIFY(graph.count_out_neighbors(s_it) == 1);
  VERIFY(graph.count_in_neighbors(e_it->first) == 1);

  RREdgeProperty &edge_prop = e_it->second.prop();
  graph.rr_paths->remove(edge_prop.get_index());

  if (graph.count_in_neighbors(s) == 0 and graph.count_out_neighbors(e) == 0) {
    // isolated vertex. Need to freeze and save its label
    graph.isolate_properties.emplace(s, std::move(edge_prop));
    graph.freeze_vertex(s);
  }

  graph.remove_edge(s_it, e_it);

  auto [out_nbr_begin, out_nbr_end] = graph.out_neighbors(e);
  for (auto out_nbr_it = out_nbr_begin; out_nbr_it != out_nbr_end;
       ++out_nbr_it) {
    graph.move_edge(e, out_nbr_it, s, out_nbr_it->first);
  }
  VERIFY(graph.count_in_neighbors(e) == 0 and
         graph.count_out_neighbors(e) == 0);
  graph.remove_nodes(e);
}

void MultiplexDBGIncreaser::collapse_short_edges_into_vertices(
    MultiplexDBG &graph) {
  for (const RRVertexType &v1 : graph) {
    const RRVertexProperty &v1p = graph.node_prop(v1);
    if (graph.count_out_neighbors(v1) == 0) {
      continue;
    }
    auto [out_it_begin, out_it_end] = graph.out_neighbors(v1);
    std::vector<EdgeIndexType> edges2collapse;
    for (auto it = out_it_begin; it != out_it_end; ++it) {
      const RRVertexType &v2 = it->first;
      const RRVertexProperty &v2p = graph.node_prop(v2);
      const RREdgeProperty &edge_property = it->second.prop();
      if (v1p.len == edge_property.size() or v2p.len == edge_property.size()) {
        VERIFY(v1p.len == v2p.len);
        VERIFY(not v1p.frozen and not v2p.frozen);
        edges2collapse.push_back(edge_property.get_index());
      }
    }
    for (const EdgeIndexType &edge_index : edges2collapse) {
      // iterator might be getting invalidated every time we collapse an edge
      // thus, we find the iterator for every edge from scratch
      auto it = [&graph, &v1, &edge_index]() {
        auto it = graph.out_neighbors(v1).first;
        while (it->second.prop().get_index() != edge_index) {
          ++it;
        }
        return it;
      }();
      collapse_edge(graph, graph.find(v1), it);
    }
  }
}

MultiplexDBGIncreaser::MultiplexDBGIncreaser(const uint64_t start_k,
                                             const uint64_t saturating_k,
                                             logging::Logger &logger,
                                             const bool debug)
    : start_k{start_k}, saturating_k{saturating_k}, logger{logger}, debug{
                                                                        debug} {
  VERIFY(saturating_k >= start_k);
}

void MultiplexDBGIncreaser::Increment(MultiplexDBG &graph) {
  if (graph.is_frozen()) {
    logger.info() << "Graph is frozen, no increase of k possible" << std::endl;
    return;
  }

  if (start_k + graph.n_iter == saturating_k) {
    logger.info() << "K is saturated, no increase of k possible" << std::endl;
    return;
  }

  // since iterators over vertexes might invalidate, first save the vertexes
  const std::vector<RRVertexType> vertexes = [&graph]() {
    std::vector<RRVertexType> vertexes;
    for (auto &v : graph) {
      auto &vertex = (RRVertexType &)v;
      vertexes.emplace_back(v);
    }
    return vertexes;
  }();

  for (const auto &vertex : vertexes) {
    process_vertex(graph, vertex);
  }
  ++graph.n_iter;

  collapse_short_edges_into_vertices(graph);
  graph.freeze_unpaired_vertices();

  if (debug) {
    graph.assert_validity();
  }
}

void MultiplexDBGIncreaser::IncrementN(MultiplexDBG &graph, uint64_t N) {
  N = std::min(N, saturating_k - start_k - graph.n_iter);
  for (int i = 0; i < N; ++i) {
    std::cout << i << "\n";
    Increment(graph);
  }
}
