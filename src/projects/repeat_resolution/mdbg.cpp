//
// Created by Andrey Bzikadze on 11/25/21.
//

#include "mdbg.hpp"

using namespace repeat_resolution;

void MultiplexDBG::freeze_isolated_loops() {
  for (const auto &vertex : *this) {
    if (count_in_neighbors(vertex) == 1 and count_out_neighbors(vertex) == 1) {
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
      VERIFY_MSG(in_nbr_begin->first == vertex,
                 "No 1in-1out vertices are allowed except loops")
      freeze_vertex(vertex);
    }
  }
}

void MultiplexDBG::assert_validity() const {
  int64_t est_max_vert_index = [this]() {
    int64_t est_max_vert_index{-1};
    for (const auto &vertex : *this) {
      est_max_vert_index = std::max(est_max_vert_index, (int64_t)vertex);
    }
    return est_max_vert_index;
  }();
  VERIFY(next_vert_index >= 1 + est_max_vert_index);

  int64_t est_max_edge_index = [this]() {
    int64_t est_max_edge_index{-1};
    for (const auto &vertex : *this) {
      auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
      for (auto it2 = out_nbr_begin; it2 != out_nbr_end; ++it2) {
        est_max_edge_index = std::max(est_max_edge_index,
                                      (int64_t)it2->second.prop().get_index());
      }
    }
    return est_max_edge_index;
  }();
  VERIFY(next_edge_index >= 1 + est_max_edge_index);

  for (const auto &vertex : *this) {
    if (count_in_neighbors(vertex) == 1 and count_out_neighbors(vertex) == 1) {
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
      VERIFY_MSG(in_nbr_begin->first == vertex,
                 "No 1in-1out vertices are allowed except loops")
      VERIFY_MSG(node_prop(vertex).frozen, "An isolated loop must be frozen");
    }
    auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
    auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
    for (auto in_it = in_nbr_begin; in_it != in_nbr_end; ++in_it) {
      for (auto out_it = out_nbr_begin; out_it != out_nbr_end; ++out_it) {
        in_it->second.prop().assert_incidence(out_it->second.prop(),
                                              node_prop(vertex).len);
      }
    }
  }
}

void MultiplexDBG::move_edge(const RRVertexType &s1, NeighborsIterator e1_it,
                             const RRVertexType &s2, const RRVertexType &e2) {
  // this method by itself does not update read paths
  add_edge_with_prop(s2, e2, std::move(e1_it->second.prop()));
  ConstIterator s1_it = find(s1);
  remove_edge(s1_it, e1_it);
}

void MultiplexDBG::merge_edges(const RRVertexType &s1, NeighborsIterator e1_it,
                               const RRVertexType &s2, NeighborsIterator e2_it,
                               const uint64_t overlap_len) {
  VERIFY_MSG(not node_prop(s2).frozen,
             "Cannot merge edges via a frozen vertex");
  RREdgeProperty &e1_prop = e1_it->second.prop();
  RREdgeProperty &e2_prop = e2_it->second.prop();
  rr_paths->merge(e1_prop.get_index(), e2_prop.get_index());
  e1_prop.merge(std::move(e2_prop), overlap_len);
  move_edge(s1, e1_it, s1, e2_it->first);
  remove_edge(find(s2), e2_it);
  if (e1_it->first != s2) {
    remove_nodes(s2);
  }
  remove_nodes(e1_it->first);
}

EdgeIndexType MultiplexDBG::add_connecting_edge(NeighborsIterator e1_it,
                                                const RRVertexType &s2,
                                                NeighborsIterator e2_it) {
  VERIFY_MSG(e1_it->first != s2, "Can only add edge b/w disconnected edges");
  const uint64_t vertex_len = node_prop(s2).len;
  RREdgeProperty &e1_prop = e1_it->second.prop();
  RREdgeProperty &e2_prop = e2_it->second.prop();
  const EdgeIndexType new_index = next_edge_index;
  ++next_edge_index;
  RREdgeProperty e_new_prop = add(e1_prop, e2_prop, vertex_len, new_index);
  rr_paths->add(e1_prop.get_index(), e2_prop.get_index(),
                e_new_prop.get_index());
  add_edge_with_prop(e1_it->first, s2, std::move(e_new_prop));
  return new_index;
}

RRVertexType MultiplexDBG::get_new_vertex(const uint64_t len) {
  RRVertexType new_vertex{next_vert_index};
  ++next_vert_index;
  RRVertexProperty property{len, false};
  add_node_with_prop(new_vertex, property);
  return new_vertex;
}

MultiplexDBG::MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges,
                           const uint64_t start_k, RRPaths *const rr_paths)
    : rr_paths{rr_paths} {
  for (const SuccinctEdgeInfo &edge : edges) {
    next_vert_index = std::max(next_vert_index, 1 + edge.start_ind);
    next_vert_index = std::max(next_vert_index, 1 + edge.end_ind);
    add_node_with_prop(edge.start_ind, edge.start_prop);
    add_node_with_prop(edge.end_ind, edge.end_prop);
    RREdgeProperty edge_property{next_edge_index, edge.seq, edge.unique};
    add_edge_with_prop(edge.start_ind, edge.end_ind, std::move(edge_property));
    ++next_edge_index;
  }

  freeze_isolated_loops();
  assert_validity();
}

MultiplexDBG::MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *const rr_paths,
                           const uint64_t start_k,
                           UniqueClassificator &classificator, bool debug,
                           const std::experimental::filesystem::path &dir,
                           logging::Logger &logger)
    : rr_paths{rr_paths} {
  const std::unordered_map<std::string, uint64_t> vert2ind = [&dbg, this]() {
    std::unordered_map<std::string, uint64_t> vert2ind;
    for (const Vertex &vertex : dbg.vertices()) {
      const std::string &id = vertex.getId();
      vert2ind.emplace(id, next_vert_index);
      ++next_vert_index;
    }
    return vert2ind;
  }();

  for (auto it = dbg.edges().begin(); it != dbg.edges().end(); ++it) {
    const Edge &edge = *it;
    const RRVertexType start_ind = vert2ind.at(edge.start()->getId());
    const RRVertexType end_ind = vert2ind.at(edge.end()->getId());
    const RRVertexProperty vertex_prop{start_k, false};
    add_node_with_prop(start_ind, vertex_prop);
    add_node_with_prop(end_ind, vertex_prop);

    std::list<char> seq = [&edge]() {
      std::string seq_str = edge.suffix(0).str();
      std::list<char> seq;
      std::move(seq_str.begin(), seq_str.end(), std::back_inserter(seq));
      return seq;
    }();

    RREdgeProperty edge_property{next_edge_index, std::move(seq),
                                 classificator.isUnique(edge)};
    add_edge_with_prop(start_ind, end_ind, std::move(edge_property));
    ++next_edge_index;
  }
  freeze_isolated_loops();
  assert_validity();
}

void MultiplexDBG::serialize_to_dot(
    const std::experimental::filesystem::path &path) const {
  graph_lite::Serializer serializer(*this);
  std::ofstream dot_os(path);
  serializer.serialize_to_dot(dot_os);
}

[[nodiscard]] bool MultiplexDBG::is_frozen() const {
  return std::all_of(begin(), end(), [this](const RRVertexType &v) {
    return node_prop(v).frozen;
  });
}
