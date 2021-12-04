//
// Created by Andrey Bzikadze on 11/25/21.
//

#include "mdbg.hpp"

using namespace repeat_resolution;

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
                                      (int64_t)it2->second.prop().GetIndex());
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
      VERIFY_MSG(node_prop(vertex).IsFrozen(),
                 "An isolated loop must be frozen");
    }
    /*
    auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
    auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
    for (auto in_it = in_nbr_begin; in_it != in_nbr_end; ++in_it) {
      for (auto out_it = out_nbr_begin; out_it != out_nbr_end; ++out_it) {
        in_it->second.prop().assert_incidence(out_it->second.prop(),
                                              node_prop(vertex).size());
      }
    }
     */
  }

  for (const auto &vertex : *this) {
    auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
    const RRVertexProperty &vertex_prop = node_prop(vertex);
    for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
      const RREdgeProperty &edge_prop = it->second.prop();
      const int64_t inner_edge_size = edge_prop.size();
      if (inner_edge_size < 0) {
        const RRVertexProperty &neighbor_prop = node_prop(it->first);
        auto lhs_it = vertex_prop.GetSeq().cend();
        std::advance(lhs_it, inner_edge_size);
        for (auto rhs_it = neighbor_prop.GetSeq().cbegin();
             lhs_it != vertex_prop.GetSeq().cend(); ++lhs_it, ++rhs_it) {
          VERIFY(*lhs_it == *rhs_it);
        }
      }
    }
  }
}

void MultiplexDBG::_spread_frost() {
  std::unordered_set<RRVertexType> prev_frozen, new_frozen;
  for (const RRVertexType &vertex : *this) {
    const RRVertexProperty vertex_prop = node_prop(vertex);
    if (vertex_prop.IsFrozen()) {
      prev_frozen.insert(vertex);
    }
  }

  auto upd_new_frozen = [this, &new_frozen](const RRVertexProperty &vertex_prop,
                                            NeighborsIterator begin,
                                            NeighborsIterator end) {
    for (auto it = begin; it != end; ++it) {
      const RREdgeProperty &edge_prop = it->second.prop();
      const RRVertexType &neighbor = it->first;
      const RRVertexProperty &neighbor_prop = node_prop(neighbor);
      if (not neighbor_prop.IsFrozen() and
          edge_prop.size() == 1 + neighbor_prop.size()) {
        freeze_vertex(neighbor);
        new_frozen.insert(neighbor);
      }
    }
  };

  while (not prev_frozen.empty()) {
    for (const RRVertexType &vertex : prev_frozen) {
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
      auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
      upd_new_frozen(node_prop(vertex), in_nbr_begin, in_nbr_end);
      upd_new_frozen(node_prop(vertex), out_nbr_begin, out_nbr_end);
    }
    prev_frozen = std::move(new_frozen);
  }
}

void MultiplexDBG::freeze_unpaired_vertices() {
  for (const RRVertexType &vertex : *this) {
    RRVertexProperty &vertex_prop = node_prop(vertex);
    if (vertex_prop.IsFrozen()) {
      continue;
    }

    auto [in_edges, out_edges] = get_neighbor_edges_indexes(vertex);
    if (in_edges.size() == 1 and out_edges.size() == 1) {
      // must be a self-loop
      VERIFY(in_edges == out_edges);
      freeze_vertex(vertex);
    } else if (in_edges.size() >= 2 and out_edges.size() >= 2) {
      auto [ac_s2e, ac_e2s] = get_edgepairs_vertex(vertex);
      for (const EdgeIndexType &edge : in_edges) {
        if (ac_s2e.find(edge) == ac_s2e.end()) {
          freeze_vertex(vertex);
          break;
        }
      }
      for (const EdgeIndexType &edge : out_edges) {
        if (ac_e2s.find(edge) == ac_e2s.end()) {
          freeze_vertex(vertex);
          break;
        }
      }
    }

    _spread_frost();
  }
}

void MultiplexDBG::MoveEdge(const RRVertexType &s1, NeighborsIterator e1_it,
                            const RRVertexType &s2, const RRVertexType &e2) {
  // this method by itself does not update read paths
  add_edge_with_prop(s2, e2, std::move(e1_it->second.prop()));
  ConstIterator s1_it = find(s1);
  remove_edge(s1_it, e1_it);
}

void MultiplexDBG::MergeEdges(const RRVertexType &s1, NeighborsIterator e1_it,
                              NeighborsIterator e2_it) {
  const RRVertexType &s2 = e1_it->first;
  VERIFY_MSG(not node_prop(s2).IsFrozen(),
             "Cannot merge edges via a frozen vertex");
  RREdgeProperty &e1_prop = e1_it->second.prop();
  RREdgeProperty &e2_prop = e2_it->second.prop();
  rr_paths->merge(e1_prop.GetIndex(), e2_prop.GetIndex());
  e1_prop.Merge(std::move(node_prop(s1)), std::move(e2_prop));
  MoveEdge(s1, e1_it, s1, e2_it->first);
  remove_edge(find(s2), e2_it);
  remove_nodes(s2);
}

EdgeIndexType MultiplexDBG::AddConnectingEdge(NeighborsIterator eleft_it,
                                              const RRVertexType &vright,
                                              NeighborsIterator eright_it) {
  const RRVertexType &vleft = eleft_it->first;
  VERIFY_MSG(vleft != vright, "Can only add edge b/w disconnected edges");
  const RRVertexProperty &vleft_prop = node_prop(vleft);
  const RRVertexProperty &vright_prop = node_prop(vright);
  const RREdgeProperty &eleft_prop = eleft_it->second.prop();
  const RREdgeProperty &eright_prop = eright_it->second.prop();
  const EdgeIndexType new_index = next_edge_index;
  ++next_edge_index;

  RREdgeProperty e_new_prop = Add(vleft_prop, vright_prop, new_index);
  rr_paths->add(eleft_prop.GetIndex(), eright_prop.GetIndex(),
                e_new_prop.GetIndex());
  add_edge_with_prop(vleft, vright, std::move(e_new_prop));
  return new_index;
}

RRVertexType MultiplexDBG::GetNewVertex(std::list<char> seq) {
  RRVertexType new_vertex{next_vert_index};
  ++next_vert_index;
  RRVertexProperty property(std::move(seq), false);
  add_node_with_prop(new_vertex, std::move(property));
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
    RREdgeProperty edge_property{next_edge_index, edge.seq, edge.infix_size,
                                 edge.unique};
    add_edge_with_prop(edge.start_ind, edge.end_ind, std::move(edge_property));
    ++next_edge_index;
  }

  freeze_unpaired_vertices();
  assert_validity();
}

/*
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

  freeze_unpaired_vertices();
  assert_validity();
}
 */

void MultiplexDBG::serialize_to_dot(
    const std::experimental::filesystem::path &path) const {
  graph_lite::Serializer serializer(*this);
  std::ofstream dot_os(path);
  serializer.serialize_to_dot(dot_os);
}

[[nodiscard]] bool MultiplexDBG::IsFrozen() const {
  return std::all_of(begin(), end(), [this](const RRVertexType &v) {
    return node_prop(v).IsFrozen();
  });
}

std::vector<EdgeIndexType>
MultiplexDBG::get_in_edges_indexes(const RRVertexType &vertex) const {
  std::vector<EdgeIndexType> indexes;
  auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
  for (auto it = in_nbr_begin; it != in_nbr_end; ++it) {
    indexes.push_back(it->second.prop().GetIndex());
  }
  return indexes;
}

bool MultiplexDBG::is_vertex_complex(const RRVertexType &vertex) const {
  const int indegree = count_in_neighbors(vertex);
  const int outdegree = count_out_neighbors(vertex);
  return indegree >= 2 and outdegree >= 2;
}

bool MultiplexDBG::is_vertex_simple(const RRVertexType &vertex) const {
  return not is_vertex_complex(vertex);
}

size_t MultiplexDBG::FullEdgeSize(ConstIterator st_v_it,
                                  NeighborsConstIterator e_it) const {
  const RRVertexType &st_v = *st_v_it;
  const RRVertexType &en_v = e_it->first;
  const RRVertexProperty &st_v_prop = node_prop(st_v_it);
  const RRVertexProperty &en_v_prop = node_prop(en_v);
  const RREdgeProperty &edge_prop = e_it->second.prop();
  int64_t inner_edge_size = edge_prop.size();
  if (inner_edge_size < 0) {
    VERIFY(st_v_prop.size() >= -inner_edge_size);
    VERIFY(en_v_prop.size() >= -inner_edge_size);
  }
  return st_v_prop.size() + inner_edge_size + en_v_prop.size();
}

std::list<char> MultiplexDBG::ExtractEdgePostStartPrefix(ConstIterator st_v_it,
                                                         NeighborsIterator e_it,
                                                         uint64_t len) {
  const RRVertexType &st_v = *st_v_it;
  const RRVertexType &en_v = e_it->first;
  const RRVertexProperty &st_v_prop = node_prop(st_v);
  const RRVertexProperty &en_v_prop = node_prop(en_v);
  RREdgeProperty &edge_prop = e_it->second.prop();
  VERIFY(len + st_v_prop.size() <= FullEdgeSize(st_v_it, e_it));

  uint64_t inner_part_len =
      std::min(len, (uint64_t)std::max(0L, edge_prop.size()));
  std::list<char> prefix = edge_prop.ExtractSeqPrefix(inner_part_len);

  uint64_t en_v_part_len = len - inner_part_len;
  VERIFY(en_v_part_len <= en_v_prop.size());
  prefix.splice(prefix.end(), en_v_prop.GetSeqPrefix(en_v_part_len,
                                                     -edge_prop.size()));
  if (en_v_part_len) {
    edge_prop.ShortenWithEmptySeq(en_v_part_len);
  }
  return prefix;
}

std::list<char> MultiplexDBG::ExtractEdgePreEndSuffix(ConstIterator en_v_it,
                                                      NeighborsIterator e_it,
                                                      uint64_t len) {
  // Edge is reversed
  const RRVertexType &st_v = e_it->first;
  const RRVertexType &en_v = *en_v_it;
  const RRVertexProperty &st_v_prop = node_prop(st_v);
  const RRVertexProperty &en_v_prop = node_prop(en_v);
  RREdgeProperty &edge_prop = e_it->second.prop();
  VERIFY(len + en_v_prop.size() <= FullEdgeSize(find(st_v), e_it));

  uint64_t inner_part_len =
      std::min(len, (uint64_t)std::max(0L, edge_prop.size()));
  uint64_t st_v_part_len = len - inner_part_len;
  VERIFY(st_v_part_len <= st_v_prop.size());
  std::list<char> suffix = st_v_prop.GetSeqSuffix(st_v_part_len,
                                                  -edge_prop.size());
  suffix.splice(suffix.end(), edge_prop.ExtractSeqSuffix(inner_part_len));
  if (st_v_part_len) {
    edge_prop.ShortenWithEmptySeq(st_v_part_len);
  }
  return suffix;
}

void MultiplexDBG::IncreaseVertex(const RRVertexType &vertex, uint64_t len) {
  const int indegree = count_in_neighbors(vertex);
  const int outdegree = count_out_neighbors(vertex);
  VERIFY((indegree == 1) != (outdegree == 1));
  if (indegree == 1) {
    NeighborsIterator edge_it = in_neighbors(vertex).first;
    std::list<char> new_seq =
        ExtractEdgePreEndSuffix(find(vertex), edge_it, len);
    node_prop(vertex).IncLeft(new_seq);
  } else {
    VERIFY(outdegree == 1);
    NeighborsIterator edge_it = out_neighbors(vertex).first;
    std::list<char> new_seq =
        ExtractEdgePostStartPrefix(find(vertex), edge_it, len);
    node_prop(vertex).IncRight(new_seq);
  }
}

std::vector<EdgeIndexType>
MultiplexDBG::get_out_edges_indexes(const RRVertexType &vertex) const {
  std::vector<EdgeIndexType> indexes;
  auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
  for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
    indexes.push_back(it->second.prop().GetIndex());
  }
  return indexes;
}

std::pair<std::vector<EdgeIndexType>, std::vector<EdgeIndexType>>
MultiplexDBG::get_neighbor_edges_indexes(const RRVertexType &vertex) const {
  return {get_in_edges_indexes(vertex), get_out_edges_indexes(vertex)};
}

std::pair<MultiplexDBG::EdgeNeighborMap, MultiplexDBG::EdgeNeighborMap>
MultiplexDBG::get_edgepairs_vertex(const RRVertexType &vertex) const {
  auto get_init_transitions =
      [this](const std::vector<EdgeIndexType> &in_edges,
             const std::vector<EdgeIndexType> &out_edges) {
        EdgeNeighborMap ac_s2e, ac_e2s;
        for (const EdgeIndexType &in_ind : in_edges) {
          for (const EdgeIndexType &out_ind : out_edges) {
            if (rr_paths->contains_pair(in_ind, out_ind)) {
              ac_s2e[in_ind].emplace(out_ind);
              ac_e2s[out_ind].emplace(in_ind);
            }
          }
        }
        return std::make_pair(ac_s2e, ac_e2s);
      };

  auto extend_transitions_single_loop =
      [this, &vertex](const std::vector<EdgeIndexType> &in_edges,
                      const std::vector<EdgeIndexType> &out_edges,
                      EdgeNeighborMap &ac_s2e, EdgeNeighborMap &ac_e2s) {
        std::vector<EdgeIndexType> loops;
        for (const EdgeIndexType &index : in_edges) {
          if (std::find(out_edges.begin(), out_edges.end(), index) !=
              out_edges.end()) {
            loops.push_back(index);
          }
        }

        if (loops.size() == 1) {
          const EdgeIndexType loop = loops.front();
          const RREdgeProperty &loop_prop =
              find_out_edge_constiterator(vertex, loop)->second.prop();
          if (loop_prop.IsUnique()) {
            if (in_edges.size() == 2) {
              const size_t loop_index = in_edges.back() == loop;
              const size_t nonloop = in_edges[loop_index ^ 1];
              ac_s2e[nonloop].emplace(loop);
              ac_e2s[loop].emplace(nonloop);
            }
            if (out_edges.size() == 2) {
              const size_t loop_index = out_edges.back() == loop;
              const size_t nonloop = out_edges[loop_index ^ 1];
              ac_s2e[loop].emplace(nonloop);
              ac_e2s[nonloop].emplace(loop);
            }
          }
        }
      };

  auto extend_transitions_all_unique =
      [this, &vertex](const std::vector<EdgeIndexType> &in_edges,
                      const std::vector<EdgeIndexType> &out_edges,
                      EdgeNeighborMap &ac_s2e, EdgeNeighborMap &ac_e2s) {
        std::vector<EdgeIndexType> unpaired_in, unpaired_out;
        for (const EdgeIndexType &index : out_edges) {
          if (ac_e2s.find(index) == ac_e2s.end()) {
            unpaired_out.push_back(index);
          }
        }
        for (const EdgeIndexType &index : in_edges) {
          if (ac_e2s.find(index) == ac_e2s.end()) {
            unpaired_in.push_back(index);
          }
        }

        bool all_in_unique =
            std::all_of(in_edges.begin(), in_edges.end(),
                        [this, &vertex](const EdgeIndexType &edge_index) {
                          return find_in_edge_constiterator(vertex, edge_index)
                              ->second.prop()
                              .IsUnique();
                        });
        bool all_out_unique =
            std::all_of(out_edges.begin(), out_edges.end(),
                        [this, &vertex](const EdgeIndexType &edge_index) {
                          return find_out_edge_constiterator(vertex, edge_index)
                              ->second.prop()
                              .IsUnique();
                        });
        if (unpaired_in.size() == 1 and unpaired_out.size() == 1 and
            (all_in_unique or all_out_unique)) {
          ac_s2e.emplace(unpaired_in.front(), unpaired_out.front());
          ac_e2s.emplace(unpaired_out.front(), unpaired_in.front());
        }
      };

  const auto [in_edges, out_edges] = get_neighbor_edges_indexes(vertex);
  auto [ac_s2e, ac_e2s] = get_init_transitions(in_edges, out_edges);
  extend_transitions_single_loop(in_edges, out_edges, ac_s2e, ac_e2s);
  extend_transitions_all_unique(in_edges, out_edges, ac_s2e, ac_e2s);

  return std::make_pair(ac_s2e, ac_e2s);
}

MultiplexDBG::NeighborsIterator
MultiplexDBG::find_in_edge_iterator(const RRVertexType &v,
                                    const EdgeIndexType &edge) {
  auto [it, end] = in_neighbors(v);
  while (it != end and it->second.prop().GetIndex() != edge) {
    ++it;
  }
  return it;
};

MultiplexDBG::NeighborsConstIterator
MultiplexDBG::find_in_edge_constiterator(const RRVertexType &v,
                                         const EdgeIndexType &edge) const {
  auto [it, end] = in_neighbors(v);
  while (it != end and it->second.prop().GetIndex() != edge) {
    ++it;
  }
  return it;
};

MultiplexDBG::NeighborsIterator
MultiplexDBG::find_out_edge_iterator(const RRVertexType &v,
                                     const EdgeIndexType &edge) {
  auto [it, end] = out_neighbors(v);
  while (it != end and it->second.prop().GetIndex() != edge) {
    ++it;
  }
  return it;
};

MultiplexDBG::NeighborsConstIterator
MultiplexDBG::find_out_edge_constiterator(const RRVertexType &v,
                                          const EdgeIndexType &edge) const {
  auto [it, end] = out_neighbors(v);
  while (it != end and it->second.prop().GetIndex() != edge) {
    ++it;
  }
  return it;
};
