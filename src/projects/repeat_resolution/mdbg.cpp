//
// Created by Andrey Bzikadze on 11/25/21.
//

#include "mdbg.hpp"

using namespace repeat_resolution;

void MultiplexDBG::AssertValidity() const {
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
        est_max_edge_index =
            std::max(est_max_edge_index, (int64_t)it2->second.prop().Index());
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
  }

  for (const auto &vertex : *this) {
    auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
    const RRVertexProperty &vertex_prop = node_prop(vertex);
    for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
      const RREdgeProperty &edge_prop = it->second.prop();
      const int64_t inner_edge_size = edge_prop.size();
      if (inner_edge_size < 0) {
        const RRVertexProperty &neighbor_prop = node_prop(it->first);
        auto lhs_it = vertex_prop.Seq().cend();
        std::advance(lhs_it, inner_edge_size);
        for (auto rhs_it = neighbor_prop.Seq().cbegin();
             lhs_it != vertex_prop.Seq().cend(); ++lhs_it, ++rhs_it) {
          VERIFY(*lhs_it == *rhs_it);
        }
      }
    }
  }

  {
    for (const RRVertexType &vertex : *this) {
      const RRVertexProperty &vertex_prop = node_prop(vertex);
      if (not vertex_prop.IsFrozen()) {
        VERIFY(vertex_prop.size() == n_iter + start_k);
      }
    }
  }

  if (contains_rc) {
    std::unordered_set<std::string> seq_vertex;
    for (const RRVertexType &vertex : *this) {
      const RRVertexProperty &vertex_prop = node_prop(vertex);
      const std::string seq = vertex_prop.Seq().ToString();
      seq_vertex.emplace(seq);
    }

    for (const RRVertexType &vertex : *this) {
      const RRVertexProperty &vertex_prop = node_prop(vertex);
      VERIFY(seq_vertex.count(vertex_prop.Seq().GetRC().ToString()) == 1);
    }

    std::unordered_set<std::string> seq_edge;
    for (const RRVertexType &vertex : *this) {
      auto [begin, end] = out_neighbors(vertex);
      for (auto it = begin; it != end; ++it) {
        const RREdgeProperty &edge_prop = it->second.prop();
        if (edge_prop.size() > 0) {
          seq_edge.emplace(
              GetEdgeSequence(find(vertex), it, false, false).ToString());
        }
      }
    }

    for (const RRVertexType &vertex : *this) {
      auto [begin, end] = out_neighbors(vertex);
      for (auto it = begin; it != end; ++it) {
        const RREdgeProperty &edge_prop = it->second.prop();
        if (edge_prop.size() > 0) {
          const std::string seq =
              GetEdgeSequence(find(vertex), it, false, false).ToString();
          VERIFY(seq_edge.count(seq) == 1);
        }
      }
    }
  }
}

void MultiplexDBG::SpreadFrost() {
  std::unordered_set<RRVertexType> prev_frozen, new_frozen;
  for (const RRVertexType &vertex : *this) {
    const RRVertexProperty &vertex_prop = node_prop(vertex);
    if (vertex_prop.IsFrozen()) {
      prev_frozen.insert(vertex);
    }
  }

  auto upd_new_frozen = [this, &new_frozen](const RRVertexType &vertex,
                                            const RRVertexProperty &vertex_prop,
                                            NeighborsIterator begin,
                                            NeighborsIterator end) {
    for (auto it = begin; it != end; ++it) {
      const RREdgeProperty &edge_prop = it->second.prop();
      const RRVertexType &neighbor = it->first;
      const RRVertexProperty &neighbor_prop = node_prop(neighbor);
      if (not neighbor_prop.IsFrozen() and
          FullEdgeSize(find(vertex), it) == 1 + neighbor_prop.size()) {
        FreezeVertex(neighbor);
        new_frozen.insert(neighbor);
      }
    }
  };

  while (not prev_frozen.empty()) {
    for (const RRVertexType &vertex : prev_frozen) {
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
      auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
      upd_new_frozen(vertex, node_prop(vertex), in_nbr_begin, in_nbr_end);
      upd_new_frozen(vertex, node_prop(vertex), out_nbr_begin, out_nbr_end);
    }
    prev_frozen = std::move(new_frozen);
  }
}

void MultiplexDBG::FreezeUnpairedVertices() {
  for (const RRVertexType &vertex : *this) {
    RRVertexProperty &vertex_prop = node_prop(vertex);
    if (vertex_prop.IsFrozen()) {
      continue;
    }

    auto [in_edges, out_edges] = GetNeighborEdgesIndexes(vertex);
    if (in_edges.size() == 1 and out_edges.size() == 1) {
      // must be a self-loop
      VERIFY(in_edges == out_edges);
      FreezeVertex(vertex);
    } else if (in_edges.size() >= 2 and out_edges.size() >= 2) {
      auto [ac_s2e, ac_e2s] = GetEdgepairsVertex(vertex);
      for (const RREdgeIndexType &edge : in_edges) {
        if (ac_s2e.find(edge) == ac_s2e.end()) {
          FreezeVertex(vertex);
          break;
        }
      }
      for (const RREdgeIndexType &edge : out_edges) {
        if (ac_e2s.find(edge) == ac_e2s.end()) {
          FreezeVertex(vertex);
          break;
        }
      }
    }

    SpreadFrost();
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
  const RREdgeIndexType e2_index = e2_prop.Index();
  rr_paths->merge(e1_prop.Index(), e2_prop.Index());
  const RRVertexProperty &v1 = node_prop(s1);
  const RRVertexProperty &v3 = node_prop(e2_it->first);
  e1_prop.Merge(std::move(node_prop(s2)), std::move(e2_prop));
  MoveEdge(s1, e1_it, s1, e2_it->first);
  remove_edge(find(s2), FindOutEdgeIterator(s2, e2_index));
  remove_nodes(s2);
}

RREdgeIndexType MultiplexDBG::AddConnectingEdge(NeighborsIterator eleft_it,
                                              const RRVertexType &vright,
                                              NeighborsIterator eright_it) {
  const RRVertexType &vleft = eleft_it->first;
  VERIFY_MSG(vleft != vright, "Can only add edge b/w disconnected edges");
  const RRVertexProperty &vleft_prop = node_prop(vleft);
  const RRVertexProperty &vright_prop = node_prop(vright);
  VERIFY(vleft_prop.size() == vright_prop.size());
  {
    auto lit = vleft_prop.Seq().cbegin();
    ++lit;
    auto rit = vright_prop.Seq().cbegin();
    // while (lit != vleft_prop.Seq().cend()) {
    //   VERIFY(*lit == *rit);
    //   ++lit, ++rit;
    // }
  }

  const RREdgeProperty &eleft_prop = eleft_it->second.prop();
  const RREdgeProperty &eright_prop = eright_it->second.prop();
  const RREdgeIndexType new_index = next_edge_index;
  ++next_edge_index;

  RREdgeProperty e_new_prop = Add(vleft_prop, vright_prop, new_index);
  rr_paths->add(eleft_prop.Index(), eright_prop.Index(), e_new_prop.Index());
  add_edge_with_prop(vleft, vright, std::move(e_new_prop));
  return new_index;
}

RRVertexType MultiplexDBG::GetNewVertex(MDBGSeq seq) {
  RRVertexType new_vertex{next_vert_index};
  ++next_vert_index;
  RRVertexProperty property(std::move(seq), false);
  add_node_with_prop(new_vertex, std::move(property));
  return new_vertex;
}

MultiplexDBG::MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges,
                           const uint64_t start_k, RRPaths *const rr_paths)
    : rr_paths{rr_paths}, start_k{start_k}, contains_rc{false} {
  for (const SuccinctEdgeInfo &edge : edges) {
    next_vert_index = std::max(next_vert_index, 1 + edge.start_ind);
    next_vert_index = std::max(next_vert_index, 1 + edge.end_ind);
    add_node_with_prop(
        edge.start_ind,
        RRVertexProperty(edge.start_prop.Seq(), edge.start_prop.IsFrozen()));
    add_node_with_prop(
        edge.end_ind,
        RRVertexProperty(edge.end_prop.Seq(), edge.end_prop.IsFrozen()));
    RREdgeProperty edge_property(next_edge_index, edge.seq, edge.infix_size,
                                 edge.unique);
    add_edge_with_prop(edge.start_ind, edge.end_ind, std::move(edge_property));
    ++next_edge_index;
  }

  FreezeUnpairedVertices();
  AssertValidity();
}

MultiplexDBG::MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *const rr_paths,
                           const uint64_t start_k,
                           UniqueClassificator &classificator, bool debug,
                           const std::experimental::filesystem::path &dir,
                           logging::Logger &logger)
    : rr_paths{rr_paths}, start_k{start_k} {
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

    const Sequence &seq = edge.suffix(0);

    Sequence prefix = seq.Prefix(start_k);
    int64_t infix_size = ((int64_t)seq.size()) - 2 * start_k;
    Sequence infix =
        infix_size > 0 ? seq.Subseq(start_k, start_k + infix_size) : Sequence();
    Sequence suffix = seq.Suffix(start_k);

    RRVertexProperty st_v_prop(MDBGSeq(prefix), false);
    add_node_with_prop(start_ind, std::move(st_v_prop));

    RRVertexProperty en_v_prop(MDBGSeq(suffix), false);
    add_node_with_prop(end_ind, std::move(en_v_prop));

    RREdgeProperty edge_prop(next_edge_index, MDBGSeq(infix), infix_size,
                             classificator.isUnique(edge));
    add_edge_with_prop(start_ind, end_ind, std::move(edge_prop));
    ++next_edge_index;
  }

  FreezeUnpairedVertices();
  AssertValidity();
}

void MultiplexDBG::SerializeToDot(
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

std::vector<RREdgeIndexType>
MultiplexDBG::GetInEdgesIndexes(const RRVertexType &vertex) const {
  std::vector<RREdgeIndexType> indexes;
  auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
  for (auto it = in_nbr_begin; it != in_nbr_end; ++it) {
    indexes.push_back(it->second.prop().Index());
  }
  return indexes;
}

bool MultiplexDBG::IsVertexComplex(const RRVertexType &vertex) const {
  const int indegree = count_in_neighbors(vertex);
  const int outdegree = count_out_neighbors(vertex);
  return indegree >= 2 and outdegree >= 2;
}

bool MultiplexDBG::IsVertexSimple(const RRVertexType &vertex) const {
  return not IsVertexComplex(vertex);
}

bool MultiplexDBG::IsVertexCanonical(const RRVertexType &vertex) const {
  const RRVertexProperty &vertex_prop = node_prop(vertex);
  return vertex_prop.IsCanonical();
}

bool MultiplexDBG::IsEdgeCanonical(ConstIterator vertex,
                                   NeighborsConstIterator e_it) const {
  const RREdgeProperty &edge_prop = e_it->second.prop();
  const MDBGSeq edge_seq = GetEdgeSequence(vertex, e_it, false, false);
  return edge_seq.IsCanonical();
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

MDBGSeq MultiplexDBG::ExtractEdgePostStartPrefix(ConstIterator st_v_it,
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
  MDBGSeq prefix = edge_prop.ExtractSeqPrefix(inner_part_len);

  uint64_t en_v_part_len = len - inner_part_len;
  VERIFY(en_v_part_len <= en_v_prop.size());
  prefix.Append(en_v_prop.GetSeqPrefix(en_v_part_len, -edge_prop.size()));
  if (en_v_part_len) {
    edge_prop.ShortenWithEmptySeq(en_v_part_len);
  }
  return prefix;
}

MDBGSeq MultiplexDBG::ExtractEdgePreEndSuffix(ConstIterator en_v_it,
                                              NeighborsIterator e_it,
                                              uint64_t len) {
  // Edge is reversed
  const RRVertexType &st_v = e_it->first;
  const RRVertexType &en_v = *en_v_it;
  const RRVertexProperty &st_v_prop = node_prop(st_v);
  const RRVertexProperty &en_v_prop = node_prop(en_v);
  RREdgeProperty &edge_prop = e_it->second.prop();
  size_t full_edge_size = FullEdgeSize(find(en_v), e_it);
  VERIFY(len + en_v_prop.size() <= full_edge_size);

  uint64_t inner_part_len =
      std::min(len, (uint64_t)std::max(0L, edge_prop.size()));
  uint64_t st_v_part_len = len - inner_part_len;
  VERIFY(st_v_part_len <= st_v_prop.size());
  MDBGSeq suffix = st_v_prop.GetSeqSuffix(st_v_part_len, -edge_prop.size());
  suffix.Append(edge_prop.ExtractSeqSuffix(inner_part_len));
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
    MDBGSeq new_seq = ExtractEdgePreEndSuffix(find(vertex), edge_it, len);
    node_prop(vertex).IncLeft(std::move(new_seq));
  } else {
    VERIFY(outdegree == 1);
    NeighborsIterator edge_it = out_neighbors(vertex).first;
    MDBGSeq new_seq = ExtractEdgePostStartPrefix(find(vertex), edge_it, len);
    node_prop(vertex).IncRight(std::move(new_seq));
  }
}

std::vector<RREdgeIndexType>
MultiplexDBG::GetOutEdgesIndexes(const RRVertexType &vertex) const {
  std::vector<RREdgeIndexType> indexes;
  auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
  for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
    indexes.push_back(it->second.prop().Index());
  }
  return indexes;
}

std::pair<std::vector<RREdgeIndexType>, std::vector<RREdgeIndexType>>
MultiplexDBG::GetNeighborEdgesIndexes(const RRVertexType &vertex) const {
  return {GetInEdgesIndexes(vertex), GetOutEdgesIndexes(vertex)};
}

std::pair<MultiplexDBG::EdgeNeighborMap, MultiplexDBG::EdgeNeighborMap>
MultiplexDBG::GetEdgepairsVertex(const RRVertexType &vertex) const {
  auto get_init_transitions =
      [this](const std::vector<RREdgeIndexType> &in_edges,
             const std::vector<RREdgeIndexType> &out_edges) {
        EdgeNeighborMap ac_s2e, ac_e2s;
        for (const RREdgeIndexType &in_ind : in_edges) {
          for (const RREdgeIndexType &out_ind : out_edges) {
            if (rr_paths->contains_pair(in_ind, out_ind)) {
              ac_s2e[in_ind].emplace(out_ind);
              ac_e2s[out_ind].emplace(in_ind);
            }
          }
        }
        return std::make_pair(ac_s2e, ac_e2s);
      };

  auto extend_transitions_single_loop =
      [this, &vertex](const std::vector<RREdgeIndexType> &in_edges,
                      const std::vector<RREdgeIndexType> &out_edges,
                      EdgeNeighborMap &ac_s2e, EdgeNeighborMap &ac_e2s) {
        std::vector<RREdgeIndexType> loops;
        for (const RREdgeIndexType &index : in_edges) {
          if (std::find(out_edges.begin(), out_edges.end(), index) !=
              out_edges.end()) {
            loops.push_back(index);
          }
        }

        if (loops.size() == 1) {
          const RREdgeIndexType loop = loops.front();
          const RREdgeProperty &loop_prop =
              FindOutEdgeConstiterator(vertex, loop)->second.prop();
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
      [this, &vertex](const std::vector<RREdgeIndexType> &in_edges,
                      const std::vector<RREdgeIndexType> &out_edges,
                      EdgeNeighborMap &ac_s2e, EdgeNeighborMap &ac_e2s) {
        std::vector<RREdgeIndexType> unpaired_in, unpaired_out;
        for (const RREdgeIndexType &index : out_edges) {
          if (ac_e2s.find(index) == ac_e2s.end()) {
            unpaired_out.push_back(index);
          }
        }
        for (const RREdgeIndexType &index : in_edges) {
          if (ac_e2s.find(index) == ac_e2s.end()) {
            unpaired_in.push_back(index);
          }
        }

        bool all_in_unique =
            std::all_of(in_edges.begin(), in_edges.end(),
                        [this, &vertex](const RREdgeIndexType &edge_index) {
                          return FindInEdgeConstiterator(vertex, edge_index)
                              ->second.prop()
                              .IsUnique();
                        });
        bool all_out_unique =
            std::all_of(out_edges.begin(), out_edges.end(),
                        [this, &vertex](const RREdgeIndexType &edge_index) {
                          return FindOutEdgeConstiterator(vertex, edge_index)
                              ->second.prop()
                              .IsUnique();
                        });
        if (unpaired_in.size() == 1 and unpaired_out.size() == 1 and
            (all_in_unique or all_out_unique)) {
          ac_s2e.emplace(unpaired_in.front(), unpaired_out.front());
          ac_e2s.emplace(unpaired_out.front(), unpaired_in.front());
        }
      };

  const auto [in_edges, out_edges] = GetNeighborEdgesIndexes(vertex);
  auto [ac_s2e, ac_e2s] = get_init_transitions(in_edges, out_edges);
  extend_transitions_single_loop(in_edges, out_edges, ac_s2e, ac_e2s);
  extend_transitions_all_unique(in_edges, out_edges, ac_s2e, ac_e2s);

  return std::make_pair(ac_s2e, ac_e2s);
}

MultiplexDBG::NeighborsIterator
MultiplexDBG::FindInEdgeIterator(const RRVertexType &v,
                                 const RREdgeIndexType &edge) {
  auto [it, end] = in_neighbors(v);
  while (it != end and it->second.prop().Index() != edge) {
    ++it;
  }
  return it;
};

MultiplexDBG::NeighborsConstIterator
MultiplexDBG::FindInEdgeConstiterator(const RRVertexType &v,
                                      const RREdgeIndexType &edge) const {
  auto [it, end] = in_neighbors(v);
  while (it != end and it->second.prop().Index() != edge) {
    ++it;
  }
  return it;
};

MultiplexDBG::NeighborsIterator
MultiplexDBG::FindOutEdgeIterator(const RRVertexType &v,
                                  const RREdgeIndexType &edge) {
  auto [it, end] = out_neighbors(v);
  while (it != end and it->second.prop().Index() != edge) {
    ++it;
  }
  return it;
};

MultiplexDBG::NeighborsConstIterator
MultiplexDBG::FindOutEdgeConstiterator(const RRVertexType &v,
                                       const RREdgeIndexType &edge) const {
  auto [it, end] = out_neighbors(v);
  while (it != end and it->second.prop().Index() != edge) {
    ++it;
  }
  return it;
}

int64_t MultiplexDBG::GetInnerEdgeSize(ConstIterator vertex,
                                       NeighborsConstIterator e_it) const {
  return e_it->second.prop().size();
}

MDBGSeq MultiplexDBG::GetEdgeSequence(ConstIterator vertex,
                                      NeighborsConstIterator e_it,
                                      bool trim_left, bool trim_right) const {
  const RRVertexProperty &vertex_prop = node_prop(vertex);
  const RRVertexProperty &neighbor_prop = node_prop(e_it->first);
  const RREdgeProperty &edge_prop = e_it->second.prop();

  if (edge_prop.size() < 0) {
    if (trim_left and trim_right) {
      return {};
    }
    if (trim_left) {
      MDBGSeq neighbor_seq = neighbor_prop.Seq();
      neighbor_seq.TrimLeft(-edge_prop.size());
      return neighbor_seq;
    }
    if (trim_right) {
      MDBGSeq vertex_seq = vertex_prop.Seq();
      vertex_seq.TrimRight(-edge_prop.size());
      return vertex_seq;
    }
  }

  MDBGSeq seq = edge_prop.Seq();
  if (not trim_left) {
    MDBGSeq vertex_seq = vertex_prop.Seq();
    seq.Prepend(std::move(vertex_seq));
  }
  if (not trim_right) {
    MDBGSeq neighbor_seq = neighbor_prop.Seq();
    seq.Append(std::move(neighbor_seq));
  }
  return seq;
}

std::vector<Contig>
MultiplexDBG::GetTrimEdges(int64_t min_inner_edge_size) const {
  const std::unordered_map<RRVertexType, bool> trim = [this]() {
    std::unordered_map<std::string, RRVertexType> seq2vertex;
    for (const RRVertexType &vertex : *this) {
      const RRVertexProperty &vertex_prop = node_prop(vertex);
      seq2vertex.emplace(vertex_prop.Seq().ToString(), vertex);
    }

    std::unordered_map<RRVertexType, bool> trim;
    for (const RRVertexType &vertex : *this) {
      const RRVertexProperty &vertex_prop = node_prop(vertex);
      if (IsVertexCanonical(vertex)) {
        const bool trim_vertex = count_out_neighbors(vertex) != 1;
        trim.emplace(vertex, trim_vertex);
        trim.emplace(seq2vertex.at(vertex_prop.Seq().GetRC().ToString()),
                     not trim_vertex);
      }
    }
    return trim;
  }();

  std::vector<Contig> contigs;
  for (const RRVertexType &vertex : *this) {
    auto vertex_it = find(vertex);
    auto [out_begin, out_end] = out_neighbors(vertex);
    for (auto it = out_begin; it != out_end; ++it) {
      // std::cout << vertex << " " << IsVertexCanonical(vertex) << " "
      //           << IsVertexCanonical(it->first) << " "
      //           << it->second.prop().Index() << " "
      //           << GetInnerEdgeSize(vertex_it, it) << " "
      //           << IsEdgeCanonical(vertex_it, it) << " " << trim.at(vertex)
      //           << " " << trim.at(it->first) << "\n";
      if (IsEdgeCanonical(vertex_it, it) and
          GetInnerEdgeSize(vertex_it, it) >= min_inner_edge_size) {
        MDBGSeq edge_list = GetEdgeSequence(vertex_it, it, trim.at(vertex),
                                            not trim.at(it->first));

        std::string edge_str = edge_list.ToString();
        contigs.emplace_back(std::move(edge_str),
                             itos(it->second.prop().Index()));
      }
    }
  }
  return contigs;
}

std::vector<Contig> MultiplexDBG::PrintTrimEdges(
    const std::experimental::filesystem::path &f) const {
  std::ofstream os;
  os.open(f);
  std::vector<Contig> edges = GetTrimEdges();
  for (const Contig &contig : edges) {
    os << ">" << contig.id << "\n" << contig.seq << "\n";
  }
  os.close();
  return edges;
}