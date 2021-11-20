//
// Created by Andrey Bzikadze on 11/10/21.
//

#pragma once

#include "error_correction/multiplicity_estimation.hpp"
#include "multiplex_dbg_topology.hpp"
#include "paths.hpp"

namespace repeat_resolution {

class MultiplexDBG
    : public graph_lite::Graph<
          /*typename NodeType=*/RRVertexType,
          /*typename NodePropType=*/RRVertexProperty,
          /*typename EdgePropType=*/RREdgeProperty,
          /*EdgeDirection direction=*/graph_lite::EdgeDirection::DIRECTED,
          /*MultiEdge multi_edge=*/graph_lite::MultiEdge::ALLOWED,
          /*SelfLoop self_loop=*/graph_lite::SelfLoop::ALLOWED,
          /*Map adj_list_spec=*/graph_lite::Map::UNORDERED_MAP,
          /*Container neighbors_container_spec=*/
          graph_lite::Container::MULTISET> {
  RRPaths *rr_paths;
  uint64_t max_edge_index{0};
  uint64_t max_vert_index{0};
  uint64_t niter{0};

  bool node_has_loop(const node_type &node) const {
    auto [begin, end] = in_neighbors(node);
    for (auto it = begin; it != end; ++it) {
      const node_type &neighbor = it->first;
      if (neighbor == node) {
        return true;
      }
    }
    return false;
  }

  void assert_validity() const {
    uint64_t est_max_vert_index = [this]() {
      uint64_t est_max_vert_index{0};
      for (const auto &vertex : *this) {
        est_max_vert_index = std::max(est_max_vert_index, vertex);
      }
      return est_max_vert_index;
    }();
    VERIFY(max_vert_index >= 1 + est_max_vert_index);

    uint64_t est_max_edge_index = [this]() {
      uint64_t est_max_edge_index{0};
      for (const auto &vertex : *this) {
        auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
        for (auto it2 = out_nbr_begin; it2 != out_nbr_end; ++it2) {
          est_max_edge_index =
              std::max(est_max_edge_index, it2->second.prop().get_index());
        }
      }
      return est_max_edge_index;
    }();
    VERIFY(max_edge_index >= 1 + est_max_edge_index);

    for (const auto &vertex : *this) {
      VERIFY(count_in_neighbors(vertex) != 1 or
             count_out_neighbors(vertex) != 1);
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

  void move_edge(const RRVertexType &s1, NeighborsIterator e1_it,
                 const RRVertexType &s2, const RRVertexType &e2) {
    // this method by itself does not update read paths
    add_edge_with_prop(s2, e2, std::move(e1_it->second.prop()));
    ConstIterator s1_it = find(s1);
    remove_edge(s1_it, e1_it);
  }

  void smart_remove_edge(ConstIterator s_it, NeighborsIterator e_it,
                         bool moving) {
    const RREdgeProperty &edge_prop = e_it->second.prop();
    rr_paths->remove(edge_prop.get_index());

    remove_edge(s_it, e_it);
    if (moving) {
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(e_it->first);
      for (auto in_nbr_it = in_nbr_begin; in_nbr_it != in_nbr_end;
           ++in_nbr_it) {
        move_edge(e_it->first, in_nbr_it, in_nbr_it->first, *s_it);
      }

      auto [out_nbr_begin, out_nbr_end] = out_neighbors(e_it->first);
      for (auto out_nbr_it = out_nbr_begin; out_nbr_it != out_nbr_end;
           ++out_nbr_it) {
        move_edge(e_it->first, out_nbr_it, *s_it, out_nbr_it->first);
      }
    }
  }

  void merge_edges(ConstIterator s1_it, NeighborsIterator e1_it,
                   ConstIterator s2_it, NeighborsIterator e2_it) {
    VERIFY_MSG(e1_it->first == *s2_it, "Can only merge adjacent edges");
    VERIFY_MSG(not node_prop(s2_it).frozen,
               "Cannot merge edges via a frozen vertex");
    const uint64_t overlap_len = node_prop(s2_it).len;
    RREdgeProperty &e1_prop = e1_it->second.prop();
    RREdgeProperty &e2_prop = e2_it->second.prop();
    rr_paths->merge(e1_prop.get_index(), e2_prop.get_index());
    e1_prop.merge(std::move(e2_prop), overlap_len);
    move_edge(*s1_it, e1_it, *s1_it, e2_it->first);
    smart_remove_edge(s2_it, e2_it, false);
  }

  void smart_add_edge(ConstIterator s1_it, NeighborsIterator e1_it,
                      ConstIterator s2_it, NeighborsIterator e2_it) {
    VERIFY_MSG(e1_it->first != *s2_it,
               "Can only add edge b/w disconnected edges");
    const uint64_t overlap_len = node_prop(s2_it).len;
    RREdgeProperty &e1_prop = e1_it->second.prop();
    RREdgeProperty &e2_prop = e2_it->second.prop();
    RREdgeProperty e_new_prop =
        add(e1_prop, e2_prop, overlap_len, max_edge_index);
    ++max_edge_index;
    rr_paths->add(e1_prop.get_index(), e2_prop.get_index(),
                  e_new_prop.get_index());
    add_edge_with_prop(e1_it->first, *s2_it, std::move(e_new_prop));
  }

  [[nodiscard]] bool is_frozen() const {
    return std::all_of(begin(), end(), [this](const RRVertexType &v) {
      return node_prop(v).frozen;
    });
  }

  RRVertexType get_new_vertex(const uint64_t len) {
    RRVertexType new_vertex{max_vert_index};
    ++max_vert_index;
    RRVertexProperty property{len, false};
    add_node_with_prop(new_vertex, property);
    return new_vertex;
  }

  void process_simple_vertex(const RRVertexType &vertex, const int indegree,
                             const int outdegree) {
    VERIFY_MSG(indegree != 1 or outdegree != 1,
               "no vertexes on nonbranching paths allowed");
    RRVertexProperty &v_prop = node_prop(vertex);
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
      auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
      for (auto it = out_nbr_begin; it != out_nbr_end; ++it) {
        RRVertexType new_vertex = get_new_vertex(v_prop.len + 1);
        move_edge(vertex, it, new_vertex, it->first);
      }
      remove_nodes(vertex); // careful: Iterator is invalidated

    } else if (indegree > 1 and outdegree == 0) {
      // "Finishing" vertex
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
      for (auto it = in_nbr_begin; it != in_nbr_end; ++it) {
        RRVertexType new_vertex = get_new_vertex(v_prop.len + 1);
        // need to construct a NeighborIterator pointing to vertex
        auto out_nbr = out_neighbors(it->first).first;
        while (out_nbr->first != vertex) {
          ++out_nbr;
        }
        move_edge(it->first, out_nbr, it->first, new_vertex);
      }
      remove_nodes(vertex); // careful: Iterator is invalidated

    } else if (indegree == 1 and outdegree > 1) {
      // TODO

    } else if (indegree > 1 and outdegree == 1) {
      // TODO
    }
  }

  void process_complex_vertex(const RRVertexType &vertex, const int indegree,
                              const int outdegree) {}

  void process_vertex(const RRVertexType &vertex) {
    const int indegree = count_in_neighbors(vertex);
    const int outdegree = count_out_neighbors(vertex);
    if (indegree >= 2 and outdegree >= 2) {
      process_complex_vertex(vertex, indegree, outdegree);
    } else {
      process_simple_vertex(vertex, indegree, outdegree);
    }
  }

  void finalize_transformation() {
    for (const RRVertexType &v1 : *this) {
      const RRVertexProperty &v1p = node_prop(v1);
      auto [out_it_begin, out_it_end] = out_neighbors(v1);
      for (auto it = out_it_begin; it != out_it_end; ++it) {
        const RRVertexType &v2 = it->first;
        const RRVertexProperty &v2p = node_prop(v2);
        const RREdgeProperty &edge_property = it->second.prop();
        if (v1p.len == edge_property.size() or
            v2p.len == edge_property.size()) {
          VERIFY(v1p.len == v2p.len);
          VERIFY(not v1p.frozen and not v2p.frozen);
          smart_remove_edge(find(v1), it, true);
          if (count_in_neighbors(v2) + count_out_neighbors(v2) == 0) {
            remove_nodes(v2);
          }
        }
      }
    }
  }

public:
  // This constructor is for testing purposes
  MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges,
               const uint64_t start_k, RRPaths *const rr_paths)
      : rr_paths{rr_paths} {
    for (const SuccinctEdgeInfo &edge : edges) {
      max_vert_index = std::max(max_vert_index, 1 + edge.start_ind);
      max_vert_index = std::max(max_vert_index, 1 + edge.end_ind);
      add_node_with_prop(edge.start_ind, edge.start_prop);
      add_node_with_prop(edge.end_ind, edge.end_prop);
      RREdgeProperty edge_property{max_edge_index, edge.seq, edge.unique};
      add_edge_with_prop(edge.start_ind, edge.end_ind,
                         std::move(edge_property));
      ++max_edge_index;
    }
    assert_validity();
  }

  MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *const rr_paths,
               const uint64_t start_k, UniqueClassificator &classificator,
               bool debug, const std::experimental::filesystem::path &dir,
               logging::Logger &logger)
      : rr_paths{rr_paths} {
    const std::unordered_map<std::string, uint64_t> vert2ind = [&dbg, this]() {
      std::unordered_map<std::string, uint64_t> vert2ind;
      for (const Vertex &vertex : dbg.vertices()) {
        const std::string &id = vertex.getId();
        vert2ind.emplace(id, max_vert_index);
        ++max_vert_index;
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

      RREdgeProperty edge_property{max_edge_index, std::move(seq),
                                   classificator.isUnique(edge)};
      add_edge_with_prop(start_ind, end_ind, std::move(edge_property));
      ++max_edge_index;
    }
    assert_validity();
  }

  void serialize_to_dot(const std::experimental::filesystem::path &path) const {
    graph_lite::Serializer serializer(*this);
    std::ofstream dot_os(path);
    serializer.serialize_to_dot(dot_os);
  }

  void inc() {
    if (is_frozen()) {
      return;
    }
    const std::vector<RRVertexType> vertexes = [this]() {
      std::vector<RRVertexType> vertexes;
      for (auto &v : *this) {
        auto &vertex = (RRVertexType &)v;
        vertexes.emplace_back(v);
      }
      return vertexes;
    }();
    for (const auto &vertex : vertexes) {
      process_vertex(vertex);
    }
    finalize_transformation();
    assert_validity();
    ++niter;
  }

  void incN(const uint64_t n_iter) {
    for (uint64_t i = 0; i < n_iter; ++i) {
      inc();
    }
  }
};

} // End namespace repeat_resolution