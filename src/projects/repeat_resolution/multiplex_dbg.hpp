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
          /*typename NodePropType=*/void,
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
        est_max_vert_index = std::max(est_max_vert_index, vertex.index);
      }
      return est_max_vert_index;
    }();
    VERIFY(max_vert_index == 1 + est_max_vert_index);

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
    VERIFY(max_edge_index == 1 + est_max_edge_index);

    for (const auto &vertex : *this) {
      VERIFY(count_in_neighbors(vertex) != 1 or
             count_out_neighbors(vertex) != 1);
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
      auto [out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
      for (auto in_it = in_nbr_begin; in_it != in_nbr_end; ++in_it) {
        for (auto out_it = out_nbr_begin; out_it != out_nbr_end; ++out_it) {
          in_it->second.prop().assert_incidence(out_it->second.prop(),
                                                vertex.len);
        }
      }
    }
  }

  void move_edge(const RRVertexType &s1, NeighborsIterator e1_it,
                 const RRVertexType &s2, const RRVertexType &e2) {
    add_edge_with_prop(s2, e2, std::move(e1_it->second.prop()));
    ConstIterator s1_it = find(s1);
    remove_edge(s1_it, e1_it);
  }

  void smart_remove_edge(ConstIterator s1_it, NeighborsIterator e1_it,
                         bool moving) {
    const RREdgeProperty &edge_prop = e1_it->second.prop();
    rr_paths->remove(edge_prop.get_index());

    if (moving) {
      auto [in_nbr_begin, in_nbr_end] = in_neighbors(e1_it->first);
      for (auto in_nbr_it = in_nbr_begin; in_nbr_it != in_nbr_end;
           ++in_nbr_it) {
        move_edge(e1_it->first, in_nbr_it, in_nbr_it->first, *s1_it);
      }

      auto [out_nbr_begin, out_nbr_end] = out_neighbors(e1_it->first);
      for (auto out_nbr_it = out_nbr_begin; out_nbr_it != out_nbr_end;
           ++out_nbr_it) {
        move_edge(e1_it->first, out_nbr_it, *s1_it, out_nbr_it->first);
      }
    }
    remove_edge(s1_it, e1_it);
  }

  // void merge_edges()

public:
  // This constructor is for testing purposes
  MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges,
               const uint64_t start_k, RRPaths *const rr_paths)
      : rr_paths{rr_paths} {
    for (const SuccinctEdgeInfo &edge : edges) {
      max_vert_index = std::max(max_vert_index, 1 + edge.start.index);
      max_vert_index = std::max(max_vert_index, 1 + edge.end.index);
      add_nodes(edge.start);
      add_nodes(edge.end);
      RREdgeProperty edge_property{ max_edge_index, edge.seq, edge.unique };
      add_edge_with_prop(edge.start, edge.end, std::move(edge_property));
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
      const uint64_t start_ind = vert2ind.at(edge.start()->getId());
      const uint64_t end_ind = vert2ind.at(edge.end()->getId());
      const RRVertexType start_vt{start_ind, start_k, false};
      const RRVertexType end_vt{end_ind, start_k, false};
      add_nodes(start_vt);
      add_nodes(end_vt);

      std::list<char> seq = [&edge]() {
        std::string seq_str = edge.suffix(0).str();
        std::list<char> seq;
        std::move(seq_str.begin(), seq_str.end(), std::back_inserter(seq));
        return seq;
      }();

      RREdgeProperty edge_property{max_edge_index, std::move(seq),
                                   classificator.isUnique(edge)};
      add_edge_with_prop(start_vt, end_vt, std::move(edge_property));
      ++max_edge_index;
    }
    assert_validity();
  }

  void serialize_to_dot(const std::experimental::filesystem::path &path) const {
    graph_lite::Serializer serializer(*this);
    std::ofstream dot_os(path);
    serializer.serialize_to_dot(dot_os);
  }
};

//
//  void resolve_graph() {
//    std::vector<node_type> nodes_to_remove;
//    for (auto it = begin(); it != end(); ++it) {
//      node_type node{*it};
//      int indegree = count_in_neighbors(node);
//      int outdegree = count_out_neighbors(node);
//      VERIFY(indegree != 1 or outdegree != 1); // a node cannot be
//      on a non - branching path if (node_has_loop(node)) {
//        // TODO. Loops need special care. For now, skip
//      }
//      if (indegree == 0 or outdegree == 0) {
//        /* isolated vertex or a (multi-)tip should be skipped
//         * multi-tip (indegree = 0, outdegree > 1 or vice versa)
//         technically can be resolved
//         * but that will not increase the length of contigs and
//         thus is meaningless.
//         * It will also invalidate iterator unless we switch to
//         ordered containers
//         */
//      } else if (indegree == 1) {
//        auto [ibegin, iend] = in_neighbors(it);
//        auto &[inode, iedge] = *ibegin;
//        auto [obegin, oend] = out_neighbors(it);
//        for (auto oit = obegin; oit != oend; ++oit) {
//          // inode --[iedge]> node --[oedge]> onode
//          auto &[onode, oedge] = *oit;
//          RREdgeProperty new_prop = iedge.prop() + oedge.prop();
//          add_edge_with_prop(inode, onode, new_prop);
//        }
//        nodes_to_remove.emplace_back(std::move(node));
//      } else if (outdegree == 1) {
//        auto [ibegin, iend] = in_neighbors(it);
//        auto [obegin, oend] = out_neighbors(it);
//        auto &[onode, oedge] = *obegin;
//        for (auto iit = ibegin; iit != iend; ++iit) {
//          // inode --[iedge]> node --[oedge]> onode
//          auto &[inode, iedge] = *iit;
//          RREdgeProperty new_prop = iedge.prop() + oedge.prop();
//          add_edge_with_prop(inode, onode, new_prop);
//        }
//        nodes_to_remove.emplace_back(std::move(node));
//      }
//    }
//    for (const node_type &node : nodes_to_remove) {
//      remove_nodes(node);
//    }
//  }
//};

} // End namespace repeat_resolution