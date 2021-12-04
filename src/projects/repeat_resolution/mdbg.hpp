//
// Created by Andrey Bzikadze on 11/10/21.
//

#pragma once

#include "error_correction/multiplicity_estimation.hpp"
#include "mdbg_topology.hpp"
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
  friend class MultiplexDBGIncreaser;
  RRPaths *rr_paths;
  uint64_t next_edge_index{0};
  uint64_t next_vert_index{0};
  uint64_t n_iter{0};

  void assert_validity() const;

  void _spread_frost();
  void freeze_unpaired_vertices();

public:
  // This constructor is for testing purposes
  MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges, uint64_t start_k,
               RRPaths *rr_paths);

  /*
  MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *rr_paths, uint64_t start_k,
               UniqueClassificator &classificator, bool debug,
               const std::experimental::filesystem::path &dir,
               logging::Logger &logger);
               */

  MultiplexDBG(const MultiplexDBG &) = delete;
  MultiplexDBG(MultiplexDBG &&) = default;
  MultiplexDBG &operator=(const MultiplexDBG &) = delete;
  MultiplexDBG &operator=(MultiplexDBG &&) = default;

  void serialize_to_dot(const std::experimental::filesystem::path &path) const;

  [[nodiscard]] bool IsFrozen() const;

  RRVertexType GetNewVertex(std::list<char> seq);

  bool is_vertex_complex(const RRVertexType &vertex) const;

  bool is_vertex_simple(const RRVertexType &vertex) const;

  void freeze_vertex(const RRVertexType &vertex) {
    RRVertexProperty &prop = node_prop(vertex);
    prop.freeze();
  }

  size_t FullEdgeSize(ConstIterator vertex, NeighborsConstIterator e_it) const;

  std::list<char> ExtractEdgePostStartPrefix(ConstIterator vertex,
                                             NeighborsIterator e_it,
                                             uint64_t len);

  std::list<char> ExtractEdgePreEndSuffix(ConstIterator vertex,
                                          NeighborsIterator e_it, uint64_t len);

  void IncreaseVertex(const RRVertexType &vertex, uint64_t len);

  void MoveEdge(const RRVertexType &s1, NeighborsIterator e1_it,
                const RRVertexType &s2, const RRVertexType &e2);

  void MergeEdges(const RRVertexType &s1, NeighborsIterator e1_it,
                  NeighborsIterator e2_it);

  EdgeIndexType AddConnectingEdge(NeighborsIterator e1_it,
                                  const RRVertexType &s2,
                                  NeighborsIterator e2_it);

  std::vector<EdgeIndexType>
  get_in_edges_indexes(const RRVertexType &vertex) const;
  std::vector<EdgeIndexType>
  get_out_edges_indexes(const RRVertexType &vertex) const;

  std::pair<std::vector<EdgeIndexType>, std::vector<EdgeIndexType>>
  get_neighbor_edges_indexes(const RRVertexType &vertex) const;

  using EdgeNeighborMap =
      std::unordered_map<EdgeIndexType, std::unordered_set<EdgeIndexType>>;
  std::pair<EdgeNeighborMap, EdgeNeighborMap>
  get_edgepairs_vertex(const RRVertexType &vertex) const;

  NeighborsIterator find_in_edge_iterator(const RRVertexType &v,
                                          const EdgeIndexType &edge);
  NeighborsConstIterator
  find_in_edge_constiterator(const RRVertexType &v,
                             const EdgeIndexType &edge) const;

  NeighborsIterator find_out_edge_iterator(const RRVertexType &v,
                                           const EdgeIndexType &edge);
  NeighborsConstIterator
  find_out_edge_constiterator(const RRVertexType &v,
                              const EdgeIndexType &edge) const;
};

} // End namespace repeat_resolution