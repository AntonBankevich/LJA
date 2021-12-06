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

  void AssertValidity() const;

  void SpreadFrost();
  void FreezeUnpairedVertices();

public:
  // This constructor is for testing purposes
  MultiplexDBG(std::vector<SuccinctEdgeInfo> &edges, uint64_t start_k,
               RRPaths *rr_paths);

  MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *rr_paths, uint64_t start_k,
               UniqueClassificator &classificator, bool debug,
               const std::experimental::filesystem::path &dir,
               logging::Logger &logger);

  MultiplexDBG(const MultiplexDBG &) = delete;
  MultiplexDBG(MultiplexDBG &&) = default;
  MultiplexDBG &operator=(const MultiplexDBG &) = delete;
  MultiplexDBG &operator=(MultiplexDBG &&) = default;

  void SerializeToDot(const std::experimental::filesystem::path &path) const;

  [[nodiscard]] bool IsFrozen() const;

  bool IsVertexComplex(const RRVertexType &vertex) const;

  bool IsVertexSimple(const RRVertexType &vertex) const;

  void FreezeVertex(const RRVertexType &vertex) { node_prop(vertex).freeze(); }

  RRVertexType GetNewVertex(std::list<char> seq);

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
  GetInEdgesIndexes(const RRVertexType &vertex) const;
  std::vector<EdgeIndexType>
  GetOutEdgesIndexes(const RRVertexType &vertex) const;

  std::pair<std::vector<EdgeIndexType>, std::vector<EdgeIndexType>>
  GetNeighborEdgesIndexes(const RRVertexType &vertex) const;

  using EdgeNeighborMap =
      std::unordered_map<EdgeIndexType, std::unordered_set<EdgeIndexType>>;
  std::pair<EdgeNeighborMap, EdgeNeighborMap>
  GetEdgepairsVertex(const RRVertexType &vertex) const;

  NeighborsIterator FindInEdgeIterator(const RRVertexType &v,
                                       const EdgeIndexType &edge);
  NeighborsConstIterator
  FindInEdgeConstiterator(const RRVertexType &v,
                          const EdgeIndexType &edge) const;

  NeighborsIterator FindOutEdgeIterator(const RRVertexType &v,
                                        const EdgeIndexType &edge);
  NeighborsConstIterator
  FindOutEdgeConstiterator(const RRVertexType &v,
                           const EdgeIndexType &edge) const;
};

} // End namespace repeat_resolution