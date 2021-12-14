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
  uint64_t start_k{1};
  bool contains_rc = true;

  static std::vector<SuccinctEdgeInfo>
  SparseDBG2SuccinctEdgeInfo(dbg::SparseDBG &dbg,
                             const UniqueClassificator &classificator);

  void AssertValidity() const;

  void SpreadFrost();
  void FreezeUnpairedVertices();

public:
  MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges, uint64_t start_k,
               RRPaths *rr_paths, bool contains_rc);

  MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *rr_paths, uint64_t start_k,
               UniqueClassificator &classificator);

  MultiplexDBG(const MultiplexDBG &) = delete;
  MultiplexDBG(MultiplexDBG &&) = default;
  MultiplexDBG &operator=(const MultiplexDBG &) = delete;
  MultiplexDBG &operator=(MultiplexDBG &&) = default;

  void SerializeToDot(const std::experimental::filesystem::path &path) const;

  [[nodiscard]] bool IsFrozen() const;

  [[nodiscard]] bool IsVertexComplex(const RRVertexType &vertex) const;

  [[nodiscard]] bool IsVertexSimple(const RRVertexType &vertex) const;

  [[nodiscard]] bool IsVertexCanonical(const RRVertexType &vertex) const;
  [[nodiscard]] bool IsEdgeCanonical(ConstIterator vertex,
                                     NeighborsConstIterator e_it) const;

  void FreezeVertex(const RRVertexType &vertex) { node_prop(vertex).Freeze(); }

  RRVertexType GetNewVertex(MDBGSeq seq);

  [[nodiscard]] size_t FullEdgeSize(ConstIterator vertex,
                                    NeighborsConstIterator e_it) const;

  MDBGSeq ExtractEdgePostStartPrefix(ConstIterator vertex,
                                     NeighborsIterator e_it, uint64_t len);

  MDBGSeq ExtractEdgePreEndSuffix(ConstIterator vertex, NeighborsIterator e_it,
                                  uint64_t len);

  void IncreaseVertex(const RRVertexType &vertex, uint64_t len);

  void MoveEdge(const RRVertexType &s1, NeighborsIterator e1_it,
                const RRVertexType &s2, const RRVertexType &e2);

  void MergeEdges(const RRVertexType &s1, NeighborsIterator e1_it,
                  NeighborsIterator e2_it);

  RREdgeIndexType AddConnectingEdge(NeighborsIterator e1_it,
                                    const RRVertexType &s2,
                                    NeighborsIterator e2_it);

  [[nodiscard]] std::vector<RREdgeIndexType>
  GetInEdgesIndexes(const RRVertexType &vertex) const;
  [[nodiscard]] std::vector<RREdgeIndexType>
  GetOutEdgesIndexes(const RRVertexType &vertex) const;

  std::pair<std::vector<RREdgeIndexType>, std::vector<RREdgeIndexType>>
  GetNeighborEdgesIndexes(const RRVertexType &vertex) const;

  using EdgeNeighborMap =
      std::unordered_map<RREdgeIndexType, std::unordered_set<RREdgeIndexType>>;
  [[nodiscard]] std::pair<EdgeNeighborMap, EdgeNeighborMap>
  GetEdgepairsVertex(const RRVertexType &vertex) const;

  NeighborsIterator FindInEdgeIterator(const RRVertexType &v,
                                       const RREdgeIndexType &edge);
  [[nodiscard]] NeighborsConstIterator
  FindInEdgeConstiterator(const RRVertexType &v,
                          const RREdgeIndexType &edge) const;

  NeighborsIterator FindOutEdgeIterator(const RRVertexType &v,
                                        const RREdgeIndexType &edge);
  [[nodiscard]] NeighborsConstIterator
  FindOutEdgeConstiterator(const RRVertexType &v,
                           const RREdgeIndexType &edge) const;

  [[nodiscard]] int64_t GetInnerEdgeSize(ConstIterator vertex,
                                         NeighborsConstIterator e_it) const;
  [[nodiscard]] MDBGSeq GetEdgeSequence(ConstIterator vertex,
                                        NeighborsConstIterator e_it,
                                        bool trim_left, bool trim_right) const;
  [[nodiscard]] std::vector<Contig>
  GetTrimEdges(int64_t min_inner_edge_size = 1000) const;
  [[nodiscard]] std::vector<Contig>
  PrintTrimEdges(const std::experimental::filesystem::path &f) const;
};

} // End namespace repeat_resolution