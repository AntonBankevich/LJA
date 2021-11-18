//
// Created by Andrey Bzikadze on 11/09/21.
//

#pragma once

#include "multiplex_dbg_topology.hpp"
#include <cctype>
#include <dbg/graph_alignment_storage.hpp>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace repeat_resolution {

using PathEdgeList = std::list<EdgeIndexType>;

struct RRPath {
  // TODO add invariant that no two mappings can have the same id
  std::string id;
  PathEdgeList edge_list;
};

bool operator==(const RRPath &lhs, const RRPath &rhs);

struct IteratorInPath {
  // We can change the mapping from here
  RRPath *const p_path{nullptr};
  const PathEdgeList::const_iterator iter{};

  IteratorInPath(RRPath *const p_path, const PathEdgeList::const_iterator &iter)
      : p_path{p_path}, iter{iter} {}
};

bool operator==(const IteratorInPath &lhs, const IteratorInPath &rhs);

struct IteratorInPathHash {
  // this assumes that iterator can be dereferenced!
  std::size_t operator()(const IteratorInPath &iter_in_mapping) const noexcept {
    std::size_t h1 = std::hash<RRPath *>{}(iter_in_mapping.p_path);
    std::size_t h2 = std::hash<const EdgeIndexType *>{}(&*iter_in_mapping.iter);
    return h1 ^ (h2 << 1);
  }
};

using IteratorInPathUSet =
    std::unordered_set<IteratorInPath, IteratorInPathHash>;
using EdgeIndex2PosMap = std::unordered_map<EdgeIndexType, IteratorInPathUSet>;

using PairEdgeIndexType = std::pair<EdgeIndexType, EdgeIndexType>;
struct PairEdgeIndexHash {
  std::size_t
  operator()(const PairEdgeIndexType &pair_edge_index) const noexcept {
    std::size_t h1 = std::hash<EdgeIndexType>{}(pair_edge_index.first);
    std::size_t h2 = std::hash<EdgeIndexType>{}(pair_edge_index.second);
    return h1 ^ (h2 << 1);
  }
};
using EdgeIndexPair2PosMap =
    std::unordered_map<PairEdgeIndexType, IteratorInPathUSet,
                       PairEdgeIndexHash>;

class RRPaths {
  std::vector<RRPath> paths;
  EdgeIndex2PosMap edge2pos;
  EdgeIndexPair2PosMap edgepair2pos;

public:
  RRPaths(std::vector<RRPath> paths, EdgeIndex2PosMap edge2pos,
          EdgeIndexPair2PosMap edgepair2pos)
      : paths{std::move(paths)}, edge2pos{std::move(edge2pos)},
        edgepair2pos{std::move(edgepair2pos)} {
    assert_validity();
  }

  void assert_validity() const;

  const std::vector<RRPath> &GetPaths() const;
  const EdgeIndex2PosMap &GetEdge2Pos() const;
  const EdgeIndexPair2PosMap &GetEdgepair2Pos() const;

  RRPaths() = default;
  RRPaths(const RRPaths &) = delete;
  RRPaths(RRPaths &&) = default;
  RRPaths &operator=(const RRPaths &) = delete;
  RRPaths &operator=(RRPaths &&) = default;

  void add(EdgeIndexType left_index, EdgeIndexType right_index,
           EdgeIndexType new_index);

  void remove(EdgeIndexType index);

  void merge(EdgeIndexType left_index, EdgeIndexType right_index);
};

class PathsBuilder {
public:
  static RRPaths FromPathVector(std::vector<RRPath> path_vec) {
    EdgeIndex2PosMap edge2pos;
    EdgeIndexPair2PosMap edgepair2pos;
    for (RRPath &path : path_vec) {
      for (auto it = path.edge_list.begin(); it != path.edge_list.end(); ++it) {
        edge2pos[*it].emplace(&path, it);
      }

      // TODO make a normal zip
      for (auto it2{path.edge_list.begin()}, it1{it2++};
           it2 != path.edge_list.end(); ++it1, ++it2) {
        edgepair2pos[std::make_pair(*it1, *it2)].emplace(&path, it1);
      }
    }
    return {std::move(path_vec), std::move(edge2pos), std::move(edgepair2pos)};
  }

  static RRPaths
  FromStorages(const std::vector<RecordStorage *> &storages,
               const std::unordered_map<std::string, size_t> &edgeid2ind) {
    std::vector<RRPath> paths;
    for (RecordStorage *const storage : storages) {
      if (storage == nullptr) {
        continue;
      }
      for (const AlignedRead &aligned_read : *storage) {
        const dbg::Path dbg_path = aligned_read.path.getPath();
        PathEdgeList edge_list;
        for (const Edge *p_edge : dbg_path) {
          EdgeIndexType edge_i = edgeid2ind.at(p_edge->getId());
          edge_list.emplace_back(edge_i);
        }
        paths.push_back({aligned_read.id, std::move(edge_list)});
      }
    }
    return FromPathVector(std::move(paths));
  }

  static RRPaths FromDBGStorages(dbg::SparseDBG &dbg,
                                 const std::vector<RecordStorage *> &storages) {
    std::unordered_map<std::string, size_t> edgeid2ind;
    size_t i = 0;
    for (auto it = dbg.edges().begin(); it != dbg.edges().end(); ++it) {
      const Edge &edge = *it;
      // TODO use it - dbg.edges.begin()
      edgeid2ind[edge.getId()] = i;
      ++i;
    }
    return PathsBuilder::FromStorages(storages, edgeid2ind);
  }
};

} // End namespace repeat_resolution