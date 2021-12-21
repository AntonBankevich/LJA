//
// Created by Andrey Bzikadze on 11/09/21.
//

#pragma once

#include "mdbg_topology.hpp"
#include <cctype>
#include <dbg/graph_alignment_storage.hpp>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace repeat_resolution {

using PathEdgeList = std::list<RREdgeIndexType>;

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

    IteratorInPath(RRPath *const p_path,
                   const PathEdgeList::const_iterator &iter)
        : p_path{p_path}, iter{iter} {}
};

bool operator==(const IteratorInPath &lhs, const IteratorInPath &rhs);

struct IteratorInPathHash {
    // this assumes that iterator can be dereferenced!
    std::size_t operator()(const IteratorInPath &iter_in_mapping) const noexcept {
        std::size_t h1 = std::hash<RRPath *>{}(iter_in_mapping.p_path);
        std::size_t
            h2 = std::hash<const RREdgeIndexType *>{}(&*iter_in_mapping.iter);
        return h1 ^ (h2 << 1);
    }
};

using IteratorInPathUSet =
std::unordered_set<IteratorInPath, IteratorInPathHash>;
using EdgeIndex2PosMap = std::unordered_map<RREdgeIndexType,
                                            IteratorInPathUSet>;

using PairEdgeIndexType = std::pair<RREdgeIndexType, RREdgeIndexType>;
struct PairEdgeIndexHash {
    std::size_t
    operator()(const PairEdgeIndexType &pair_edge_index) const noexcept {
        std::size_t h1 = std::hash<RREdgeIndexType>{}(pair_edge_index.first);
        std::size_t h2 = std::hash<RREdgeIndexType>{}(pair_edge_index.second);
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

    void Add(RREdgeIndexType left_index, RREdgeIndexType right_index,
             RREdgeIndexType new_index);

    void Remove(RREdgeIndexType index);

    void Merge(RREdgeIndexType left_index, RREdgeIndexType right_index);

    [[nodiscard]] bool ContainsPair(const RREdgeIndexType &lhs,
                                    const RREdgeIndexType &rhs) const;

    [[nodiscard]] std::vector<std::pair<RREdgeIndexType, RREdgeIndexType>>
    GetActiveTransitions() const;

    void ExportActiveTransitions(const std::experimental::filesystem::path &path) const;
};

class PathsBuilder {
 public:
    static RRPaths FromPathVector(std::vector<RRPath> path_vec);

    static RRPaths
    FromStorages(const std::vector<RecordStorage *> &storages,
                 const std::unordered_map<std::string, size_t> &edgeid2ind);

    static RRPaths FromDBGStorages(dbg::SparseDBG &dbg,
                                   const std::vector<RecordStorage *> &storages);
};

} // End namespace repeat_resolution