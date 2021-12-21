//
// Created by Andrey Bzikadze on 11/09/21.
//

#include "paths.hpp"
using namespace repeat_resolution;

bool repeat_resolution::operator==(const RRPath &lhs, const RRPath &rhs) {
    return lhs.id==rhs.id and lhs.edge_list==rhs.edge_list;
}

bool repeat_resolution::operator==(const IteratorInPath &lhs,
                                   const IteratorInPath &rhs) {
    return lhs.p_path==rhs.p_path and lhs.iter==rhs.iter;
}

void RRPaths::assert_validity() const {
    for (const auto &pair_ind_pos : edge2pos) {
        auto index = pair_ind_pos.first;
        for (const auto &iter_in_mapping : pair_ind_pos.second) {
            VERIFY(*(iter_in_mapping.iter)==index);
        }
    }
    for (const auto &pair_pairind_pos : edgepair2pos) {
        auto pair_index = pair_pairind_pos.first;
        for (auto iter_in_mapping : pair_pairind_pos.second) {
            VERIFY(*(iter_in_mapping.iter)==pair_index.first);
            auto next_iter = [&iter_in_mapping]() {
              auto next_iter = iter_in_mapping.iter;
              ++next_iter;
              return next_iter;
            }();
            VERIFY(*next_iter==pair_index.second);
        }
    }
}

void RRPaths::Add(RREdgeIndexType left, RREdgeIndexType right,
                  RREdgeIndexType new_index) {
    if (edgepair2pos.find({left, right})==edgepair2pos.end()) {
        return;
    }
    for (const IteratorInPath &iter_in_path : edgepair2pos.at({left, right})) {
        const PathEdgeList::const_iterator left_iter{iter_in_path.iter};
        const auto right_iter = [&iter_in_path]() {
          auto iter = iter_in_path.iter;
          return ++iter;
        }();
        RRPath *const &p_path = iter_in_path.p_path;
        const auto new_iter = p_path->edge_list.insert(right_iter, new_index);
        edge2pos[new_index].emplace(p_path, new_iter);
        edgepair2pos[{left, new_index}].emplace(p_path, left_iter);
        edgepair2pos[{new_index, right}].emplace(p_path, new_iter);
    }
    edgepair2pos.erase({left, right});
}

void RRPaths::Remove(RREdgeIndexType index) {
    if (edge2pos.find(index)==edge2pos.end()) {
        return;
    }
    for (const IteratorInPath &iter_in_path : edge2pos.at(index)) {
        RRPath *const &p_path = iter_in_path.p_path;
        const PathEdgeList::const_iterator &iter = iter_in_path.iter;

        const auto next = [&iter]() {
          auto next = iter;
          return ++next;
        }();

        const bool is_begin = iter==p_path->edge_list.begin();
        const bool is_end = next==p_path->edge_list.end();

        if (not is_end) {
            PairEdgeIndexType mid_next{*iter, *next};
            edgepair2pos.at(mid_next).erase(iter_in_path);
            if (edgepair2pos.at(mid_next).empty()) {
                edgepair2pos.erase(mid_next);
            }
        }
        if (not is_begin) {
            const auto prev = [&iter]() {
              auto prev = iter;
              return --prev;
            }();

            PairEdgeIndexType prev_mid{*prev, *iter};
            IteratorInPath prev_iter_in_path{p_path, prev};
            VERIFY(edgepair2pos.at(prev_mid).find(prev_iter_in_path)!=
                edgepair2pos.at(prev_mid).end());
            edgepair2pos.at(prev_mid).erase(prev_iter_in_path);
            if (edgepair2pos.at(prev_mid).empty()) {
                edgepair2pos.erase(prev_mid);
            }

            PairEdgeIndexType prev_next{*prev, *next};
            if (not is_end) {
                edgepair2pos[prev_next].emplace(prev_iter_in_path);
            }
        }
        p_path->edge_list.erase(iter_in_path.iter);
    }
    edge2pos.erase(index);
}

void RRPaths::Merge(RREdgeIndexType left_index, RREdgeIndexType right_index) {
    // Taken from https://stackoverflow.com/a/3792615
    if (edge2pos.find(right_index)==edge2pos.end()) {
        return;
    }
    IteratorInPathUSet &pos_right_index = edge2pos.at(right_index);
    for (auto it = pos_right_index.begin(); it!=pos_right_index.end();
        /* blank */) {
        const IteratorInPath &iter_in_path = *it;
        const PathEdgeList::const_iterator &list_iter{iter_in_path.iter};
        RRPath *const &p_path{iter_in_path.p_path};

        bool is_begin =
            iter_in_path.iter==iter_in_path.p_path->edge_list.begin();
        if (is_begin) {
            const auto list_next_iter = [&list_iter]() {
              auto next = list_iter;
              return ++next;
            }();
            bool is_end = list_next_iter==p_path->edge_list.end();

            p_path->edge_list.pop_front();
            p_path->edge_list.emplace_front(left_index);

            edge2pos[left_index].emplace(p_path, p_path->edge_list.begin());

            if (not is_end) {
                PairEdgeIndexType right_next{right_index, *list_next_iter};
                edgepair2pos.at(right_next).erase(iter_in_path);
                if (edgepair2pos.at(right_next).empty()) {
                    edgepair2pos.erase(right_next);
                }
                edgepair2pos[{left_index, *list_next_iter}].emplace(
                    p_path, p_path->edge_list.begin());
            }

            pos_right_index.erase(it++);
        } else {
            ++it;
        }
    }
    Remove(right_index);
}

const std::vector<RRPath> &RRPaths::GetPaths() const { return paths; }
const EdgeIndex2PosMap &RRPaths::GetEdge2Pos() const { return edge2pos; }
const EdgeIndexPair2PosMap &RRPaths::GetEdgepair2Pos() const {
    return edgepair2pos;
}

[[nodiscard]] bool RRPaths::ContainsPair(const RREdgeIndexType &lhs,
                                         const RREdgeIndexType &rhs) const {
    return edgepair2pos.find(std::make_pair(lhs, rhs))!=edgepair2pos.end();
}

[[nodiscard]] std::vector<std::pair<RREdgeIndexType, RREdgeIndexType>>
RRPaths::GetActiveTransitions() const {
    std::vector<std::pair<RREdgeIndexType, RREdgeIndexType>> transitions;
    for (const auto &[transition, pos] : edgepair2pos) {
        transitions.push_back(transition);
    }
    return transitions;
}

void RRPaths::ExportActiveTransitions(const std::experimental::filesystem::path &path) const {
    std::vector<std::pair<RREdgeIndexType, RREdgeIndexType>>
        transitions = GetActiveTransitions();
    std::ofstream os(path);
    for(const auto &transition : transitions) {
        os << transition.first << " " << transition.second << "\n";
    }
}

RRPaths PathsBuilder::FromPathVector(std::vector<RRPath> path_vec) {
    EdgeIndex2PosMap edge2pos;
    EdgeIndexPair2PosMap edgepair2pos;
    for (RRPath &path : path_vec) {
        for (auto it = path.edge_list.begin(); it!=path.edge_list.end();
             ++it) {
            edge2pos[*it].emplace(&path, it);
        }

        // TODO make a normal zip
        for (auto it2{path.edge_list.begin()}, it1{it2++};
             it2!=path.edge_list.end(); ++it1, ++it2) {
            edgepair2pos[std::make_pair(*it1, *it2)].emplace(&path, it1);
        }
    }
    return {std::move(path_vec), std::move(edge2pos),
            std::move(edgepair2pos)};
}

RRPaths
PathsBuilder::FromStorages(const std::vector<RecordStorage *> &storages,
                           const std::unordered_map<std::string,
                                                    size_t> &edgeid2ind) {
    std::vector<RRPath> paths;
    auto path2edge_list = [&edgeid2ind](const dbg::Path &dbg_path) {
      PathEdgeList edge_list;
      for (const dbg::Edge *p_edge : dbg_path) {
          RREdgeIndexType edge_i = edgeid2ind.at(p_edge->getId());
          edge_list.emplace_back(edge_i);
      }
      return edge_list;
    };
    for (RecordStorage *const storage : storages) {
        if (storage==nullptr) {
            continue;
        }
        for (const AlignedRead &aligned_read : *storage) {
//            if(not aligned_read.valid()) {
//                continue;
//            }
//            const VertexRecord
//                &rec = storage->getRecord(aligned_read.path.start());
//            size_t cnt = rec.countStartsWith(aligned_read.path.cpath());
//            VERIFY_MSG(cnt >= 1, "This function assumes that suffixes are stored for complete reads")
//            if (rec.countStartsWith(aligned_read.path.cpath()) >= 2) {
//                continue;
//            }
            dbg::Path path = aligned_read.path.getPath();
            if (path.size()==0) {
                continue;
            }
            paths.push_back({'+' + aligned_read.id, path2edge_list(path)});
            paths.push_back({'-' + aligned_read.id,
                             path2edge_list(path.RC())});
        }
    }
    return FromPathVector(std::move(paths));
}

RRPaths PathsBuilder::FromDBGStorages(dbg::SparseDBG &dbg,
                                      const std::vector<RecordStorage *> &storages) {
    std::unordered_map<std::string, size_t> edgeid2ind;
    size_t i = 0;
    for (auto it = dbg.edges().begin(); it!=dbg.edges().end(); ++it) {
        const dbg::Edge &edge = *it;
        // TODO use it - dbg.edges.begin()
        edgeid2ind[edge.getId()] = i;
        ++i;
    }
    return PathsBuilder::FromStorages(storages, edgeid2ind);
}
