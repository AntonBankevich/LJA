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

void RRPaths::SplitByTransition(const RREdgeIndexType &lhs, const RREdgeIndexType &rhs) {
    if (edgepair2pos.find(std::make_pair(lhs, rhs))!=edgepair2pos.end()) {
        std::cerr << "Split by transition " << lhs << " " << rhs << std::endl;
        std::cerr << "Path size: " << paths.size() << std::endl;
        std::vector<RRPath> new_paths;
        std::list<RRPath*> lrpaths;
        std::unordered_set<std::string> paths_with_transition;
        for (const IteratorInPath &iter_in_path : edgepair2pos.at({lhs, rhs})) {
            RRPath *const &p_path = iter_in_path.p_path;
            lrpaths.push_back(p_path);
            paths_with_transition.insert(iter_in_path.p_path->id);
        }
        for (auto p_path: lrpaths) {
            auto &path = *p_path;
            std::vector<PathEdgeList::iterator> split_pairs;
            PathEdgeList::iterator left_iter{path.edge_list.begin()};
            for (auto it = path.edge_list.begin(); it!=path.edge_list.end(); ++it) {
                if (*it == lhs && std::next(it, 1) != path.edge_list.end() && *std::next(it, 1) == rhs) {
                    split_pairs.push_back(left_iter);
                    left_iter = std::next(it, 1);
                }
            }

            split_pairs.push_back(left_iter);
            if (split_pairs.size() > 1) {
                for (auto it = split_pairs.at(1); it!=path.edge_list.end(); ++it) {
                    edge2pos[*it].erase(IteratorInPath(&path, it));
                }
                // TODO make a normal zip
                for (auto it2{split_pairs.at(1)}, it1{it2++}; it2!=path.edge_list.end(); ++it1, ++it2) {
                    edgepair2pos[std::make_pair(*it1, *it2)].erase(IteratorInPath(&path, it1));
                }
                std::reverse(split_pairs.begin(), split_pairs.end());
                for (int i = 0; i < split_pairs.size() - 1; ++ i){
                    PathEdgeList::iterator left_iter = split_pairs[i];
                    PathEdgeList new_edge_list;
                    new_edge_list.splice(new_edge_list.begin(), path.edge_list, left_iter, path.edge_list.end());
                    std::cerr << " Add path: " << path.id + 'R' << std::endl;
                    std::string ids_lst = "";
                    for (auto it = path.edge_list.begin(); it!= path.edge_list.end(); ++it) {
                        ids_lst += std::to_string(*it) + ",";
                    }
                    ids_lst += " ";
                    for (auto it = new_edge_list.begin(); it!=new_edge_list.end(); ++it) {
                        ids_lst += std::to_string(*it) + ",";
                    }
                    std::cerr << " " << ids_lst << std::endl;
                    new_paths.push_back({path.id + 'R', new_edge_list});
                }
            }
        }
	    edgepair2pos.erase({lhs, rhs});
        if (new_paths.size() > 0) {
             for (const auto &path:new_paths) {
                 paths.insert(std::pair<int, RRPath>(paths.size(), path));
                 RRPath &new_path = paths[paths.size() - 1];
            	 for (auto it = new_path.edge_list.begin(); it!=new_path.edge_list.end(); ++it) {
                	edge2pos[*it].emplace(&new_path, it);
            	 }
                    // TODO make a normal zip
            	 for (auto it2{new_path.edge_list.begin()}, it1{it2++};
                        	it2!=new_path.edge_list.end(); ++it1, ++it2) {
                	edgepair2pos[std::make_pair(*it1, *it2)].emplace(&new_path, it1);
            	 }
             }
	     //this->PrintRRPaths();
        }
    }
}

void RRPaths::SplitByTransitions(std::vector<std::pair<RREdgeIndexType, RREdgeIndexType>>
                          &split_pairs) {
    for (auto split_pair: split_pairs) {
        this->SplitByTransition(split_pair.first, split_pair.second);
    }
}

const std::vector<RRPath> &RRPaths::GetPaths() const {
     std::vector<RRPath> paths_vec;
     for (const auto&[key, path]: paths) {
          paths_vec.push_back(path);
     }
     return std::move(paths_vec);
}

void RRPaths::PrintRRPaths() const {
    int cnt = 0;
    int n = 100000;
    // for (const auto &[transition, pos] : edgepair2pos) {
    //      std::cerr << transition.first << " " << transition.second << " " << pos.size() << std::endl;
    //      cnt ++;
    //     if (cnt > n) break;
    // }

    cnt = 0;
    for (const auto &[key,ont_path] : paths) {
        if (ont_path.id.find("e10aa8e58d62") != std::string::npos) {
            std::cerr << ont_path.id << std::endl;
            std::string ids = "";
            for (auto const &edge: ont_path.edge_list) {
                ids += std::to_string(edge) + ",";
            }
            std::cerr << " " << ids << std::endl;
            cnt ++;
            if (cnt > n) break;
        }
    }
    // cnt = 10;
    // for (const auto &[edge, val] : edge2pos) {
    //     std::cerr << edge << " " << val.size() << ": " << std::endl;
    //     int cnt_in = 0;
    //     for (const IteratorInPath &iter_in_path : val) {
    //         std::cerr << " " << iter_in_path.p_path->id << "\n";
    //         cnt_in ++;
    //         if (cnt_in > 10) break;
    //     }
    //     cnt ++;
    //     if (cnt > n) break;
    // }
}

const EdgeIndex2PosMap &RRPaths::GetEdge2Pos() const { return edge2pos; }
const EdgeIndexPair2PosMap &RRPaths::GetEdgepair2Pos() const {
    return edgepair2pos;
}

[[nodiscard]] bool RRPaths::ContainsPair(const RREdgeIndexType &lhs,
                                         const RREdgeIndexType &rhs) const {
    return edgepair2pos.find(std::make_pair(lhs, rhs))!=edgepair2pos.end();
}

[[nodiscard]] int RRPaths::PairCount(const RREdgeIndexType &lhs,
                                         const RREdgeIndexType &rhs) const {
    if (edgepair2pos.find(std::make_pair(lhs, rhs)) == edgepair2pos.end()) {
        return 0;
    } else {
        return edgepair2pos.at({lhs, rhs}).size();
    }
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
    std::unordered_map<int, RRPath> path_map;
    for (RRPath &path : path_vec) {
        path_map.insert(std::pair<int, RRPath> (path_map.size(), path));
    }
    for (auto  &[key, path] : path_map) {
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
    return {std::move(path_map), std::move(edge2pos),
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

nano::GraphContig ExtractAlignment(const std::string &ln, std::unordered_map<std::string, std::string> &edgeid2edgeid_rc){
    std::vector<std::string> params;
    std::istringstream iss(ln);
    std::string s;
    char delim = '\t';
    while (std::getline(iss, s, delim)) {
        params.push_back(s);
    }
    return nano::GraphContig(params, edgeid2edgeid_rc);
}

RRPaths
PathsBuilder::FromGAF(dbg::SparseDBG &dbg,
                           const std::experimental::filesystem::path &batch_gaf) {
    std::unordered_map<std::string, size_t> edgeid2ind;
    std::unordered_map<std::string, std::string> edgeid2edgeid_rc;
    size_t i = 0;
    for (auto it = dbg.edges().begin(); it!=dbg.edges().end(); ++it) {
        const dbg::Edge &edge = *it;
        // TODO use it - dbg.edges.begin()
        edgeid2ind[edge.getId()] = i;
        edgeid2edgeid_rc[edge.getId()] = edge.rc().getId();
        ++i;
    }

    std::vector<RRPath> paths;
    std::ifstream is_cut;
    is_cut.open(batch_gaf);
    std::string ln;
    std::unordered_map<std::string, std::vector<nano::GraphContig>> alignments;
    while (std::getline(is_cut, ln)) {
        nano::GraphContig gcontig = ExtractAlignment(ln, edgeid2edgeid_rc);
        if (gcontig.qEnd - gcontig.qStart > 0.9*gcontig.qLen) {
            if (alignments.count(gcontig.query) == 0) {
                if (alignments.count(gcontig.query) == 0) {
                    alignments[gcontig.query] = std::vector<nano::GraphContig>();
                }
                alignments[gcontig.query].push_back(gcontig);
            } else {
                nano::GraphContig &prevAln = alignments.at(gcontig.query)[0];
                if (prevAln.nMatches * (gcontig.qEnd - gcontig.qStart)
                            < gcontig.nMatches * (prevAln.qEnd - prevAln.qStart)) {
                    alignments.at(gcontig.query)[0] = gcontig;
                }
            }
        }
    }
    is_cut.close();

    for (const auto &[key,val]: alignments) {
        std::vector<std::string> path_ids = val.at(0).path;
        PathEdgeList edge_list;
        PathEdgeList edge_list_rc;
        for (const string &cur_id: path_ids){
            edge_list.emplace_back(edgeid2ind.at(cur_id));
            edge_list_rc.emplace_front(edgeid2ind.at(edgeid2edgeid_rc.at(cur_id)));
        }
        paths.push_back({'+' + key, edge_list});
        paths.push_back({'-' + key, edge_list_rc});
        // std::cerr << std::endl << std::endl << key << std::endl;
        // for (const auto eid: edge_list) {
        //     std::cerr << eid << ",";
        // }
        // std::cerr << "\n";
        // for (const auto eid: edge_list_rc) {
        //     std::cerr << eid << ",";
        // }
    }
    return FromPathVector(std::move(paths));
}
