//
// Created by Tatiana Dvorkina on 05.07.2022.
//

#include "edlib/edlib.h"

#include "TipResolver.h"

using namespace nano;

bool TipResolver::AlignedOnTipVertex(const std::vector<GraphContig> &r_alignments) {
    for (auto gcontig: r_alignments) {
        for (auto const e_str: gcontig.path) {
            const multigraph::Edge *edge = GetEdgebyStr(e_str, mg_);
            if (tipvertices_.count(edge->end->id) == 1 || tipvertices_.count(edge->start->id) == 1) {
                return true;
            }
        }
    }
    return false;
}

void TipResolver::LoadPaths(std::unordered_map<std::string, std::vector<GraphContig>> &alignments){
    for (auto const &[key, val]: alignments) {
        if (AlignedOnTipVertex(val)) {
            std::cerr << "Tipped path " << key << std::endl;
            tipReads_[key] = val;
        }
    }
    return;
}

std::unordered_set<std::string> TipResolver::IdentifyConnectedVertices(const std::vector<GraphContig> &r_alignments){
    std::unordered_set<std::string> connected_vertices;
    for (auto const gcontig: r_alignments) {
        if (gcontig.qEnd - gcontig.qStart > 0.9*gcontig.qLen) {
            gcontig.PrintContig();
            std::cerr << gcontig.qEnd - gcontig.qStart << " " << 0.3*gcontig.qLen << std::endl;
            for (auto const e_str: gcontig.path) {
                const multigraph::Edge *edge = GetEdgebyStr(e_str, mg_);
                if (tipvertices_.count(edge->start->id) == 1 && tipvertices_.count(edge->end->rc->id) == 1) {
                    std::string key = std::to_string(edge->start->id) + "_" + std::to_string(edge->end->id);
                    // std::cerr << edge->getId() << std::endl;
                    // std::cerr << key << std::endl;
                    connected_vertices.insert(key);
                    std::string key_rc = std::to_string(edge->rc->start->id) + "_" + std::to_string(edge->rc->end->id);
                    //std::cerr << key_rc << std::endl;
                    connected_vertices.insert(key_rc);
                }
            }
        }
    }

    for (int i = 0; i < r_alignments.size() - 1; ++ i) {
        for (int j =  i + 1; j < r_alignments.size(); ++ j) {
            const GraphContig gcontig_in = r_alignments[i];
            const GraphContig gcontig_out = r_alignments[j];
            if (gcontig_in.qEnd - gcontig_in.qStart <= 0.3*gcontig_in.qLen || 
               gcontig_out.qEnd - gcontig_out.qStart <= 0.3*gcontig_out.qLen) continue;
            if (gcontig_in.qEnd <= gcontig_out.qStart) {
                const multigraph::Edge *s_edge = GetEdgebyStr(gcontig_in.path[gcontig_in.path.size() - 1], mg_);
                const multigraph::Edge *e_edge = GetEdgebyStr(gcontig_out.path[0], mg_);
                if (s_edge->getId() != e_edge->getId() && 
                    s_edge->end->outDeg() == 0 && e_edge->start->inDeg() == 0 &&
                    tipvertices_.count(s_edge->start->id)*tipvertices_.count(e_edge->end->rc->id) == 1) {
                    gcontig_in.PrintContig();
                    gcontig_out.PrintContig();
                    std::string key = std::to_string(s_edge->start->id) + "_" + std::to_string(e_edge->end->id);
                    std::cerr << s_edge->getId() << " " << e_edge->getId() << std::endl;
                    std::cerr << key << std::endl;
                    connected_vertices.insert(key);
                    std::string key_rc = std::to_string(e_edge->rc->start->id) + "_" + std::to_string(s_edge->rc->end->id);
                    std::cerr << key_rc << std::endl;
                    connected_vertices.insert(key_rc);
                }
            }
        }
    }
    VERIFY(connected_vertices.size() == 2 || connected_vertices.size() == 0)
    return connected_vertices;
}

std::pair<int, int> TipResolver::GetBestEdgePair(const std::string &key, int v1, int v2, bool rc) {
    std::string read_str = tipReads_[key][0].read_str.str();
    if (rc) {
        std::cerr << "Try rc" << std::endl;
        read_str = Sequence(tipReads_[key][0].read_str.str(), true).str();
    }
    int v1_len = mg_.vertices.at(v1).size();
    const int SHIFT = 10;
    
    int min_len_in = mg_.vertices.at(v1).outgoing.at(0)->size();
    bool in_equal = false;
    for (auto e: mg_.vertices.at(v1).outgoing) {
        std::cerr << e->getId() << std::endl;
        min_len_in = std::min((int) e->size(), min_len_in);
    }
    int ALNSIZE = std::min(min_len_in, 20000);
    int best_in_score = min_len_in;
    int best_in_edge = -1;
    for (auto e: mg_.vertices.at(v1).outgoing) {
        std::string edge_str = e->getSeq().str().substr(min_len_in - ALNSIZE, ALNSIZE);
        EdlibAlignResult result = edlibAlign(edge_str.c_str(), edge_str.size(),
                            read_str.c_str(), read_str.size(),
                            edlibNewAlignConfig(int(0.3*edge_str.size()), EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));
        if (result.status == EDLIB_STATUS_OK) {       
            std::cerr << e->getId() << " " << edge_str.size() << " " 
                      << result.editDistance << std::endl;
        } else {
            std::cerr << e->getId() << " " << edge_str.size() << " NOT OK " << std::endl;
        }
        if (result.status == EDLIB_STATUS_OK && result.editDistance > -1 
                    && result.editDistance <= best_in_score) {
            if (result.editDistance == best_in_score) {
                in_equal = true;
            } else {
                in_equal = false;
            }
            best_in_score = result.editDistance;
            best_in_edge = e->getId();
        }
        edlibFreeAlignResult(result);
    }
    std::cerr << key << " In " << min_len_in << " edge=" << best_in_edge << std::endl;
    int v2_len = mg_.vertices.at(v2).size();
    int min_len_out = mg_.vertices.at(v2).rc->outgoing.at(0)->size();

    for (auto e: mg_.vertices.at(v2).rc->outgoing) {
        std::cerr << e->getId() << std::endl;
        min_len_out = std::min((int) e->size(), min_len_out);
    }
    ALNSIZE = std::min(min_len_out, 20000);
    int best_out_score = min_len_out;
    int best_out_edge = -1;
    bool out_equal = false;
    for (auto e: mg_.vertices.at(v2).rc->outgoing) {
        //std::string edge_str = Sequence(e->getSeq().str().substr(v2_len - SHIFT, min_len_out), true).str();
        std::string edge_str_rc = Sequence(e->getSeq().str().substr(min_len_out - ALNSIZE, ALNSIZE), true).str();
        EdlibAlignResult result_rc = edlibAlign(edge_str_rc.c_str(), edge_str_rc.size(),
                            read_str.c_str(), read_str.size(),
                            edlibNewAlignConfig(int(0.3*edge_str_rc.size()), EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));
        if (result_rc.status == EDLIB_STATUS_OK) {       
            std::cerr << e->rc->getId() << " " << edge_str_rc.size() << " " << 
                         result_rc.editDistance << std::endl;
        } else {
            std::cerr << e->rc->getId() << " " << edge_str_rc.size() << " NO OK" << std::endl;
        }
        if (result_rc.status == EDLIB_STATUS_OK  && 
            result_rc.editDistance > -1 && result_rc.editDistance < best_out_score) {
            if (result_rc.editDistance == best_out_score) {
                out_equal = true;
            } else {
                out_equal = false;
            }
            best_out_score = result_rc.editDistance;
            best_out_edge = e->rc->getId();
        }
        edlibFreeAlignResult(result_rc);
    }

    std::cerr << key << " Out " << min_len_out << " edge=" << best_out_edge << std::endl;
    if (out_equal || in_equal) {
        best_in_edge = -1;
        best_out_edge = -1;
    }
    return std::pair<int, int>(best_in_edge, best_out_edge);
}


std::unordered_map<int, std::unordered_map<int,std::vector<std::string> > > TipResolver::GetTransitions(const std::vector<std::string> &count, int v1, int v2) {
    std::unordered_map<int, std::unordered_map<int,std::vector<std::string> > > transitions;
    for (auto const key: count) {
        int best_in_edge = -1, best_out_edge = -1;
        bool rc = false;
        std::tie(best_in_edge, best_out_edge) = GetBestEdgePair(key, v1, v2, rc);
        if (best_in_edge == -1 || best_out_edge == -1) {
            rc = true;
            std::tie(best_in_edge, best_out_edge) = GetBestEdgePair(key, v1, v2, rc);
        }
        if (best_in_edge != -1 && best_out_edge != -1) {
            if (transitions.count(best_in_edge) == 0) {
                transitions[best_in_edge] = std::unordered_map<int,std::vector<std::string> >();
            }
            if (transitions[best_in_edge].count(best_out_edge) == 0) {
                transitions[best_in_edge][best_out_edge] = std::vector<std::string>();
            }

            int best_in_edge_rc = mg_.edges.at(best_in_edge).rc->getId();
            int best_out_edge_rc = mg_.edges.at(best_out_edge).rc->getId();
            if (transitions.count(best_out_edge_rc) == 0) {
                transitions[best_out_edge_rc] = std::unordered_map<int,std::vector<std::string> >();
            }
            if (transitions[best_out_edge_rc].count(best_in_edge_rc) == 0) {
                transitions[best_out_edge_rc][best_in_edge_rc] = std::vector<std::string>();
            }
            if (rc) {
                transitions[best_in_edge][best_out_edge].push_back(key + "_rc");
                transitions[best_out_edge_rc][best_in_edge_rc].push_back(key);
            } else {
                transitions[best_in_edge][best_out_edge].push_back(key);
                transitions[best_out_edge_rc][best_in_edge_rc].push_back(key + "_rc");
            }
        }
    }
    return transitions;
}

std::string TipResolver::GetConsensus(std::vector<std::string> &read_lst
                                      , multigraph::Edge *e_in
                                      , multigraph::Edge *e_out) {
    for (auto key: read_lst) {
        std::string pure_key = key;
        std::cerr << key << " " << key.substr(key.size() - 3, 3) << std::endl;
        multigraph::Edge *cur_e_in = e_in;
        multigraph::Edge *cur_e_out = e_out;
        std::cerr << key << " " << key.substr(key.size() - 3, 3) << " " << cur_e_in->getId() << " " << cur_e_out->getId() << std::endl;
        if (key.substr(key.size() - 3, 3) == "_rc") {
            pure_key = key.substr(0, key.size() - 3);
            cur_e_in = e_out->rc;
            cur_e_out = e_in->rc;
        } 
        std::cerr << key << " " << key.substr(key.size() - 3, 3) << " " << cur_e_in->getId() << " " << cur_e_out->getId() << std::endl;
        std::string read_str = tipReads_[pure_key][0].read_str.str();
        int min_len = 10000;
        std::string edge_str = cur_e_in->getSeq().str().substr(cur_e_in->size() - min_len, min_len);
        EdlibAlignResult result_in = edlibAlign(edge_str.c_str(), edge_str.size(),
                            read_str.c_str(), read_str.size(),
                            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));
        int end_in_pos = result_in.endLocations[0];
        for (int i = 0; i < result_in.numLocations; ++ i) {
            std::cerr << result_in.endLocations[i] << std::endl;
        }
        edge_str = cur_e_out->getSeq().str().substr(0, min_len);
        EdlibAlignResult result_out = edlibAlign(edge_str.c_str(), edge_str.size(),
                            read_str.c_str(), read_str.size(),
                            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));
        int start_out_pos = result_out.startLocations[0];
        for (int i = 0; i < result_out.numLocations; ++ i) {
            std::cerr << result_out.endLocations[i] << std::endl;
        }
        std::cerr << key << " " << key.substr(key.size() - 3, 3) << " " << cur_e_in->getId() << " " << cur_e_out->getId()
            << " Positions: " << end_in_pos << " " << start_out_pos << " " << read_str.size() << std::endl;
        
        std::string read_substr = "";
        if (end_in_pos >= start_out_pos){
            int dif = end_in_pos - start_out_pos;
            std::cerr << "WARNING END_IN_POS >= START_OUT_POS " << dif << std::endl;
            read_substr = e_in->getSeq().str() + 
                            e_out->getSeq().str().substr(dif, e_out->getSeq().str().size() - dif);
        } else {
            read_substr = e_in->getSeq().str() +
                            read_str.substr(end_in_pos, start_out_pos - end_in_pos) +
                            e_out->getSeq().str();
        }
        edlibFreeAlignResult(result_in);
        edlibFreeAlignResult(result_out);
        std::cerr << "Readlen " << read_substr.size() << std::endl;
        return read_substr;
    }
    //TODO
    return "";
}

std::unordered_map<int, int> TipResolver::ResolveWithPaths(){
    std::unordered_map<int, std::unordered_map<int, std::vector<std::string> >> tips_stats;
    std::cerr << "Start to build matrix\n";
    for (auto const &[key, val]: tipReads_) {
        std::cerr << key << std::endl;
        std::unordered_set<std::string> connected_vertices = IdentifyConnectedVertices(val);
        std::unordered_set<std::string> used_connections;
        for (std::string s: connected_vertices) {
            if (used_connections.count(s) == 0) {
                used_connections.insert(s);
                std::string delimiter = "_";
                size_t pos = s.find(delimiter);
                int v1_id = std::atoi(s.substr(0, pos).c_str());
                s.erase(0, pos + delimiter.length());
                int v2_id = std::atoi(s.c_str());
                std::pair<int, int> connected_pair = std::pair<int, int>(v1_id, v2_id);
                if (tips_stats.count(connected_pair.first) == 0) {
                    tips_stats[connected_pair.first] = std::unordered_map<int, std::vector<std::string>>();
                }
                if (tips_stats[connected_pair.first].count(connected_pair.second) == 0) {
                    tips_stats[connected_pair.first][connected_pair.second] = std::vector<std::string>();
                }
                tips_stats[connected_pair.first][connected_pair.second].push_back(key);
                std::string s_rc = std::to_string(mg_.vertices.at(v2_id).rc->id) + "_" + std::to_string(mg_.vertices.at(v1_id).rc->id);
                used_connections.insert(s_rc);
                int v1_id_rc = mg_.vertices.at(v2_id).rc->id;
                int v2_id_rc = mg_.vertices.at(v1_id).rc->id;
                connected_pair = std::pair<int, int>(v1_id_rc, v2_id_rc);
                if (tips_stats.count(connected_pair.first) == 0) {
                    tips_stats[connected_pair.first] = std::unordered_map<int, std::vector<std::string>>();
                }
                if (tips_stats[connected_pair.first].count(connected_pair.second) == 0) {
                    tips_stats[connected_pair.first][connected_pair.second] = std::vector<std::string>();
                }
                tips_stats[connected_pair.first][connected_pair.second].push_back(key);
            }
        }
    }

    std::cerr << "Start resolving\n";
    int MIN_CONNECTION = 1;
    std::unordered_map<int, int> old2new;
    std::set<int> used_vertices;
    for (auto const &[v1, val]: tips_stats) {
        for (auto const &[v2, count]: val) {
            std::cerr << "Resolve " << v1 << " " << v2 << " " << count.size() << std::endl;
            if ((mg_.vertices.at(v1).outgoing.size() == 1 || mg_.vertices.at(v2).rc->outgoing.size() == 1) ||
               (mg_.vertices.at(v1).outgoing.size() == 2 && mg_.vertices.at(v2).rc->outgoing.size() == 2)) {
                int num_connections = mg_.vertices.at(v1).outgoing.size()*mg_.vertices.at(v2).rc->outgoing.size();
                std::cerr << "Has appropriate number of out/in-edges " << num_connections << std::endl;
                if (count.size() > MIN_CONNECTION && 
                    (used_vertices.count(v1) == 0 || used_vertices.count(v2) == 0)) {
                    std::unordered_map<int, std::unordered_map<int,std::vector<std::string>> > transitions = GetTransitions(count, v1, v2);
                    
                    std::vector<std::pair<std::pair<int, int>, int>> transitions_lst;

                    int aligned_reads_num = 0;
                    for (auto const &[e_in, val]: transitions) {
                        for (auto const &[e_out, cnt]: val) {
                            std::cerr << "Trans " << e_in << " " << e_out << " " << cnt.size() << std::endl;
                            aligned_reads_num += cnt.size();
                            transitions_lst.push_back(std::pair<std::pair<int, int>, int>
                                            (std::pair<int, int>(e_in, e_out), cnt.size()));
                        }
                    }
                    std::sort(transitions_lst.begin(), transitions_lst.end(),
                                [](std::pair<std::pair<int, int>, int> a, 
                                   std::pair<std::pair<int, int>, int> b) { return a.second > b.second;});

                    std::cerr << aligned_reads_num << " " << num_connections << " " << aligned_reads_num/num_connections << std::endl;
                    std::unordered_set<int> used;
                    if (mg_.vertices.at(v1).outgoing.size() == 1 || mg_.vertices.at(v2).rc->outgoing.size() == 1) {
                        for (auto const &[e_in_id, val]: transitions) {
                            for (auto const &[e_out_id, cnt]: val) {
                                if (mg_.edges.count(e_in_id) > 0 && mg_.edges.count(e_out_id) > 0
                                    //&& e_in_id != e_out_id
                                    && mg_.edges.at(e_in_id).end->outDeg() == 0 
                                    && mg_.edges.at(e_out_id).rc->end->outDeg() == 0) {
                                    std::cerr << "Consider " << e_in_id << " " << e_out_id 
                                                             << " " <<  mg_.edges.at(e_in_id).end->outDeg()
                                                             << " " << mg_.edges.at(e_out_id).rc->end->outDeg()
                                                             << " " << cnt.size() << std::endl;
                                    if (cnt.size() > MIN_CONNECTION) {
                                        std::cerr << "Consider1 " << e_in_id << " " << e_out_id << " " << cnt.size() << std::endl;
                                        if (mg_.vertices.at(v1).outgoing.size() > 1 && used.count(e_in_id) > 0 ||
                                            mg_.vertices.at(v2).rc->outgoing.size() > 1 && used.count(e_out_id) > 0)
                                            continue;
                                        std::cerr << "Consider2 " << e_in_id << " " << e_out_id << " " << cnt.size() << std::endl;
                                        multigraph::Edge *e_in = &mg_.edges.at(e_in_id);
                                        multigraph::Edge *e_out = &mg_.edges.at(e_out_id);
                                        used.insert(e_in_id);
                                        used.insert(e_out_id);
                                        used.insert(e_in->rc->getId());
                                        used.insert(e_out->rc->getId());
                                        std::cerr << "Add to used " << e_in_id << " " << e_out_id << " "
                                                                    << e_in->rc->getId() << " " << e_out->rc->getId() << std::endl;
                                        std::string consensus = GetConsensus(transitions[e_in_id][e_out_id], e_in, e_out);
                                        Sequence new_edge_seq = Sequence(consensus);
                                        const multigraph::Edge *new_edge = &mg_.addEdge(*e_in->start, *e_out->end, new_edge_seq);
                                        old2new[e_in_id] = new_edge->getId();
                                        old2new[e_out_id] = new_edge->getId();
                                        old2new[e_in->rc->getId()] = new_edge->rc->getId();
                                        old2new[e_out->rc->getId()] = new_edge->rc->getId();
                                        mg_.internalRemoveEdge(e_in_id);
                                        mg_.internalRemoveEdge(e_out_id);
                                        std::cerr << e_in_id << " " << e_out_id << " " << new_edge->getId() << " " << new_edge->size() << std::endl;
                                    }
                                }
                            }
                        }
                    } else {
                        for (auto trans:transitions_lst) {
                            int e_in_id = trans.first.first;
                            int e_out_id = trans.first.second;
                            int cnt = trans.second;
                            if (mg_.edges.count(e_in_id) > 0 && mg_.edges.count(e_out_id) > 0
                                //&& e_in_id != e_out_id
                                && mg_.edges.at(e_in_id).end->outDeg() == 0 
                                && mg_.edges.at(e_out_id).rc->end->outDeg() == 0) {
                                std::cerr << "Consider " << e_in_id << " " << e_out_id 
                                                             << " " <<  mg_.edges.at(e_in_id).end->outDeg()
                                                             << " " << mg_.edges.at(e_out_id).rc->end->outDeg()
                                                             << " " << cnt << std::endl;
                                if (cnt > MIN_CONNECTION) {
                                    std::cerr << "Consider1 " << e_in_id << " " << e_out_id << " " << cnt << std::endl;
                                    if (used.count(e_in_id) > 0 || used.count(e_out_id) > 0) continue;
                                    std::cerr << "Consider2 " << e_in_id << " " << e_out_id << " " << cnt << std::endl;
                                    multigraph::Edge *e_in = &mg_.edges.at(e_in_id);
                                    multigraph::Edge *e_out = &mg_.edges.at(e_out_id);
                                    used.insert(e_in_id);
                                    used.insert(e_out_id);
                                    used.insert(e_in->rc->getId());
                                    used.insert(e_out->rc->getId());
                                    std::cerr << "Add to used " << e_in_id << " " << e_out_id << " "
                                                                << e_in->rc->getId() << " " << e_out->rc->getId() << std::endl;
                                    std::string consensus = GetConsensus(transitions[e_in_id][e_out_id], e_in, e_out);
                                    Sequence new_edge_seq = Sequence(consensus);
                                    const multigraph::Edge *new_edge = &mg_.addEdge(*e_in->start, *e_out->end, new_edge_seq);
                                    old2new[e_in_id] = new_edge->getId();
                                    old2new[e_out_id] = new_edge->getId();
                                    old2new[e_in->rc->getId()] = new_edge->rc->getId();
                                    old2new[e_out->rc->getId()] = new_edge->rc->getId();
                                    mg_.internalRemoveEdge(e_in_id);
                                    mg_.internalRemoveEdge(e_out_id);
                                    std::cerr << e_in_id << " " << e_out_id << " " << new_edge->getId() << " " << new_edge->size() << std::endl;
                                }
                            }
                        }
                    }
                }
                used_vertices.insert(v1);
                used_vertices.insert(v2);
                std::cerr << "Add vertex to used: " << v1 << " " << v2 << std::endl;
                const multigraph::Vertex *v1_rc = mg_.vertices.at(v1).rc;
                const multigraph::Vertex *v2_rc = mg_.vertices.at(v2).rc;
                used_vertices.insert(v1_rc->id);
                used_vertices.insert(v2_rc->id);
                std::cerr << "Add vertex to used: " << v1_rc->id << " " << v2_rc->id << std::endl;
            }
        }
    }
    return old2new;
}

void TipResolver::ExtractTipVertices() {
    tipvertices_.clear();
    for (auto const &[key, val]: mg_.edges) {
        if (val.end->outgoing.size() == 0) {
            tipvertices_.insert(val.start->id);
        }
    }
    std::cerr << "Tip vertices num " << tipvertices_.size() << std::endl; 
}