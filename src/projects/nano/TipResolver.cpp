//
// Created by Tatiana Dvorkina on 05.07.2022.
//

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
            //std::cerr << "Tipped path " << key << std::endl;
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
                    std::string key = std::to_string(s_edge->start->id) + "_" + std::to_string(e_edge->end->id);
                    // std::cerr << s_edge->getId() << " " << e_edge->getId() << std::endl;
                    // std::cerr << key << std::endl;
                    connected_vertices.insert(key);
                    std::string key_rc = std::to_string(e_edge->rc->start->id) + "_" + std::to_string(s_edge->rc->end->id);
                    //std::cerr << key_rc << std::endl;
                    connected_vertices.insert(key_rc);
                }
            }
        }
    }
    VERIFY(connected_vertices.size() == 2 || connected_vertices.size() == 0)
    return connected_vertices;
}

void TipResolver::ResolveWithPaths(){
    std::unordered_map<int, std::unordered_map<int, int>> tips_stats;
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
                    tips_stats[connected_pair.first] = std::unordered_map<int, int>();
                }
                if (tips_stats[connected_pair.first].count(connected_pair.second) == 0) {
                    tips_stats[connected_pair.first][connected_pair.second] == 0;
                }
                ++ tips_stats[connected_pair.first][connected_pair.second];
                std::string s_rc = std::to_string(mg_.vertices.at(v2_id).rc->id) + "_" + std::to_string(mg_.vertices.at(v1_id).rc->id);
                used_connections.insert(s_rc);
                int v1_id_rc = mg_.vertices.at(v2_id).rc->id;
                int v2_id_rc = mg_.vertices.at(v1_id).rc->id;
                connected_pair = std::pair<int, int>(v1_id_rc, v2_id_rc);
                if (tips_stats.count(connected_pair.first) == 0) {
                    tips_stats[connected_pair.first] = std::unordered_map<int, int>();
                }
                if (tips_stats[connected_pair.first].count(connected_pair.second) == 0) {
                    tips_stats[connected_pair.first][connected_pair.second] == 0;
                }
                ++ tips_stats[connected_pair.first][connected_pair.second];
            }
        }
    }

    for (auto const &[v1, val]: tips_stats) {
        for (auto const &[v2, count]: val) {
            std::cerr << v1 << " " << v2 << " " << count << std::endl;
            for (auto e: mg_.vertices.at(v1).outgoing) {
                std::cerr << e->getId() << std::endl;
            }
            for (auto e: mg_.vertices.at(v2).rc->outgoing) {
                std::cerr << e->getId() << std::endl;
            }
        }
    }

    // int MIN_CONNECTION = 1;
    // std::unordered_map<int, int> old2new;
    // for (auto const &[v1, val]: tips_stats) {
    //     for (auto const &[v2, count]: val) {
    //         if (count > MIN_CONNECTION) {
    //             int tip_start_cnt = 0;
    //             for (auto const e: mg_.vertices.at(v1).outgoing) {
    //                 if (e->end->outDeg() == 0) ++ tip_start_cnt; 
    //             }
    //             int tip_end_cnt = 0;
    //             for (auto e: mg_.vertices.at(v2).rc->outgoing) {
    //                 if (e->end->outDeg() == 0) ++ tip_end_cnt; 
    //             }
    //             if (tip_start_cnt == 1 || tip_end_cnt == 1) {
                    
    //             }
    //         }
    //     }
    // }
    // for (auto const &[e1_id, val]: tip2tip) {
    //     for (auto const &[e2_id, connections]: val) {
    //         if (connections.size() > MIN_CONNECTION) {
    //             //TODO Make proper consensus
    //             const multigraph::Edge *e1 = &mg_.edges.at(e1_id);
    //             const multigraph::Edge *e2 = &mg_.edges.at(e2_id);
    //             Sequence new_edge_seq = e1->getSeq() + Sequence(connections[0]) + e2->getSeq();
    //             const multigraph::Edge *new_edge = &mg_.addEdge(*e1->start, *e2->end, new_edge_seq);
    //             old2new[e1_id] = new_edge->getId();
    //             old2new[e2_id] = new_edge->getId();
    //             old2new[e1->rc->getId()] = new_edge->rc->getId();
    //             old2new[e2->rc->getId()] = new_edge->rc->getId();
    //             mg_.internalRemoveEdge(e1_id);
    //             mg_.internalRemoveEdge(e2_id);
    //             std::cerr << e1_id << " " << e2_id << " " << new_edge->getId() << std::endl;
    //         }
    //     }
    // }
    // return old2new;
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