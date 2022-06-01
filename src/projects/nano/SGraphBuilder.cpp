//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#include <omp.h>
#include <common/omp_utils.hpp>

#include "edlib/edlib.h"
#include "SGraphBuilder.h"

using namespace nano;

void SGraphBuilder::PrintSGraph() {
    for (auto const &[e1, edges]: sgraph_) {
        for (auto const &[e2, val]: edges) {
            std::cerr << e1 << " " << e2 << " " << val << std::endl;
        }
    }
}

void SGraphBuilder::SaveSGraph(const std::experimental::filesystem::path &sgraph_filename) {
    std::ofstream out_file;
    out_file.open(sgraph_filename);
    for (auto const &[e1, edges]: sgraph_) {
        for (auto const &[e2, val]: edges) {
            std::string subpath = "";
            for (auto const edge_id: val.second) {
                subpath += std::to_string(edge_id) + ",";
            }
            out_file << e1 << "\t" << e2 << "\t" << val.first <<
                              "\t" << subpath.substr(0, subpath.size()-1) << std::endl;
        }
    }
    out_file.close();
    std::cerr << "SGraph saved to " << sgraph_filename << std::endl;
}


void SGraphBuilder::LoadAlignments(const std::unordered_map<std::string, nano::GraphContig> &alignments,
                                   const size_t threads) {
    std::vector<std::string> names;
    for (auto const &[key, val]: alignments) {
        names.push_back(key);
    }
    #pragma omp parallel for num_threads(threads)
    for (auto const name: names) {
        const std::string key = name;
        const nano::GraphContig &val = alignments.at(key);
        int prev_index = -1;
        for (int i = 0; i < val.path.size(); ++ i) {
            const std::string &edge_id = val.path[i];
            if (uedges_.count(edge_id.substr(0, edge_id.size() - 1)) == 1) {
                if (prev_index != -1) {
                    std::string e1 = val.path[prev_index];
                    std::string e2 = edge_id;
                    std::vector<std::string> subpath = {val.path.begin() + prev_index, val.path.begin() + i + 1};
                    //std::cerr << "Check " << e1 << " " << e2 << " " << key << std::endl;
                    val.PrintContig();
                    if (by_vertex_) {
                        std::tie(e1, e2) = ResolveByVertex(subpath, val.read_str);
                        //std::cerr << " Add " << e1 << " " << e2 <<  " " << key << std::endl;
                        //std::cerr << std::endl;
                    }
#pragma omp critical(add_edge)
                    if (e1 != "" && e2 != "") {
                        AddEdge(e1, e2, subpath);
                    }
                }
                prev_index = i;
            }
        }
    }
}

void SGraphBuilder::LoadSGraphEdges(const std::experimental::filesystem::path &sgraph_filename) {
    sgraph_.clear();
    std::ifstream in_file;
    in_file.open(sgraph_filename);
    std::string ln;
    while (std::getline(in_file, ln)) {
        //std::cerr << ln << std::endl;
        std::istringstream iss(ln);
        std::string s;
        char delim = '\t';
        std::vector<std::string> params;
        while (std::getline(iss, s, delim)) {
            params.push_back(s);
        }
        int e1 = atoi(params[0].c_str());
        int e2 = atoi(params[1].c_str());
        int cnt = atoi(params[2].c_str());
        std::vector<int> edges;
        if (params.size() > 3) {
            std::istringstream subpath_iss(params[3]);
            //std::cerr << params[3] << std::endl;
            delim = ',';
            while (std::getline(subpath_iss, s, delim)) {
                edges.push_back(atoi(s.c_str()));
            }
        }
        sgraph_[e1][e2] = std::pair<int, std::vector<int>> (cnt, edges);
        //std::cerr << e1 << " " << e2 << " " << cnt << " " << edges.size() << std::endl;
    }
    in_file.close();
    std::cerr << "SGraph uploaded from " << sgraph_filename << std::endl;
}

const multigraph::Edge *GetEdgebyStr(const std::string &edge_id_str, const multigraph::MultiGraph &mg_) {
    int edge_id = atoi(edge_id_str.substr(0, edge_id_str.size() - 1).c_str());
    edge_id = edge_id_str[edge_id_str.size() - 1] == '-'? -edge_id: edge_id;
    return &mg_.edges.at(edge_id);
}

void SGraphBuilder::AddEdge(const std::string &prev_edge_id_str, const std::string &edge_id_str, const std::vector<std::string> &subpath) {
    const multigraph::Edge *prev_edge_id = GetEdgebyStr(prev_edge_id_str, mg_);
    const multigraph::Edge *edge_id = GetEdgebyStr(edge_id_str, mg_);
    std::vector<int> subpath_edge_id;
    for (int i = 1; i < subpath.size() - 1;  ++ i) {
        const std::string &edge_id_str = subpath[i];
        subpath_edge_id.push_back(GetEdgebyStr(edge_id_str, mg_)->getId());
    }
    if (sgraph_.count(prev_edge_id->getId()) == 0) {
        sgraph_[prev_edge_id->getId()] = std::unordered_map<int, std::pair<int, std::vector<int>>>();
    }
    if (sgraph_[prev_edge_id->getId()].count(edge_id->getId()) == 0) {
        sgraph_[prev_edge_id->getId()][edge_id->getId()].first = 0;
    }
    ++ sgraph_[prev_edge_id->getId()][edge_id->getId()].first;
    sgraph_[prev_edge_id->getId()][edge_id->getId()].second = subpath_edge_id;

    const multigraph::Edge *prev_edge_id_rc = prev_edge_id->rc;
    const multigraph::Edge *edge_id_rc = edge_id->rc;
    std::vector<int> subpath_edge_id_rc;
    for (auto it = subpath_edge_id.rbegin(); it != subpath_edge_id.rend(); ++it) {
        subpath_edge_id_rc.push_back(mg_.edges.at(*it).rc->getId());
    }

    if (sgraph_.count(edge_id_rc->getId()) == 0) {
        sgraph_[edge_id_rc->getId()] = std::unordered_map<int, std::pair<int, std::vector<int>>>();
    }
    if (sgraph_[edge_id_rc->getId()].count(prev_edge_id_rc->getId()) == 0) {
        sgraph_[edge_id_rc->getId()][prev_edge_id_rc->getId()].first = 0;
    }
    ++ sgraph_[edge_id_rc->getId()][prev_edge_id_rc->getId()].first;
    sgraph_[edge_id_rc->getId()][prev_edge_id_rc->getId()].second = subpath_edge_id_rc;
}

const multigraph::Edge* SGraphBuilder::GetEdgeByNuc(const char nuc, const multigraph::Edge *edge, bool is_inEdge) {
    if (is_inEdge) {
        std::unordered_map<char, char> nuc_map = {
                {'A', 'T'},
                {'C', 'G'},
                {'T', 'A'},
                {'G', 'C'},
        };
        const multigraph::Vertex *v = edge->end->rc;
        for (auto const out_edge: v->outgoing) {
            int v_len = v->size();
            char out_edge_nuc = out_edge->getSeq().str()[v_len];
            if (out_edge_nuc == nuc_map[nuc]) {
                return out_edge->rc;
            }
        }
    } else {
        const multigraph::Vertex *v = edge->start;
        for (auto const out_edge: v->outgoing) {
            int v_len = v->size();
            char out_edge_nuc = out_edge->getSeq().str()[v_len];
            if (out_edge_nuc == nuc) {
                return out_edge;
            }
        }
    }
    return nullptr;
}

std::pair<std::string, std::string> SGraphBuilder::ResolveByVertex(const std::vector<std::string> &subpath, const Contig &read_str) {
    std::string path_str = "";
    if (subpath.size() == 2) {
        const multigraph::Edge *cur_graph_edge = GetEdgebyStr(subpath[0], mg_);
        path_str = cur_graph_edge->end->seq.str();
    } else {
        for (int i = 1; i < subpath.size() - 1; ++ i) {
            const std::string &cur_edge = subpath[i];
            const multigraph::Edge *cur_graph_edge = GetEdgebyStr(cur_edge, mg_);
            int v_len = cur_graph_edge->end->size();
            path_str += cur_graph_edge->getSeq().str().substr(0, cur_graph_edge->size() - v_len);
            //std::cerr << path_str.size() << " " << cur_graph_edge->getId() << " " << v_len << std::endl;
        }
        const std::string &cur_edge = subpath[subpath.size() - 1];
        const multigraph::Edge *cur_graph_edge = GetEdgebyStr(cur_edge, mg_);
        int v_len = cur_graph_edge->start->size();
        path_str += cur_graph_edge->getSeq().str().substr(0, v_len);
        //std::cerr << path_str.size() << " " << cur_graph_edge->getId() << " " << v_len << std::endl;
    }
    EdlibAlignResult result = edlibAlign(path_str.c_str(), path_str.size(),
                                         read_str.str().c_str(), read_str.str().size(),
                                         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    std::string in_edge_id = "";
    std::string out_edge_id = "";
    if (result.status == EDLIB_STATUS_OK) {
        int start_pos = result.startLocations[0];
        int end_pos = result.endLocations[0];
        char prev_nuc = read_str.str()[start_pos - 1];
        char post_nuc = read_str.str()[end_pos + 1];
//        std::cerr << path_str.size() << " " << read_str.str().size() << " "
//                 << start_pos << " " << end_pos << " "
//                 << result.editDistance << " nuc1:" << prev_nuc
//                 << " nuc2:" << post_nuc << " " << result.numLocations << std::endl;
        const multigraph::Edge *path_in_edge = GetEdgebyStr(subpath[0], mg_);
        const multigraph::Edge *in_edge = GetEdgeByNuc(prev_nuc, path_in_edge, true);
        const multigraph::Edge *path_out_edge = GetEdgebyStr(subpath[subpath.size() - 1], mg_);
        const multigraph::Edge *out_edge = GetEdgeByNuc(post_nuc, path_out_edge, false);
        if (in_edge != nullptr) {
            in_edge_id = in_edge->getLabel() + (in_edge->isCanonical()? "+" : "-");
        }
        if (out_edge != nullptr) {
            out_edge_id = out_edge->getLabel() + (out_edge->isCanonical()? "+" : "-");
        }
    }
    edlibFreeAlignResult(result);
    return std::pair<std::string, std::string>(in_edge_id, out_edge_id);
}

void SGraphBuilder::LoadUEdges(const std::experimental::filesystem::path &unique_edges) {
    uedges_.clear();
    std::ifstream is_cut;
    is_cut.open(unique_edges);
    std::string ln;
    while (std::getline(is_cut, ln)) {
        uedges_.insert(ln);
    }
    std::cerr << "Unique loaded" << std::endl;
    bool changed = true;
    while (changed) {
        changed = false;
        std::unordered_set<std::string> to_remove;
        for (auto e_str: uedges_) {
            const multigraph::Edge *edge = GetEdgebyStr(e_str + "+", mg_);
            bool is_unique = true;
            for (auto e_out: edge->start->outgoing) {
                std::string e_out_str = std::to_string(abs(e_out->getId()));
                if (uedges_.count(e_out_str) == 0) {
                    is_unique = false;
                }
            }
            for (auto e_in: edge->rc->start->outgoing) {
                std::string e_in_str = std::to_string(abs(e_in->getId()));
                if (uedges_.count(e_in_str) == 0) {
                    is_unique = false;
                }
            }
            if (not is_unique) {
                std::cerr << "erase from unique " << e_str << std::endl;
                to_remove.insert(e_str);
                changed = true;
            }
        }
        for (auto e: to_remove) {
            uedges_.erase(e);
        }
    }
    std::cerr << "Unique filtered" << std::endl;
}
