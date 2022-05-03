//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#include "edlib/edlib.h"
#include "SGraphBuilder.h"

using namespace nano;

void SGraphBuilder::PrintSgraph() {
    for (auto const &[e1, edges]: sgraph_) {
        for (auto const &[e2, val]: edges) {
            std::cerr << e1 << " " << e2 << " " << val << std::endl;
        }
    }
}

void SGraphBuilder::LoadAlignments(const std::unordered_map<std::string, nano::GraphContig> &alignments) {
    for (auto const &[key, val]: alignments) {
        int prev_index = -1;
        for (int i = 0; i < val.path.size(); ++ i) {
            const std::string &edge_id = val.path[i];
            if (uedges_.count(edge_id.substr(0, edge_id.size() - 1)) == 1) {
                if (prev_index != -1) {
                    std::string e1 = val.path[prev_index];
                    std::string e2 = edge_id;
                    std::cerr << "Check " << e1 << " " << e2 << std::endl;
                    if (by_vertex_) {
                        std::vector<std::string> subpath = {val.path.begin() + prev_index, val.path.begin() + i + 1};
                        std::tie(e1, e2) = ResolveByVertex(subpath, val.read_str);
                        std::cerr << "Add " << e1 << " " << e2 << std::endl;
                    }
                    if (e1 != "" && e2 != "") {
                        AddEdge(e1, e2);
                    }
                }
                prev_index = i;
            }
        }
    }
}

std::string RC(const std::string &edge_id) {
    if (edge_id[edge_id.size() - 1] == '-') {
        return edge_id.substr(0, edge_id.size() - 1) + '+';
    }
    return edge_id.substr(0, edge_id.size() - 1) + '-';
}

void SGraphBuilder::AddEdge(const std::string &prev_edge_id, const std::string &edge_id) {
    if (sgraph_.count(prev_edge_id) == 0) {
        sgraph_[prev_edge_id] = std::unordered_map<std::string, int>();
    }
    if (sgraph_[prev_edge_id].count(edge_id) == 0) {
        sgraph_[prev_edge_id][edge_id] = 0;
    }
    ++ sgraph_[prev_edge_id][edge_id];

    std::string prev_edge_id_rc = RC(prev_edge_id);
    std::string edge_id_rc = RC(edge_id);
    if (sgraph_.count(edge_id_rc) == 0) {
        sgraph_[edge_id_rc] = std::unordered_map<std::string, int>();
    }
    if (sgraph_[edge_id_rc].count(prev_edge_id_rc) == 0) {
        sgraph_[edge_id_rc][prev_edge_id_rc] = 0;
    }
    ++ sgraph_[edge_id_rc][prev_edge_id_rc];
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

const multigraph::Edge *GetEdgebyStr(const std::string &edge_id_str, const multigraph::MultiGraph &mg_) {
    int edge_id = atoi(edge_id_str.substr(0, edge_id_str.size() - 1).c_str());
    edge_id = edge_id_str[edge_id_str.size() - 1] == '-'? -edge_id: edge_id;
    return &mg_.edges.at(edge_id);
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
            if (cur_edge[cur_edge.size() - 1] == '-') {
                cur_graph_edge = cur_graph_edge->rc;
            }
            int v_len = cur_graph_edge->end->size();
            path_str += cur_graph_edge->getSeq().str().substr(0, cur_graph_edge->size() - v_len);
        }
        const std::string &cur_edge = subpath[subpath.size() - 1];
        const multigraph::Edge *cur_graph_edge = GetEdgebyStr(cur_edge, mg_);
        if (cur_edge[cur_edge.size() - 1] == '-') {
            cur_graph_edge = cur_graph_edge->rc;
        }
        int v_len = cur_graph_edge->start->size();
        path_str += cur_graph_edge->getSeq().str().substr(0, v_len);
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
        std::cerr << path_str.size() << " " << read_str.str().size() << " "
                 << start_pos << " " << end_pos << " "
                 << result.editDistance << " nuc1:" << prev_nuc
                 << " nuc2:" << post_nuc << " " << result.numLocations << std::endl;
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
}
