//
// Created by Tatiana Dvorkina on 21.06.2022.
//

#include "edlib/edlib.h"

#include "VertexExtender.h"

using namespace nano;

bool VertexExtender::IsNearlyDimerString(const std::string &str) {
    std::string compressed_str;
    compressed_str.push_back(str[0]);
    compressed_str.push_back(str[1]);
    std::unordered_map<std::string, int> dimers;
    dimers[compressed_str] = 1;
    size_t i = 2;
    while (i - 1 < str.size()) {
        if (! (compressed_str[compressed_str.size() - 2] == str[i] &&
               compressed_str[compressed_str.size() - 1] == str[i + 1])) {
            compressed_str.push_back(str[i]);
            ++ i;
        } else {
            std::string dimer;
            dimer.push_back(str[i]);
            i += 2;
        }
    }
    return compressed_str.size()*2 < str.size();
}

bool VertexExtender::IsDimerString(const std::string &str, int sz) {
    for (int i = 0; i < sz; ++ i){
        if (str[2*i] != str[0] || str[2*i + 1] != str[1]) {
            return false;
        }
    }
    return true;
}

int VertexExtender::IsDimerStructure(const int v_id, multigraph::MultiGraph &mg){
    const int PREFIX_SIZE = 100;
    const int DIMER_LEN = 2;
    const int SUFFIX_LEN = 20;

    multigraph::Edge *e1 = mg.vertices.at(v_id).outgoing[0];
    multigraph::Edge *e2 = mg.vertices.at(v_id).outgoing[1];
    int v_len = mg.vertices.at(v_id).size();
    std::cerr << "Consider " << e1->getId() << " " << e2->getId() << " " << v_len << std::endl;
    std::string e1_str = e1->getSeq().str().substr(v_len,
                                                   std::min(PREFIX_SIZE,
                                                            (int) (e1->getSeq().size() - v_len) ));
    std::string e2_str = e2->getSeq().str().substr(v_len,
                                                   std::min(PREFIX_SIZE,
                                                            (int) (e2->getSeq().size() - v_len) ));
    std::cerr << e1_str << std::endl;
    std::cerr << e2_str << std::endl;
    EdlibAlignResult result = edlibAlign(e1_str.c_str(), e1_str.size(),
                                         e2_str.c_str(), e2_str.size(),
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        std::cerr << result.alignmentLength << " " << int(result.alignment[0]) << " " << int(result.alignment[1]) << std::endl;
        std::cerr << result.editDistance << std::endl;
        std::string vertex_end = e1->getSeq().str().substr(v_len - SUFFIX_LEN, SUFFIX_LEN);
        std::reverse(vertex_end.begin(), vertex_end.end());
        std::cerr << vertex_end << std::endl;
        if (result.alignment[0] == result.alignment[1] && int(result.alignment[0]) == 1) {
            std::cerr << vertex_end[1] << vertex_end[0] << " " << e1_str[0] << e1_str[1] << std::endl;
            std::cerr << IsDimerString(vertex_end, DIMER_LEN) << std::endl;
            if (vertex_end[1] == e1_str[0] &&
                vertex_end[0] == e1_str[1] && IsDimerString(vertex_end, DIMER_LEN)) {
                edlibFreeAlignResult(result);
                return 1;
            }
        } else if (result.alignment[0] == result.alignment[1] && int(result.alignment[0]) == 2) {
            std::cerr << vertex_end[1] << vertex_end[0] << " " << e2_str[0] << e2_str[1] << std::endl;
            std::cerr << IsDimerString(vertex_end, DIMER_LEN) << std::endl;
            if (vertex_end[1] == e2_str[0] &&
                vertex_end[0] == e2_str[1] && IsDimerString(vertex_end, DIMER_LEN)) {
                edlibFreeAlignResult(result);
                return 0;
            }
        }
    }
    edlibFreeAlignResult(result);
    return -1;
}

std::unordered_map<int, int> VertexExtender::ExtendVertices(multigraph::MultiGraph &mg){
    std::unordered_map<int, int> vertices_to_extend;
    for (auto const&[key, val]: mg.vertices) {
        if (mg.vertices.at(key).outgoing.size() == 2) {
            int isdimer = IsDimerStructure(key, mg);
            std::cerr << isdimer << std::endl;
            if (isdimer > -1) vertices_to_extend.insert(std::pair<int, int> (key, isdimer));
        }
    }
    std::cerr << "Vertex extend num: " << vertices_to_extend.size() << std::endl;

    std::unordered_map<int, int> new_unique;
    std::unordered_map<int, int> new_unique_op;
    for (auto v_id: vertices_to_extend) {
        int index = v_id.second == 2? 0: 1;

        int dimer_sz = 0;
        std::string vertex_str = mg.vertices.at(v_id.first).rc->seq.str();
        int i = 0;
        while (i + 3 < vertex_str.size() &&
               vertex_str[i] == vertex_str[i + 2] && vertex_str[i + 1] == vertex_str[i + 3]) {
            i += 2;
            ++ dimer_sz;
        }
        std::cerr << "Dimer:" << mg.vertices.at(v_id.first).outgoing[0]->getId()
                  << " " << mg.vertices.at(v_id.first).outgoing[1]->getId() << " " << dimer_sz << std::endl;


        const multigraph::Edge *e = mg.vertices.at(v_id.first).outgoing[index];
        int v_len = mg.vertices.at(v_id.first).size();
        std::string dimer;
        dimer.push_back(e->getSeq().str()[v_len - 2]);
        dimer.push_back(e->getSeq().str()[v_len - 1]);
        std::cerr << v_id.first << " " << v_id.second << " " << e->getId() << " " << e->size()
                  << " " << dimer << std::endl;
        Sequence new_edge_seq = Sequence(e->getSeq().str().substr(0, v_len)
                                         + dimer
                                         + e->getSeq().str().substr(v_len, e->getSeq().str().size() - v_len));
        //std::cerr << v_id.first << " " << v_id.second << " " << e->getId() << " " << new_edge_seq.size() << std::endl;
        int new_edge_id = mg.addEdge(*(e->start), *(e->end), new_edge_seq).getId();
        //std::cerr << e->getId() << " " << new_edge_id << " " << mg.edges.at(new_edge_id).size() << std::endl;
        if (new_unique.count(e->getId()) == 0 && new_unique_op.count(e->getId()) == 0) {
            new_unique[e->getId()] = new_edge_id;
            new_unique_op[new_edge_id] = e->getId();
            new_unique[e->rc->getId()] = mg.edges.at(new_edge_id).rc->getId();
            new_unique_op[mg.edges.at(new_edge_id).rc->getId()] = e->rc->getId();
        } else if (new_unique.count(e->getId()) == 1) {
            new_unique[e->getId()] = new_edge_id;
            new_unique_op[new_edge_id] = e->getId();
            new_unique[e->rc->getId()] = mg.edges.at(new_edge_id).rc->getId();
            new_unique_op[mg.edges.at(new_edge_id).rc->getId()] = e->rc->getId();
        } else if (new_unique_op.count(e->getId()) == 1) {
            new_unique[new_unique_op[e->getId()]] = new_edge_id;
            new_unique_op[new_edge_id] = new_unique_op[e->getId()];

            new_unique[new_unique_op[e->getId()]*(-1)] = new_edge_id*(-1);
            new_unique_op[new_edge_id*(-1)] = new_unique_op[e->getId()]*(-1);
        }
        mg.internalRemoveEdge(e->getId());
    }

    for (auto v_id: vertices_to_extend) {
        if (mg.vertices.count(v_id.first) > 0 && mg.vertices.at(v_id.first).outgoing.size() == 2) {
            std::unordered_map<int, int> cur_new_unique = mg.extendVertex(v_id.first);
            for (auto ue: cur_new_unique) {
                std::cerr << "After extension " << ue.first << " " << ue.second
                          << " " << mg.edges.at(ue.second).size() << std::endl;
                if (new_unique.count(ue.first) == 0 && new_unique_op.count(ue.first) == 0) {
                    new_unique[ue.first] = ue.second;
                    new_unique_op[ue.second] = ue.first;
                } else if (new_unique.count(ue.first) == 1) {
                    if (new_unique[ue.first] != -1) {
                        new_unique[ue.first] = ue.second;
                        new_unique_op[ue.second] = ue.first;
                    }
                } else if (new_unique_op.count(ue.first) == 1) {
                    new_unique[new_unique_op[ue.first]] = ue.second;
                    new_unique_op[ue.second] = new_unique_op[ue.first];
                }
            }
        }
    }
    std::cerr << std::endl;
    for (const auto&[key, val]: new_unique) {
        std::cerr << key << " " << val << std::endl;
    }

    for (auto v_id: vertices_to_extend) {
        if (mg.vertices.count(v_id.first) > 0) {
            mg.compressVertex(v_id.first);
        }
    }
    return new_unique;
}