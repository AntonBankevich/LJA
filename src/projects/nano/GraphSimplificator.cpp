//
// Created by Tatiana Dvorkina on 04.05.2022.
//

#include "GraphSimplificator.h"

using namespace nano;

void GraphSimplificator::Simplify(multigraph::MultiGraph &mg, const std::experimental::filesystem::path &dir){
    std::unordered_set<int> to_delete;
    std::ofstream os_cut;
    const std::experimental::filesystem::path &final_fasta = dir / "final_test.fasta";
    os_cut.open(final_fasta);
    for (auto const&[e1_id, e2_id]: chains_) {
        if (inside_edges_.count(e1_id) == 1) continue;
        std::cerr << "Cur chain " << e1_id << std::endl;
        std::vector<int> cur_chain;
        cur_chain.push_back(e1_id);
        const multigraph::Edge *cur_edge = &mg.edges.at(e1_id);
        bool addEdge = true;
        if (to_delete.count(cur_edge->rc->getId()) == 1) addEdge = false;
        to_delete.insert(cur_edge->getId());
        Sequence new_edge_seq = cur_edge->getSeq();
        int cur_len = new_edge_seq.size();
        int cur_e_id = e1_id;
        while (chains_.count(cur_e_id) > 0) {
            std::vector<int> nonunique_subpath = std::vector<int>();
            if (sgraph_.count(cur_e_id) > 0 && sgraph_.at(cur_e_id).count(chains_[cur_e_id]) > 0) {
                nonunique_subpath = sgraph_.at(cur_e_id).at(chains_[cur_e_id]).second;
            }
            for (auto subpath_edge_id: nonunique_subpath) {
                cur_chain.push_back(subpath_edge_id);
                const multigraph::Edge *subpath_edge = &mg.edges.at(subpath_edge_id);
                to_delete.insert(subpath_edge->getId());
                int v_len = subpath_edge->start->size();
                new_edge_seq = new_edge_seq +
                                Sequence(subpath_edge->getSeq().str().substr(v_len, subpath_edge->size() - v_len));
                cur_len += subpath_edge->size() - v_len;
            }
            cur_chain.push_back(chains_[cur_e_id]);
            cur_e_id = chains_[cur_e_id];
            const multigraph::Edge *cur_edge = &mg.edges.at(cur_e_id);
            int v_len = cur_edge->start->size();
            new_edge_seq = new_edge_seq + Sequence(cur_edge->getSeq().str().
                                                      substr(v_len, cur_edge->size() - v_len));
            cur_len += cur_edge->size() - v_len;
            to_delete.insert(cur_edge->getId());
        }
        std::string to_print = "";
        for (auto it: cur_chain) {
            to_print += std::to_string(it) + " ";
        }
        std::cerr << "New edge: " << to_print << std::endl;
        if (addEdge) {
            const multigraph::Edge *new_edge = &mg.addEdge(*mg.edges.at(cur_chain[0]).start,
                       *mg.edges.at(cur_chain[cur_chain.size() - 1]).end, new_edge_seq);
            os_cut << ">" << new_edge->getId() << "\n" << new_edge_seq << "\n";
            std::cerr << "New edge id: " << new_edge->getId() << " " << new_edge->size() << std::endl;
        }
    }
    os_cut.close();
    std::unordered_set<int> compress_vertices;
    for (auto it: to_delete) {
        if (mg.edges.count(it) > 0) {
            std::cerr << it << " ";
            compress_vertices.insert(mg.edges.at(it).start->id);
            compress_vertices.insert(mg.edges.at(it).end->id);
            mg.internalRemoveEdge(it);
        }
    }
    for (auto it: compress_vertices) {
        if (mg.vertices.count(it) > 0) {
            mg.compressVertex(it);
        }
    }
    std::cerr << std::endl;
}

void GraphSimplificator::ResolveWithMajor(multigraph::MultiGraph &mg) {
    std::vector<std::unordered_set<int>> in_edges;
    std::vector<std::unordered_set<int>> out_edges;
//    std::unordered_map<int, std::unordered_set<int>> parallel_edges;
//    for (auto const &[e1_id, val1]: sgraph_) {
//        parallel_edges.insert(std::pair<int, std::unordered_set<int>>(e1_id, std::unordered_set<int>()));
//        std::cerr << "Cur edges " << e1_id << std::endl;
//        for (auto const &[e2_id, num2]: val1) {
//            const multigraph::Edge *edge2 = &mg.edges.at(e2_id);
//            for (auto const &[e3_id, num3]: sgraph_.at(edge2->rc->getId())) {
//                const multigraph::Edge *edge = &mg.edges.at(e3_id);
//                std::cerr << " Insert " << edge->rc->getId() << std::endl;
//                parallel_edges[e1_id].insert(edge->rc->getId());
//            }
//        }
//    }
//    for (auto const &[e1_id, val1]: sgraph_) {
//        in_edges.push_back(std::unordered_set<int>());
//        std::cerr << "e1 " << e1_id << " " << parallel_edges.count(e1_id) << std::endl;
//        for (auto const e_id: parallel_edges[e1_id]) {
//            const multigraph::Edge *edge1 = &mg.edges.at(e_id);
//            for (auto const edge: edge1->end->rc->outgoing) {
//                in_edges[in_edges.size() - 1].insert(edge->rc->getId());
//            }
//        }
//        out_edges.push_back(std::unordered_set<int>());
//        for (auto const &[e2_id, num]: val1){
//            int e2_id_rc = (mg.edges.at(e2_id)).rc->getId();
//            std::cerr << " e2 " << e2_id << " " << parallel_edges.count(e2_id_rc) << std::endl;
//            for (auto const e_id: parallel_edges[e2_id_rc]) {
//                const multigraph::Edge *edge2 = &mg.edges.at(e_id);
//                for (auto const edge: edge2->rc->start->outgoing) {
//                    std::cerr << "  e3 " << edge->getId() << std::endl;
//                    out_edges[out_edges.size() - 1].insert(edge->getId());
//                }
//            }
//        }
//    }
    std::cerr << "Start resolving" << std::endl;
    for (auto const &[key, val]: mg.vertices){
        if (val.outgoing.size() == 2 && val.rc->outgoing.size() == 2) {
            bool is_uedges = true;
            for (auto const e: val.outgoing)
                if (sgraph_.count(e->getId()) == 0) {
                    is_uedges = false;
                    break;
                }
            for (auto const e: val.rc->outgoing)
                if (sgraph_.count(e->getId()) == 0) {
                    is_uedges = false;
                    break;
                }
            if (!is_uedges) continue;
            in_edges.push_back(std::unordered_set<int>());
            out_edges.push_back(std::unordered_set<int>());
            for (auto const e: val.outgoing) out_edges[out_edges.size() - 1].insert(e->getId());
            for (auto const e: val.rc->outgoing) in_edges[in_edges.size() - 1].insert(e->rc->getId());
        }
    }

    for (int i = 0; i < in_edges.size(); ++ i) {
        std::unordered_set<int> &cur_in_edges = in_edges[i];
        std::unordered_set<int> &cur_out_edges = out_edges[i];
        if (cur_in_edges.size() == 2 && cur_out_edges.size() == 2) {
            std::vector<std::pair<std::pair<int, int>, int>> transitions;
            for (int in_edge_id: cur_in_edges) {
                for (int out_edge_id: cur_out_edges) {
                    std::pair<int, int> edge_pair(in_edge_id, out_edge_id);
                    if (sgraph_.count(in_edge_id) > 0 && sgraph_.at(in_edge_id).count(out_edge_id) > 0) {
                        transitions.push_back(
                                std::pair<std::pair<int, int>, int>(edge_pair,
                                                                    sgraph_.at(in_edge_id).at(out_edge_id).first));
                    } else {
                        transitions.push_back(std::pair<std::pair<int, int>, int>(edge_pair, 0));
                    }
                }
            }
            std::sort(transitions.begin(), transitions.end(),
                      [](std::pair<std::pair<int, int>, int> a, std::pair<std::pair<int, int>, int> b) {
                          return a.second > b.second;
                      });
            std::cerr << "Transition" << std::endl;
            for (auto trans: transitions) {
                std::cerr << trans.first.first << " " << trans.first.second << " " << trans.second << std::endl;
            }
            if (transitions[0].second == 0) continue;

            int j = 1;
            while (j < transitions.size() &&
                  (transitions[0].first.first == transitions[j].first.first
                  || transitions[0].first.second == transitions[j].first.second)) ++j;

            if (j < transitions.size()) {
                std::pair<int, int> best_transition = transitions[0].first;
                chains_[best_transition.first] = best_transition.second;
                std::cerr << "First " << best_transition.first << "->" << best_transition.second << std::endl;
                inside_edges_.insert(best_transition.second);
                chains_[transitions[j].first.first] = transitions[j].first.second;
                inside_edges_.insert(transitions[j].first.second);
                std::cerr << "Second " << transitions[j].first.first << "->" << transitions[j].first.second
                          << std::endl;
            }
        }
    }
    for (auto it: chains_) {
        std::cerr << it.first << " " << it.second << std::endl;
    }
}