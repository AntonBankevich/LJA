#include "trio.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <common/string_utils.hpp>

#include <unordered_set>
#include <utility>
#include <vector>
#include <iostream>
#include <array>
#include <algorithm>
#include <lja/multi_graph.hpp>
using std::vector;
using std::pair;
using std::array;
using std::sort;
using logging::Logger;

using std::cout;
using std::cerr;

struct HaplotypeStats {
    char haplotype;
    std::string label;
//order: p, m, as in file;
    std::array<int, 2> decisive_strips;
    std::array<int, 2> decisive_counts;
    int total_kmers;

//        32         m       0       273     28      390     22      24      112906  17
//#s->seq[i].name, type, s->cnt[i].sc[0], s->cnt[i].sc[1],c[0<<2|2], c[2<<2|0], c[0<<2|1], c[1<<2|0], s->cnt[i].nk, c[0])
//TODO: arr 6, 7, 9
    HaplotypeStats(std::string s) {
        std::vector<std::string> tokens = ::split(s);
        haplotype = tokens[1][0];
        label = tokens[0];
        decisive_strips = std::array<int,2>{stoi(tokens[2]), stoi(tokens[3])};
        decisive_counts = std::array<int,2>{stoi(tokens[4]), stoi(tokens[5])};
        total_kmers = stoi(tokens[8]);
    }
    HaplotypeStats():haplotype('u'),label(""), decisive_counts{0,0}, decisive_strips{0,0}, total_kmers(0) {}
    bool is_undefined() {
        return (haplotype != 'm' and haplotype != 'p');
    }
};

typedef  std::unordered_map<std::string, HaplotypeStats> haplo_map_type;
//primitive clean graph

void cleanGraph(multigraph::MultiGraph &graph) {
    bool changed = true;
    size_t MAX_TIP_LENGTH = 1000000;
    size_t tips = 0;
    size_t bulges = 0;
    while (changed) {
        changed = false;
        std::vector<int> eids;
        for (auto p: graph.edges) {
            eids.push_back(p.first);
        }
        for (auto eid: eids){
            if (graph.edges.find(eid) == graph.edges.end())
                continue;
            if (graph.edges[eid]->start->outDeg() == 2) {
                auto first_e = graph.edges[eid]->start->outgoing[0];
                auto second_e = graph.edges[eid]->start->outgoing[1];
                if (first_e->end == second_e->end && first_e != second_e->rc) {
                    graph.deleteEdgeById(first_e->getId());
                    changed = true;
                    bulges ++;
                }
            } else if (graph.edges[eid]->isTip()) {
                if (graph.edges[eid]->size() < MAX_TIP_LENGTH) {
                    graph.deleteEdgeById(eid);
                    changed = true;
                    tips ++;
                }
            }
        }
    }
    std::cout << "Deleted tips "<< tips << " Bulges " << bulges << endl;
}

std::unordered_map<std::string, std::string> getBulgeLabels(multigraph::MultiGraph &graph) {
    std::set<int> used;
    std::unordered_map<std::string, std::string> res;
    for (auto p : graph.vertices) {
        int vid = p.first;
        multigraph::Vertex* v = p.second;
        if (v->outDeg() == 2) {
            multigraph::Edge *e1 = v->outgoing[0];
            multigraph::Edge *e2 = v->outgoing[1];
            if (e1->end == e2->end) {
                if (used.find(e1->getId()) == used.end() && used.find(e2->getId()) == used.end()) {
                    used.insert(e1->getId());
                    used.insert(e2->getId());
                    res[e1->getLabel()] = e2->getLabel();
                    res[e2->getLabel()] = e1->getLabel();
                }
            }
        }
    }
    return res;
}

void updateFixableHaplotypes(std::unordered_map<std::string, std::string> & bulges, haplo_map_type &haplostats) {
    size_t count = 0;
    for (auto p: bulges) {
        std::string e_label = p.first;
        std::string alt_label = p.second;
//        cerr << e_label << " " << alt_label << " " << haplostats[e_label].haplotype << " " << haplostats[alt_label].haplotype << endl;
        if (haplostats.find(e_label) != haplostats.end() && haplostats[e_label].is_undefined()) {
            if (haplostats[alt_label].haplotype == 'm') {
                haplostats[e_label].haplotype = 'p';
                count++;
            }
            if (haplostats[alt_label].haplotype == 'p') {
                haplostats[e_label].haplotype = 'm';
                count++;
            }
        }
    }
    std::cout << "Updated " << count << " fixable bulges";
}

void updateAmbiguousHaplotypes(std::unordered_map<std::string, std::string> & bulges, haplo_map_type &haplotypes, multigraph::MultiGraph &graph) {
    size_t count = 0;
//size_t short_length = 100
//Currently removed
    for (auto p: bulges) {
        std::string e_label = p.first;
        std::string alt_label = p.second;
        if (haplotypes[e_label].is_undefined() && haplotypes[alt_label].is_undefined()) {
            auto top_h = haplotypes[e_label];
            auto bottom_h = haplotypes[alt_label];
            size_t pat_count = top_h.decisive_counts[0] * bottom_h.decisive_counts[1];
            size_t mat_count = top_h.decisive_counts[1] * bottom_h.decisive_counts[0];
            size_t top_total = top_h.decisive_counts[0] + top_h.decisive_counts[1];
            size_t bottom_total = bottom_h.decisive_counts[0] + bottom_h.decisive_counts[1];
            char decision = 'a';
            if (mat_count > pat_count || (top_total == 0 && bottom_h.decisive_counts[0] >  bottom_h.decisive_counts[1])
                || (bottom_total == 0 && top_h.decisive_counts[1] > top_h.decisive_counts[0]))
                decision = 'm';
            else if (mat_count < pat_count || (top_total == 0 && bottom_h.decisive_counts[1] >  bottom_h.decisive_counts[0])
                     || (bottom_total == 0 && top_h.decisive_counts[0] > top_h.decisive_counts[1]))
                decision = 'p';
            if (decision == 'm') {
                haplotypes[e_label].haplotype = 'm';
                haplotypes[alt_label].haplotype = 'p';
                count ++;
            }
            if (decision == 'p') {
                haplotypes[e_label].haplotype = 'p';
                haplotypes[alt_label].haplotype = 'm';
                count ++;
            }
        }
    }
    std::cout << "Updated " << count << " ambiguous bulges";
}

void removeHaplotype(haplo_map_type &haplotypes, multigraph::MultiGraph &graph, char haplo_to_remove, logging::Logger &logger) {
    size_t removed = 0;
    size_t bridges = 0;
    size_t removed_len = 0;
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<int> eids;
        for (auto p: graph.edges) {
            eids.push_back(p.first);
        }
        for (auto eid:eids){
            if (graph.edges.find(eid) == graph.edges.end())
                continue;
            auto label = graph.edges[eid]->getLabel();
            if (haplotypes.find(label) != haplotypes.end()) {
                if (haplotypes[label].haplotype == haplo_to_remove) {
                    if (graph.edges[eid]->isBridge()) {
                        bridges ++;
                        continue;
                    }
                    removed_len += graph.edges[eid]->size();
                    graph.deleteEdgeById(eid);
                    logger.trace() << "removing " << eid << endl;
                    removed ++;
                    changed = true;
                }
            }
        }
    }
    logger.info() << "Saved " << bridges << "bridges\n";
    logger.info() << "Removed " << removed << " edges of haplo " << haplo_to_remove  << " total len " << removed_len << endl;
}




std::experimental::filesystem::path simplifyHaplo(logging::Logger &logger, size_t threads,
                                                  const std::experimental::filesystem::path &output_file,
                                                  const std::experimental::filesystem::path &diplo_graph,
                                                  const std::experimental::filesystem::path &haployak,
                                                  const char haplotype) {
    multigraph::MultiGraph mmg;
    mmg.LoadGFA(diplo_graph);
    multigraph::MultiGraph mg = mmg.DBG();
    mg.printEdgeGFA("tsat.gfa",true);

    auto bulges = getBulgeLabels(mg);
    logger.info() << "got "<< bulges.size() << " bulges\n";
//to_separate_fucntion
    haplo_map_type haplotypes;
    string s;
    std::ifstream haplo_file(haployak);
    while (std::getline(haplo_file, s)){
        HaplotypeStats h(s);
        haplotypes[h.label] = h;
    }

    updateFixableHaplotypes(bulges, haplotypes);
    cout << "updated simple\n";
    updateAmbiguousHaplotypes(bulges, haplotypes, mg);
    cout << "updated complex\n";
    for (auto p: haplotypes)
        logger.debug() << p.second.label << " " << p.second.haplotype << endl;
    removeHaplotype(haplotypes, mg, haplotype, logger);
    cout << "removed \n";
    mg.printEdgeGFA("before_clean.gfa", true);
    cleanGraph(mg);
    cout << "cleaned \n";
    mg.printEdgeGFA("after_clean.gfa", true);
    mg.printEdgeGFA(output_file, true);
    return output_file;
}
