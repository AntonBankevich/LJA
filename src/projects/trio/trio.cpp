#include "trio.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <common/string_utils.hpp>
#include "dbg/dbg_construction.hpp"

#include <unordered_set>
#include <utility>
#include <vector>
#include <iostream>
#include <array>
#include <algorithm>

#include <lja/multi_graph.hpp>
#include <lja/subdataset_processing.hpp>

#include <lja/pipeline.hpp>

using std::vector;
using std::pair;
using std::array;
using std::sort;
using logging::Logger;

using std::cout;
using std::cerr;


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
            std::cout<< eid << endl;
            if (graph.edges[eid]->start->outDeg() == 2) {
                std::cout << "is being deleted as bulge\n";
                auto first_e = graph.edges[eid]->start->outgoing[0];
                auto second_e = graph.edges[eid]->start->outgoing[1];
                if (first_e->end == second_e->end && first_e != second_e->rc) {
                    std::cout << "is deleted as bulge\n";

                    graph.deleteEdgeById(first_e->getId());
                    changed = true;
                    bulges ++;
                }
            } else if (graph.edges[eid]->isTip()) {
                std::cout << "is being deleted as tip\n";
                if (graph.edges[eid]->size() < MAX_TIP_LENGTH) {
                    std::cout << "is deleted as tip\n";
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
                                                  const char haplotype,  const std::experimental::filesystem::path &corrected_reads, io::Library & reads,   const std::experimental::filesystem::path &dir) {
    multigraph::MultiGraph mmg;
    mmg.LoadGFA(diplo_graph, true);
    multigraph::MultiGraph mg = mmg.DBG();

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
    std::string out_name = "corrected_";
    size_t k = 5001;
    out_name+=other_haplo(haplotype);
    std::experimental::filesystem::path out_dir = dir / out_name;
    ensure_dir_existance(out_dir);

    updateFixableHaplotypes(bulges, haplotypes);
    updateAmbiguousHaplotypes(bulges, haplotypes, mg);
    for (auto p: haplotypes)
        logger.debug() << p.second.label << " " << p.second.haplotype << endl;
    removeHaplotype(haplotypes, mg, haplotype, logger);
    cout << "removed \n";
    mg.printEdgeGFA(out_dir / "before_clean.gfa", true);
    cleanGraph(mg);
    cout << "cleaned \n";
    mg.printEdgeGFA(out_dir / "after_clean.gfa", true);
    mg.printEdgeGFA(output_file, true);

//printing alignments and contigs, should be refactored
    std::string out_aligns = out_name; out_aligns += ".alignments";
    std::string out_contigs = out_name; out_contigs += ".fasta";
    io::Library ref_lib;
    multigraph::LJAPipeline pipeline (ref_lib);
    std::vector<std::experimental::filesystem::path> uncompressed_results =
           pipeline.PolishingPhase(logger, threads, out_dir, out_dir, output_file,
                           corrected_reads, reads, StringContig::max_dimer_size / 2, k, false, true);
/*
    hashing::RollingHash hasher(k, 239);
    SparseDBG dbg = DBGPipeline(logger, hasher, w, reads, dir, threads);
    dbg.fillAnchors(w, logger, threads);
    size_t extension_size = std::max<size_t>(k * 5 / 2, 3000);
    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, false);
    io::SeqReader reader(reads);
    readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);

    std::vector<Contig> contigs = mmg.getCutEdges();
    PrintAlignments(logger, threads, contigs, readStorage, k, dir / out_aligns );
    readStorage.printFasta(logger, dir / out_contigs);
    */
    return output_file;
}
