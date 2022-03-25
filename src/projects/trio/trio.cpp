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

using namespace trio;
using std::vector;
using std::pair;
using std::array;
using std::sort;
using logging::Logger;
using std::cout;
using std::cerr;

HaplotypeRemover::HaplotypeRemover(logging::Logger &logger, multigraph::MultiGraph &mg,
                                            const std::experimental::filesystem::path &haployak, const char haplotype,
                                            const std::experimental::filesystem::path &out_dir) : logger_(logger), mg(mg),
                                                                                                  haplotype_(haplotype), out_dir(out_dir) {

    auto bulges = getBulgeLabels();
    logger_.info() << "got " << bulges.size() << " bulges\n";
    haplo_map_type haplotypes;
    string s;
    std::ifstream haplo_file(haployak);
    while (std::getline(haplo_file, s)) {
        HaplotypeStats h(s);
        haplotypes[h.label] = h;
    }
    std::string out_name = "haplotype_";
    ensure_dir_existance(out_dir);
}

void HaplotypeRemover::deleteEdgeHaplo(int eid) {
    logger_.debug() << "Removing " << eid << endl;
    auto to_merge = mg.deleteEdgeById(eid);
    for (auto p: to_merge){
        HaplotypeStats new_haplo(haplotypes[p.second[0]]);
        new_haplo.label = p.first;
        haplotypes.insert(std::make_pair(p.first, new_haplo));
        for (size_t i = 1; i < p.second.size(); i++) {
            if (haplotypes[p.first].haplotype != haplotypes[p.second[i]].haplotype) {
                logger_.trace() << "Merging different haplotypes " << haplotypes[p.first].label <<
                " " << haplotypes[p.second[i]].label << endl;
            }
            haplotypes[p.first].appendKmerStats(haplotypes[p.second[i]]);
        }
    }
}

void HaplotypeRemover::cleanGraph() {
    bool changed = true;
    size_t tips = 0;
    size_t bulges = 0;
    while (changed) {
        changed = false;
        std::vector<int> eids;
        for (auto p: mg.edges) {
            eids.push_back(p.first);
        }
        for (auto eid: eids){
            if (mg.edges.find(eid) == mg.edges.end())
                continue;
            logger_.debug() << "considering " << eid << " label " << mg.edges[eid]->getLabel() << endl;
            if (mg.edges[eid]->isTip()) {
                logger_.debug() << "is being deleted as tip\n";
                if (mg.edges[eid]->size() < MAX_TIP_LENGTH) {
                    logger_.debug() << "is deleted as tip\n";
                    deleteEdgeHaplo(eid);
                    changed = true;
                    tips ++;
                }
            } else if (mg.edges[eid]->start->outDeg() == 2) {
                logger_.debug() << "is being deleted as bulge\n";
                auto first_e = mg.edges[eid]->start->outgoing[0];
                auto second_e = mg.edges[eid]->start->outgoing[1];
                if (first_e->end == second_e->end && first_e != second_e->rc
                && first_e->size() < BULGE_MULTIPLICATIVE_CUTOFF * second_e->size() &&  second_e->size() < BULGE_MULTIPLICATIVE_CUTOFF * first_e->size()) {
                    logger_.debug() << "is deleted as bulge\n";
                    char decision = AssignBulge(haplotypes[first_e->getLabel()], haplotypes[second_e->getLabel()]);
                    if (decision == haplotype_)
                        deleteEdgeHaplo(first_e->getId());
                    else
                        deleteEdgeHaplo(second_e->getId());
                    changed = true;
                    bulges ++;
                }
            }

        }
    }
    logger_.info() << "Deleted tips "<< tips << " Bulges " << bulges << endl;
}

std::unordered_map<std::string, std::string> HaplotypeRemover::getBulgeLabels() {
    std::set<int> used;
    std::unordered_map<std::string, std::string> res;
    for (auto p : mg.vertices) {
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

void HaplotypeRemover::updateFixableHaplotypes() {
    size_t count = 0;
    for (auto p: bulges) {
        std::string e_label = p.first;
        std::string alt_label = p.second;
        if (haplotypes.find(e_label) != haplotypes.end() && haplotypes[e_label].is_undefined()) {
            if (haplotypes[alt_label].haplotype == 'm') {
                haplotypes[e_label].haplotype = 'p';
                count++;
            }
            if (haplotypes[alt_label].haplotype == 'p') {
                haplotypes[e_label].haplotype = 'm';
                count++;
            }
        }
    }
    logger_.info() << "Updated " << count << " fixable bulges";
}

void HaplotypeRemover::updateAmbiguousHaplotypes() {
    size_t count = 0;
    for (auto p: bulges) {
        std::string e_label = p.first;
        std::string alt_label = p.second;
        if (haplotypes[e_label].is_undefined() && haplotypes[alt_label].is_undefined()) {
            auto top_h = haplotypes[e_label];
            auto bottom_h = haplotypes[alt_label];
            char decision = AssignBulge(top_h, bottom_h);
            if (decision != 'a') {
                haplotypes[e_label].haplotype = decision;
                haplotypes[alt_label].haplotype = other_haplo(decision);
                count ++;
            }
        }
    }
    logger_.info() << "Updated " << count << " ambiguous bulges";
}

void HaplotypeRemover::removeHaplotype() {
    size_t removed = 0;
    size_t bridges = 0;
    size_t removed_len = 0;
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<int> eids;
        for (auto p: mg.edges) {
            eids.push_back(p.first);
        }
        for (auto eid:eids){
            if (mg.edges.find(eid) == mg.edges.end())
                continue;
            auto label = mg.edges[eid]->getLabel();
            if (haplotypes.find(label) != haplotypes.end()) {
                if (haplotypes[label].haplotype == haplotype_) {
                    if (mg.edges[eid]->isSimpleBridge() && mg.edges[eid]->size() < SAVED_BRIDGE_CUTOFF) {
                        bridges ++;
                        logger_.info() << "Skipping edge " << eid << " as bridge\n";
                        continue;
                    }
                    removed_len += mg.edges[eid]->size();
                    deleteEdgeHaplo(eid);

                    logger_.trace() << "removing " << eid  << " label " << label << endl;
                    removed ++;
                    changed = true;
                } else { 
                    logger_.trace() << "skipping edge label" << label <<" "<< haplotypes[label].haplotype << endl;
                }
            } else {
                logger_.trace() << "skipping edge NOT FOUND label" << label << endl;
                
            }
        }
    }
    logger_.info() << "Saved " << bridges << "bridges\n";
    logger_.info() << "Removed " << removed << " edges of haplo " << haplotype_  << " total len " << removed_len << endl;
}

void HaplotypeRemover::process() {
    updateFixableHaplotypes();
    updateAmbiguousHaplotypes();
    removeHaplotype();
    logger_.debug() << "removed \n";
    mg.printEdgeGFA(out_dir / "before_clean.gfa", true);
    cleanGraph();
    mg.printEdgeGFA(out_dir / "mdbg.hpc.gfa", true);

}

//TODO:: хорошо бы чтоб все кроме этого не обращалось с файлами
std::experimental::filesystem::path trio::simplifyHaplo(logging::Logger &logger, size_t threads,
                                                  const std::experimental::filesystem::path &output_file,
                                                  const std::experimental::filesystem::path &diplo_graph,
                                                  const std::experimental::filesystem::path &haployak,
                                                  const char haplotype,  const std::experimental::filesystem::path &corrected_reads,
                                                  io::Library & reads,   const std::experimental::filesystem::path &dir) {
    multigraph::MultiGraph mmg;
    mmg.LoadGFA(diplo_graph, true);
//TODO:: it would be cool not to create twice
    multigraph::MultiGraph mg = mmg.DBG();
    std::string out_name = "haplotype_";
    out_name += other_haplo(haplotype);
    std::experimental::filesystem::path out_dir = dir / out_name;
    HaplotypeRemover hr(logger, mg, haployak, haplotype, out_dir);
    hr.process();
    mg.printEdgeGFA(output_file, true);

//printing alignments and contigs, should be refactored
    std::string out_aligns = out_name; out_aligns += ".alignments";
    std::string out_contigs = out_name; out_contigs += ".fasta";
    io::Library ref_lib;
    pipeline::LJAPipeline pipeline (ref_lib);
//TODO:: get rid of this magic const
    size_t k = 5001;
    std::vector<std::experimental::filesystem::path> uncompressed_results =
           pipeline.PolishingPhase(logger, threads, out_dir, out_dir, output_file,
                           corrected_reads, reads, StringContig::max_dimer_size / 2, k, false, true);


    return output_file;
}
