#include "trio.hpp"
#include <dbg/multi_graph.hpp>
#include "dbg/assembly_graph.hpp"
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <unordered_set>
#include <utility>
#include <vector>
#include <iostream>
#include <array>
#include <algorithm>
#include <dbg/graph_algorithms.hpp>

using namespace trio;
using std::vector;
using std::pair;
using std::array;
using std::sort;
using logging::Logger;
using std::cout;
using std::cerr;
using namespace multigraph;

HaplotypeRemover::HaplotypeRemover(logging::Logger &logger, size_t threads, multigraph::MultiGraph &mg,
                                            const std::experimental::filesystem::path &haployak, const Haplotype haplotype,
                                            const std::experimental::filesystem::path &out_dir, const size_t saved_bridge_cutoff)
                                            : logger_(logger), mg(mg), haplotype_(haplotype), out_dir(out_dir),
                                            saved_bridge_cutoff(saved_bridge_cutoff) {

    string s;
    std::ifstream haplo_file(haployak);
    while (std::getline(haplo_file, s)) {
        HaplotypeStats h(s);
        haplotype_info[h.label.front()] = h;
    }
//    std::string out_name = "haplotype_";
    ensure_dir_existance(out_dir);
}

void HaplotypeRemover::deleteEdgeHaplo(EdgeId eid) {
    logger_.debug() << "Removing " << eid <<std::endl;
    eid->getStart().removeEdge(*eid);
//    mg.internalRemoveEdge(*eid);
}


void HaplotypeRemover::compressAllVertices() {
    size_t all_count = 0;
    std::vector<GraphPath> paths = ag::AllUnbranchingPaths<MGTraits>(logger_, threads, mg);
    for (GraphPath &path: paths) {
        std::vector<Edge::id_type> ids;
        std::vector<Edge::id_type> rcids;
        for(Edge &edge : path.edges()) {
            ids.emplace_back(edge.getInnerId());
        }
        GraphPath rc = path.RC();
        for(Edge &edge : rc.edges()) {
            ids.emplace_back(edge.getInnerId());
        }
        HaplotypeStats new_haplo(haplotype_info[path.getEdge(0).getInnerId()]);
        for(size_t i = 1; i < path.size(); i++) {
            Edge::id_type eid = path.getEdge(i).getInnerId();
            if (new_haplo.haplotype != haplotype_info[eid].haplotype) {
                logger_.trace() << "Merging different haplotypes " << path.str() <<
                                " " << i << " " << new_haplo.haplotype << " " << haplotype_info[eid].haplotype <<std::endl;
            }
            new_haplo.appendKmerStats(haplotype_info[eid]);
        }
        Edge &merged = ag::CompressPath(path);
        new_haplo.label = std::move(ids);
        haplotype_info[merged.getInnerId()] = new_haplo;
        new_haplo.label = std::move(rcids);
        haplotype_info[merged.rc().getInnerId()] = new_haplo;
    }
    logger_.info() << "Compressed " << paths.size() << " unbranching paths" <<std::endl;
}

void HaplotypeRemover::cleanGraph() {
    bool changed = true;
    size_t tips = 0;
    size_t bulges_cnt = 0;
    size_t iter = 0;
    while (changed) {
        changed = false;
        logger_.info() << "Iteration " << ++iter << " of graph cleaning." <<std::endl;
        std::unordered_set<EdgeId> to_delete;
        for (Edge &edge : mg.edgesUnique()) {
            logger_.debug() << "considering " << edge.getId() << " label " << edge.getLabel() <<std::endl;
            if (edge.getFinish().outDeg() == 0 || edge.getStart().inDeg() == 0) {
                logger_.debug() << "is being deleted as tip\n";
                if (edge.fullSize() < MAX_TIP_LENGTH) {
                    logger_.debug() << "is deleted as tip\n";
                    to_delete.insert(edge.getId());
                    to_delete.insert(edge.rc().getId());
                    changed = true;
                    tips ++;
                }
            } else if (edge.getStart().outDeg() == 2) {
                logger_.debug() << "is being deleted as bulge\n";
                MGEdge & first_e = edge.getStart().front();
                MGEdge & second_e = edge.getStart().back();
                if (first_e.getFinish() == second_e.getFinish() && first_e != second_e.rc()
                    && first_e.fullSize() < BULGE_MULTIPLICATIVE_CUTOFF * second_e.fullSize() && second_e.fullSize() < BULGE_MULTIPLICATIVE_CUTOFF *
                                                                                                                       first_e.fullSize()) {
                    logger_.debug() << "is deleted as bulge\n";
                    Haplotype decision = AssignBulge(haplotype_info[first_e.getInnerId()], haplotype_info[second_e.getInnerId()]);
                    if (decision == haplotype_) {
                        to_delete.insert(first_e.getId());
                        to_delete.insert(first_e.rc().getId());
                    }
                    else {
                        to_delete.insert(second_e.getId());
                        to_delete.insert(second_e.rc().getId());
                    }
                    changed = true;
                    bulges_cnt ++;
                }
            }
        }
        for(EdgeId eid : to_delete) {
            if(eid->isCanonical())
                deleteEdgeHaplo(eid);
        }
        compressAllVertices();
    }
    logger_.info() << "Deleted tips " << tips << " Bulges " << bulges_cnt << std::endl;
}

std::vector<std::pair<EdgeId, EdgeId>> HaplotypeRemover::getBulgeLabels() {
    std::vector<std::pair<EdgeId, EdgeId>> res;
    for (MGVertex &v : mg.vertices()) {
        if (v.outDeg() == 2) {
            multigraph::MGEdge &e1 = v.front();
            multigraph::MGEdge &e2 = v.back();
            MGVertex &f = e1.getFinish().rc();
            if (e1.getFinish() == e2.getFinish() && (v <= f || f.outDeg() > 2)) {
                res.emplace_back(e1.getId(), e2.getId());
            }
        }
    }
    return res;
}

void HaplotypeRemover::updateFixableHaplotypes(const std::vector<std::pair<EdgeId, EdgeId>> &bulges) {
    size_t count = 0;
    for (const auto &p: bulges) {
        Edge &e1 = *p.first;
        Edge &e2 = *p.second;
        VERIFY(haplotype_info.find(e1.getInnerId()) != haplotype_info.end());
        VERIFY(haplotype_info.find(e2.getInnerId()) != haplotype_info.end());
        Haplotype &h1 = haplotype_info[e1.getInnerId()].haplotype;
        Haplotype &h2 = haplotype_info[e2.getInnerId()].haplotype;
        if(is_defined_haplo(h1) && !is_defined_haplo(h2)) {
            h2 = other_haplo(h1);
            count++;
        }
        if(is_defined_haplo(h2) && !is_defined_haplo(h1)) {
            h1 = other_haplo(h2);
            count++;
        }
    }
    logger_.info() << "Updated " << count << " fixable bulges";
}

void HaplotypeRemover::updateAmbiguousHaplotypes(const std::vector<std::pair<EdgeId, EdgeId>> &bulges) {
    size_t count = 0;
    for (const auto &p: bulges) {
        Edge &e1 = *p.first;
        Edge &e2 = *p.second;
        VERIFY(haplotype_info.find(e1.getInnerId()) != haplotype_info.end());
        VERIFY(haplotype_info.find(e2.getInnerId()) != haplotype_info.end());
        HaplotypeStats &top = haplotype_info[e1.getInnerId()];
        HaplotypeStats &bottom = haplotype_info[e2.getInnerId()];
        if (top.is_undefined() && bottom.is_undefined()) {
            Haplotype decision = AssignBulge(top, bottom);
            if (decision != Haplotype::Shared) {
                top.haplotype = decision;
                bottom.haplotype = other_haplo(decision);
                count ++;
            }
        }
    }
    logger_.info() << "Updated " << count << " ambiguous bulges";
}

bool isTip(const Edge &edge) {
    return (edge.getStart().inDeg() == 0  || edge.getFinish().outDeg() == 0);
}


bool isSimpleBridge(const Edge &edge) {
    if (isTip(edge))
        return false;
    for (MGEdge &alt_e: edge.getStart()) {
        if (alt_e != edge && !isTip(alt_e))
            return false;
    }
    for (MGEdge &alt_e: edge.rc().getStart()) {
        if (alt_e != edge.rc() && !isTip(alt_e))
            return false;
    }
    return true;
}

void HaplotypeRemover::removeHaplotype() {
    size_t removed = 0;
    size_t bridges = 0;
    size_t removed_len = 0;
    bool changed = true;
    while (changed) {
        std::unordered_set<EdgeId> to_delete;
        changed = false;
        for (Edge &edge : mg.edgesUnique()){
            VERIFY(haplotype_info.find(edge.getInnerId()) != haplotype_info.end());
            if (haplotype_info[edge.getInnerId()].haplotype == haplotype_) {
                if (isSimpleBridge(edge) && edge.fullSize() < saved_bridge_cutoff) {
                    bridges ++;
                    logger_.info() << "Skipping getEdge " << edge.getId() << " as bridge\n";
                    continue;
                }
                removed_len += edge.fullSize();
                to_delete.insert(edge.getId());
                to_delete.insert(edge.rc().getId());
                logger_.trace() << "removing " << edge.getId()  << " label ";
                for(Edge::id_type id : haplotype_info[edge.getInnerId()].label) {
                    logger_ << id << "_";
                }
                logger_ << std::endl;
                removed ++;
                changed = true;
            } else {
                logger_.trace() << "skipping edge " << edge.getId() << " " << haplotype_info[edge.getInnerId()].haplotype <<std::endl;
            }
        }
        for(EdgeId eid : to_delete) {
            deleteEdgeHaplo(eid);
        }
        compressAllVertices();
    }
    logger_.info() << "Saved " << bridges << "bridges\n";
    logger_.info() << "Removed " << removed << " edges of haplo " << haplotype_  << " total truncLen " << removed_len <<std::endl;
}

void HaplotypeRemover::process() {
    std::vector<std::pair<EdgeId, EdgeId>> bulges = getBulgeLabels();
    logger_ << "Detected " << bulges.size() << " bulges" << std::endl;
    updateFixableHaplotypes(bulges);
    updateAmbiguousHaplotypes(bulges);
    removeHaplotype();
    logger_.debug() << "removed \n";
    MultiGraphHelper::printEdgeGFA(mg, out_dir / "before_clean.gfa", true);
    cleanGraph();
    MultiGraphHelper::printEdgeGFA(mg, out_dir / "mdbg.hpc.gfa", true);

}


