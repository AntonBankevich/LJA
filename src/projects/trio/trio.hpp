#pragma once


#include <array>
#include <vector>
#include <unordered_map>
#include <sequences/seqio.hpp>
#include "common/logging.hpp"
#include <common/string_utils.hpp>
#include "dbg/multi_graph.hpp"
#include "haplo_stats.hpp"

namespace trio {
using namespace multigraph;
struct HaplotypeRemover {
    std::unordered_map<Edge::id_type, HaplotypeStats> haplotype_info;
    MultiGraph &mg;
//    haplo_map_type haplotypes;
    logging::Logger &logger_;
    size_t threads;
    Haplotype haplotype_;
    static const size_t MAX_TIP_LENGTH = 1000000;
    static constexpr double BULGE_MULTIPLICATIVE_CUTOFF = 1.2;
//Bridges of wrong haplotype longer that this cutoff are deleted, shorter are saved;
    const size_t saved_bridge_cutoff;
    std::experimental::filesystem::path out_dir;

    HaplotypeRemover(logging::Logger &logger, size_t threads, multigraph::MultiGraph &mg,
                     const std::experimental::filesystem::path &haployak, const Haplotype haplotype,
                     const std::experimental::filesystem::path &out_dir, const size_t saved_bridge_cutoff);

    void updateEdgeHaplotype(Edge &edge, Haplotype h) {
        haplotype_info[edge.getInnerId()].haplotype = h;
        haplotype_info[edge.rc().getInnerId()].haplotype = h;
    }

    void process();

    void deleteEdgeHaplo(multigraph::EdgeId eid);

    void compressAllVertices();

    void cleanGraph();

    std::vector<std::pair<EdgeId, EdgeId>> getBulgeLabels();

    void updateFixableHaplotypes(const std::vector<std::pair<EdgeId, EdgeId>> &bulges);

    void updateAmbiguousHaplotypes(const std::vector<std::pair<EdgeId, EdgeId>> &bulges);

    void removeHaplotype();
};

}
