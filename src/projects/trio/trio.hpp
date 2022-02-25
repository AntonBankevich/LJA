#pragma once


#include <array>
#include <vector>
#include <unordered_map>
#include <sequences/seqio.hpp>
#include "common/logging.hpp"
#include <common/string_utils.hpp>
#include "lja/multi_graph.hpp"
#include "lja/subdataset_processing.hpp"



struct HaplotypeRemover {
    multigraph::MultiGraph &mg;
    haplo_map_type haplotypes;
    logging::Logger &logger_;
    char haplotype_;
    std::experimental::filesystem::path out_dir;
    std::unordered_map<std::string, std::string> bulges;
    HaplotypeRemover(logging::Logger &logger, multigraph::MultiGraph &mg,
                     const std::experimental::filesystem::path &haployak, const char haplotype, const std::experimental::filesystem::path &out_dir):logger_(logger), mg(mg), haplotype_(haplotype), out_dir(out_dir) {

        auto bulges = getBulgeLabels();
        logger_.info() << "got "<< bulges.size() << " bulges\n";
//to_separate_fucntion
        haplo_map_type haplotypes;
        string s;
        std::ifstream haplo_file(haployak);
        while (std::getline(haplo_file, s)){
            HaplotypeStats h(s);
            haplotypes[h.label] = h;
        }
        std::string out_name = "haplotype_";
        ensure_dir_existance(out_dir);
    }


    void process();
    void deleteEdgeHaplo(int eid);
    void cleanGraph();
    std::unordered_map<std::string, std::string> getBulgeLabels();
    void updateFixableHaplotypes();
    void updateAmbiguousHaplotypes();
    void removeHaplotype();
};

std::experimental::filesystem::path simplifyHaplo(logging::Logger &logger, size_t threads,
                                           const std::experimental::filesystem::path &output_file,
                                           const std::experimental::filesystem::path &diplo_graph,
                                           const std::experimental::filesystem::path &haployak,
                                           const char haplotype, const std::experimental::filesystem::path &corrected_reads, io::Library & reads,const std::experimental::filesystem::path &dir);

