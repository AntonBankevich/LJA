//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#ifndef DR_READSALIGNER_H
#define DR_READSALIGNER_H

#include <dbg/multi_graph.hpp>
#include <nano/GraphContig.hpp>


namespace nano {

    class ReadsAlignerGA {
    public:
        ReadsAlignerGA(const multigraph::MultiGraph &mg) : mg_(mg) {}

        void Align(const std::unordered_map<std::string, Contig> &sequences,
                   const std::experimental::filesystem::path &graph,
                   const std::experimental::filesystem::path &output_dir,
                   const size_t threads,
                   const int batch_num);

        std::unordered_map<std::string, std::vector<GraphContig>> ExtractPaths(
                const std::experimental::filesystem::path &output_dir,
                int batch_num);

    private:
        std::experimental::filesystem::path SaveBatch(const std::unordered_map<std::string, Contig> &sequences,
                                                      const std::experimental::filesystem::path &output_dir,
                                                      int batch_num);

        GraphContig ExtractAlignment(const std::string &ln,
                                     const std::unordered_map<std::string, Contig> &sequences);
        std::unordered_map<std::string, std::vector<GraphContig>> ExtractBestPaths(const std::experimental::filesystem::path &batch_gaf,
                                                                                   const std::unordered_map<std::string, Contig> &sequences);

        const multigraph::MultiGraph &mg_;

    };

}
#endif //DR_READSALIGNER_H