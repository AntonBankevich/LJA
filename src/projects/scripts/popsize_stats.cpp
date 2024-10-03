//
// Created by enzo on 7/29/24.
//
#include <dbg/graph_algorithms.hpp>
#include <common/logging.hpp>
#include <error_correction/diploidy_analysis.hpp>
#include <alignment/ksw_aligner.hpp>

int main(int argc, char* argv[]){
    std::vector<std::experimental::filesystem::path> lib;
    for (int i = 1; i < argc; ++i) {
        lib.push_back(argv[i]);
    }
    int threads = 20;
    logging::Logger logger;
    hashing::RollingHash hasher(3000);
    dbg::SparseDBG dbg = dbg::LoadDBGFromEdgeSequences(logger, threads, lib, hasher);
    std::cout << "graph was read";
    dbg::BulgePathFinder finder(dbg, -1.0);
    KSWAligner kswAligner(1, 5, 5, 3);
    logger << "paths found: " << finder.paths.size() << std::endl;
    int path_id = 1;
    omp_set_num_threads(threads);
    struct edgerec {
        bool bulge;
        size_t first;
        size_t second;
    };
    std::ofstream os("bulge_paths.stats");
    for(dbg::BulgePath<dbg::DBGTraits> &path : finder.paths) {
        if(path.size() == 1)
            continue;
        std::vector<edgerec> output;
        os << "#path" << std::endl;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(path, output, kswAligner)
        for(int i = 0; i < path.size(); ++i) {
            bool bulge;
            const std::pair<ag::BaseEdge<dbg::DBGTraits> *, ag::BaseEdge<dbg::DBGTraits> *> &edge_pair = path[i];
            ag::BaseEdge<dbg::DBGTraits> &edge1 = *edge_pair.first;
            ag::BaseEdge<dbg::DBGTraits> &edge2 = *edge_pair.second;
            if (edge1 == edge2) bulge = false;
            output[i] = {bulge, edge1.truncSize(), edge2.truncSize()};
        }
        for (auto rec : output) {
            os << (rec.bulge ? "bulge" : "homo") << rec.first << '\t' << rec.second << std::endl;
        }
    }
    os.close();
}
