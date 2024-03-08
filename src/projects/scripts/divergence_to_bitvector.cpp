#include <common/pipeline_tools.hpp>
#include <dbg/graph_algorithms.hpp>
#include <error_correction/diploidy_analysis.hpp>
#include <alignment/ksw_aligner.hpp>
#include <dbg/dbg_construction.hpp>
#include <dbg/visualization.hpp>
#include "dbg/graph_printing.hpp"


class DivergenceToBitVector : public Stage {
public:
    DivergenceToBitVector() : Stage(AlgorithmParameters({"kmer-size=", "window-size="}, {"reads"}, ""), {"contigs"}, {"bitvector"}) {
    }

    virtual std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                  const std::experimental::filesystem::path &dir, bool debug,
                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override{
        std::experimental::filesystem::path contigs_path = parameterValues.getListValue("contigs").front();
        size_t k = std::stoull(parameterValues.getValue("kmer-size"));
        size_t w = std::stoull(parameterValues.getValue("window-size"));
        hashing::RollingHash hasher(k);
        logger << "starting DBG build\n";
        dbg::SparseDBG dbg = DBGPipeline(logger, hasher, w, {contigs_path}, dir, threads);
        printDot(dir / "graph.dot", dbg::Component(dbg));
        logger << "graph built";
//        dbg::printGFA()
        BulgePathFinder finder(dbg, -1.0);
        KSWAligner kswAligner(1, 5, 5, 3);
        logger << "paths found: " << finder.paths.size() << std::endl;
        std::ofstream outfile;
        outfile.open(dir / "divergence.psmcfa");
        int path_id = 1;
        omp_set_num_threads(threads);
        for(BulgePath &path : finder.paths){
            if(path.size() == 1)
                continue;
            std::vector<std::vector<bool>> pathDivergence(path.size());
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(path, pathDivergence, kswAligner)
            for(int i = 0; i < path.size(); ++i) {
                const std::pair<dbg::Edge *, dbg::Edge *> &edge_pair = path[i];
                dbg::Edge &edge1 = *edge_pair.first;
                dbg::Edge &edge2 = *edge_pair.second;
                if(edge1 == edge2) {
                    pathDivergence[i] = std::vector<bool>(edge1.truncSize(), false);
                } else {
                    AlignmentForm alignment = kswAligner.globalAlignment(edge1.truncSeq().str(), edge2.truncSeq().str());
                    pathDivergence[i] = alignmentToBitvector(alignment);
                }
            }
            pathDivergenceToFile(outfile, path_id, pathDivergence);
            ++path_id;
        }
        outfile.close();
        std::unordered_map<std::string, std::experimental::filesystem::path> results;
        results["bitvector"] = dir / "divergence.psmcfa";
        return std::move(results);
    }

    void pathDivergenceToFile(std::ofstream &os, int path_id,
                              std::vector<std::vector<bool>> pathDivergence){
        os << ">path_" << path_id << std::endl;
        for (const auto &edgeDivergence: pathDivergence){
            for (auto window : edgeDivergence){
                os << (window ? 'K' : 'T');
            }
        }
        os << std::endl;
    }

    std::vector<bool> alignmentToBitvector(const AlignmentForm &alignment){
        int i = 0;
        std::vector<bool> result;
        bool windowDivergent = false;
        for(const AlignmentForm::AlignmentColumn &alnColumn : alignment.columns()){
            if ((i % 100 == 0) && (i>0)) {
                result.push_back(windowDivergent);
                windowDivergent = false;
            }
            if (alnColumn.event != 'M') {
                windowDivergent = true;}
            ++i;
        }
        if (i % 100 != 0){result.push_back(windowDivergent);}
        return result;
    }
};


int main(int argc, char **argv) {
    DivergenceToBitVector phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads", "k=kmer-size", "w=window-size"});
    LoggedProgram consensus("DivergenceToBitVector", std::move(phase), std::move(parser), "Starting divergence processing", "Finished bit vector construction");
    consensus.run(oneline::initialize<std::string, char*>(argv, argv + argc));
    return 0;
}