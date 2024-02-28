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
        dbg::SparseDBG dbg = DBGPipeline(logger, hasher, w, {contigs_path}, dir, threads);
        std::cout << "graph built" << std::endl;
        printDot(dir / "graph.dot", dbg::Component(dbg));
//        dbg::printGFA()
        BulgePathFinder finder(dbg, -1.0);
        KSWAligner kswAligner(1, 5, 5, 3);
        AlignmentForm diploAln;
        logger << "paths found: " << finder.paths.size() << std::endl;
        for(BulgePath &path : finder.paths) {
            if(path.size() == 1)
                continue;
            for(auto &edge_pair : path) {
                dbg::Edge &edge1 = *edge_pair.first;
                dbg::Edge &edge2 = *edge_pair.second;
                if(edge1 == edge2) {
                    diploAln += AlignmentForm({CigarPair('M', edge1.truncSize())});
                } else {
                    diploAln += kswAligner.globalAlignment(edge1.truncSeq().str(), edge2.truncSeq().str());
                }
            }
        }
        int i = 0;
        bool windowDivergent = false;
        std::vector<bool> divergence;
        for(AlignmentForm::AlignmentColumn alnColumn : diploAln.columns()) {
            if ((i % 100 == 0) && (i>0)) {
                divergence.push_back(windowDivergent);
                windowDivergent = false;
            }
            if (alnColumn.event != 'M') {
                logger << "not match\n";
                windowDivergent = true;}
            ++i;
        }
        if (i % 100 != 0){divergence.push_back(windowDivergent);}

        int T = 0;
        int K = 0;
        std::ofstream os;
        os.open(dir / "divergence.psmcfa");
        os << "> testseq\n";
        for (auto window : divergence){
            if (window){
                os << 'T';
                ++T;
            } else {
                os << 'K';
                ++K;
            }
        }
        os.close();
        logger << "divergent: " << T << "\nsimilar: " << K << std::endl; 
        std::unordered_map<std::string, std::experimental::filesystem::path> results;
        results["bitvector"] = dir / "divergence.psmcfa";
        return std::move(results);
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