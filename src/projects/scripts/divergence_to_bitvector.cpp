#include <common/pipeline_tools.hpp>
#include <dbg/graph_algorithms.hpp>
#include <error_correction/diploidy_analysis.hpp>
#include <alignment/ksw_aligner.hpp>
#include <dbg/dbg_construction.hpp>
#include <assembly_graph/visualization.hpp>
#include "dbg/graph_printing.hpp"
using namespace dbg;

struct edgerec {
    bool bulge;
    size_t first;
    size_t second;
    std::string first_id;
    std::string second_id;
};

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
        logger << "graph built";
        ObjInfo<dbg::Vertex> vertexInfo = VertexPrintStyles<dbg::DBGTraits>::defaultDotInfo();
        ObjInfo<dbg::Edge> edgeInfo = EdgePrintStyles<dbg::DBGTraits>::defaultDotInfo();
        Printer<DBGTraits> printer(vertexInfo, edgeInfo);
        printer.printGFA(dir / "graph.gfa", dbg::Component(dbg));
        printer.printDot(dir / "graph.dot", dbg::Component(dbg));
        BulgePathFinder finder(dbg, -1.0);
        KSWAligner kswAligner(1, 5, 5, 3);
        logger << "paths found: " << finder.paths.size() << std::endl;
        std::ofstream outfile;
        std::ofstream outstats(dir/ "bulge_paths.stats");
        outfile.open(dir / "divergence.psmcfa");
        int path_id = 1;
        omp_set_num_threads(threads);
        for(BulgePath<dbg::DBGTraits> &path : finder.paths){
            if(path.size() == 1)
                continue;
            std::vector<std::vector<bool>> pathDivergence(path.size());
            std::vector<edgerec> statrecs(path.size());
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(path, pathDivergence, statrecs, kswAligner)
            for(int i = 0; i < path.size(); ++i) {
                const std::pair<dbg::EdgeId, dbg::EdgeId> &edge_pair = path[i];
                dbg::Edge &edge1 = *edge_pair.first;
                dbg::Edge &edge2 = *edge_pair.second;
                if(edge1 == edge2) {
                    pathDivergence[i] = std::vector<bool>(edge1.truncSize(), false);
                    statrecs[i] = {false, edge1.truncSize(), edge1.truncSize(),
                                     edge1.getInnerId().str(), edge1.getInnerId().str()};
                } else {
                    AlignmentForm alignment = kswAligner.globalAlignment(edge1.truncSeq().str(), edge2.truncSeq().str());
                    pathDivergence[i] = alignmentToBitvector(alignment);
                    statrecs[i] = {true, edge1.truncSize(), edge2.truncSize(),
                                     edge1.getInnerId().str(), edge2.getInnerId().str()};
                }
            }
            pathDivergenceToFile(outfile, path_id, pathDivergence);
            writeStats(outstats, path_id, statrecs);
            ++path_id;
        }
        outfile.close();
        outstats.close();
        std::unordered_map<std::string, std::experimental::filesystem::path> results;
        results["bitvector"] = dir / "divergence.psmcfa";
        return std::move(results);
    }

    void writeStats(std::ofstream &os, int path_id, std::vector<edgerec> &statrecs) {
        os << "#path_" << path_id << "\n";
        for(const auto &rec: statrecs) {
            os << (rec.bulge ? "bulge\t" : "single\t");
            os << rec.first_id << '\t' << rec.first << '\t';
            os << rec.second_id << '\t' << rec.second << '\n';
        }
    }

    void pathDivergenceToFile(std::ofstream &os, int path_id,
                              std::vector<std::vector<bool>> pathDivergence){
        std::vector<bool> pd;
        int pos = 0;
        bool window = false;
        os << ">path_" << path_id << std::endl;
        for (const auto &segment : pathDivergence) {
            for(const auto &bp : segment) {
                if(pos % 100 == 0) {
                    os << (window ? 'K' : 'T');
                    window = false;
                }
                if(bp) window = true;
                pos++;
            }
        }
        if (pos % 100 != 0) os << (window ? 'K' : 'T');
        os << std::endl;
    }

    std::vector<bool> alignmentToBitvector(const AlignmentForm &alignment){
        std::vector<bool> result;
        for(const AlignmentForm::AlignmentColumn &alnColumn : alignment.columns()) {
            result.emplace_back(alnColumn.event != 'M');
        }
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