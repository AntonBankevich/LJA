#include "decision_rules.hpp"
#include "supregraph.hpp"
#include "multiplexer.hpp"
#include "read_storage.hpp"
#include "converter.hpp"
#include "unique_vertex_storage.hpp"
#include <dbg/dbg_construction.hpp>
#include <dbg/graph_algorithms.hpp>
#include <common/pipeline_tools.hpp>
#include <dbg/dbg_read_alignment_storage.hpp>
#include <error_correction/multiplicity_estimation.hpp>
#include <dbg/aln_reads_reader.hpp>

using namespace spg;

namespace ag {
    template<class Traits>
    inline void printEdge(std::ostream &os, const typename Traits::Edge &edge, const std::string &extra_label = "",
                          const std::string &color = "black") {
        const typename Traits::Vertex &end = edge.getFinish();
        os << "\"" << edge.getStart().getId() << "\" -> \"" << end.getId() <<
           "\" [label=\"" << edge.getInnerId() << " " << edge.firstNucl() << " " << edge.truncSize() << "\"";
        if (!extra_label.empty()) {
            os << " labeltooltip=\"" << extra_label << "\"";
//        os << "\\n"<<extra_label;
        }
        os << " color=\"" + color + "\"]\n";
    }

    template<class Traits>
    inline void printDot(std::ostream &os, AssemblyGraph<Traits> &component,
                         std::function<std::string(Edge &)> labeler = [](Edge &){return "";},
                         std::function<std::string(Vertex &)> vertex_colorer = [](Vertex &){return "white";}) {
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_set<dbg::VertexId> extended;
        for (typename Traits::Vertex &vert: component.vertices()) {
            std::string color = "white";
            os << vert.getId();
            os << " [style=filled fillcolor=\"" + vertex_colorer(vert) + "\"";
            if (vert.size() < 10)
                os << " label=" << vert.getSeq();
            else
                os << " label=\"" << vert.getId() << " " << vert.size() << "\"";
            os << "]\n";
        }
        for (typename Traits::Edge &edge: component.edges()) {
            printEdge<Traits>(os, edge, labeler(edge));
        }
        os << "}\n";
    }

    template<class Traits>
    inline void printDot(const std::experimental::filesystem::path &path, AssemblyGraph<Traits> &component,
                         std::function<std::string(Edge &)> labeler = [](Edge &){return "";},
                         std::function<std::string(Vertex &)> vertex_colorer = [](Vertex &){return "white";}) {
        std::ofstream os;
        os.open(path);
        printDot(os, component, std::move(labeler), vertex_colorer);
        os.close();
    }
}

std::unordered_map<std::string, std::experimental::filesystem::path>
RunMultiplexing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, size_t k, size_t w,
                  const io::Library &graph_gfa, const io::Library &reads_file, bool debug) {
    hashing::RollingHash hasher(k, 239);
    logger.info() << "Loading graph" << std::endl;
    dbg::SparseDBG dbg = graph_gfa.empty() ? DBGPipeline(logger, hasher, w, reads_file, dir, threads) :
            dbg::LoadDBGFromEdgeSequences(logger, threads, graph_gfa, hasher);
//    IdIndex<dbg::Vertex> index(dbg.vertices().begin(), dbg.vertices().end());
    logger.info() << "Loading reads" << std::endl;
    dbg::SeqReader reader(reads_file, logger, threads);
    dbg::DBGAlignedReadStorage dbg_storage(logger, threads, dbg,
                                           dbg::AlignReads(logger, threads, reader.begin(), reader.end(), dbg, w),
                                           true);
    dbg_storage.trackSuffixes(logger, threads, dbg, 0, 10000000);
    logger.info() << "Converting graph" << std::endl;
    spg::SPGConverter<dbg::DBGTraits> converter;
    spg::SupreGraph graph = converter.convert(dbg);
    ag::LoggingListener<SPGTraits> modificationLogger(graph, logger.getLoggerStream(logging::LogLevel::trace));
    logger.info() << "Printing initial graph" << std::endl;
    ag::printDot(dir/"supregraph1.dot", graph);
    logger.info() << "Converting reads" << std::endl;
    ag::AlignedReadStorage<SPGTraits> reads(logger, threads, graph, converter.convertLib(dbg_storage, graph));
    logger.info() << "Constructing path index" << std::endl;
    ag::SuffixTracker<SPGTraits> suffixTracker(reads, graph, 0, 1000000000);
    suffixTracker.fillFromStorage(logger, threads);
    logger.info() << "Reconstructing uniqueness" << std::endl;
    UniqueClassificator classificator(dbg, dbg_storage, 0, false, false);
    classificator.classify(logger, 40000, dir / "mult");
    UniqueVertexStorage unique_storage(graph);
    for(dbg::Edge &edge : dbg.edgesUnique()) {
        if(edge.isOuter() && classificator.isUnique(edge))
            unique_storage.add(converter.map(edge));
    }
    logger.info() << "Multiplexing" << std::endl;
//    ChainRule rule(path_index, 4000);
    AndreyRule rule(suffixTracker, unique_storage);
    spg::Multiplexer multiplexer(graph, reads, rule, 200000);
//    multiplexer.fullMultiplex(logger, threads);
    size_t cnt = 0;
    while(!multiplexer.finished()) {
        multiplexer.multiplex(logger, threads);
        ag::printDot(dir/("supregraph" + itos(cnt, 3) + ".dot"), graph, suffixTracker.labeler(), unique_storage.getColorer("white", "green"));
        cnt++;
    }
    graph.removeMarked();
    logger.info() << "Printing final graph" << std::endl;
    ag::printDot(dir/"supregraph_final.dot", graph);
    return {{"supregraph_initial", dir / "supregraph1.dot"}, {"supregraph_final", dir / "supregraph_final.dot"}};
}

class SupreGraphPhase : public Stage {
public:
    SupreGraphPhase() : Stage(AlgorithmParameters(
            {"k-mer-size=5001", "window=500"},
            {}, ""), {"graph", "reads"}, {"supregraph_initial", "supregraph_final"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoull(parameterValues.getValue("k-mer-size"));
        size_t w = std::stoull(parameterValues.getValue("window"));
        return RunMultiplexing(logger, threads, dir, k, w,input.at("graph"), input.at("reads"), debug);
    }
};
int main(int argc, char **argv) {
    SupreGraphPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads", "k=k-mer-size", "w=window"}, {});
    LoggedProgram multiplexing("Multiplexing", std::move(phase), std::move(parser),
                               "Starting multiplexing procedure", "Finished multiplexing procedure");
    multiplexing.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}