#include "decision_rules.hpp"
#include "supregraph.hpp"
#include "multiplexer.hpp"
#include "read_storage.hpp"
#include "converter.hpp"
#include "unique_vertex_storage.hpp"
#include <dbg/graph_algorithms.hpp>
#include <common/pipeline_tools.hpp>
#include <dbg/graph_alignment_storage.hpp>
#include <error_correction/multiplicity_estimation.hpp>

using namespace spg;

namespace ag {
    template<class Traits>
    inline void printEdge(std::ostream &os, const typename Traits::Edge &edge, const std::string &extra_label = "",
                          const std::string &color = "black") {
        const typename Traits::Vertex &end = edge.getFinish();
        os << "\"" << edge.getStart().getId() << "\" -> \"" << end.getId() <<
           "\" [label=\"" << edge.getInnerId() << " " << edge.nuclLabel() << " " << edge.truncSize() << "\"";
        if (!extra_label.empty()) {
            os << " labeltooltip=\"" << extra_label << "\"";
//        os << "\\n"<<extra_label;
        }
        os << " color=\"" + color + "\"]\n";
    }

    template<class Traits>
    inline void printDot(std::ostream &os, const AssemblyGraph<Traits> &component) {
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_set<dbg::VertexId> extended;
        for (const typename Traits::Vertex &vert: component.vertices()) {
            std::string color = "white";
            os << vert.getId();
            os << " [style=filled fillcolor=\"" + color + "\"";
            if (vert.size() < 10)
                os << " label=" << vert.getSeq();
            else
                os << " label=\"" << vert.getId() << " " << vert.size() << "\"";
            os << "]\n";
        }
        for (const typename Traits::Edge &edge: component.edges()) {
            printEdge<Traits>(os, edge);
        }
        os << "}\n";
    }

    template<class Traits>
    inline void printDot(const std::experimental::filesystem::path &path, const AssemblyGraph<Traits> &component) {
        std::ofstream os;
        os.open(path);
        printDot(os, component);
        os.close();
    }
}

std::unordered_map<std::string, std::experimental::filesystem::path>
RunMultiplexing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                  const io::Library &graph_gfa, const std::experimental::filesystem::path &reads_file, bool debug) {
    hashing::RollingHash hasher(5001, 239);
    logger.info() << "Loading graph" << std::endl;
    dbg::SparseDBG dbg = dbg::LoadDBGFromEdgeSequences(logger, threads, graph_gfa, hasher);
    IdIndex<dbg::Vertex> index(dbg.vertices().begin(), dbg.vertices().end());
    logger.info() << "Loading reads" << std::endl;
    dbg::ReadAlignmentStorage storage(dbg, 0, 1000000, true, false, true);
    dbg::ReadAlignmentStorage extra_reads(dbg, 0, 1000000, false, false, false);
    ag::LoadAllReads<dbg::DBGTraits>(reads_file, {&storage, &extra_reads}, index);
    logger.info() << "Converting graph" << std::endl;
    spg::SPGConverter<dbg::DBGTraits> converter;
    spg::SupreGraph graph = converter.convert(dbg);
    logger.info() << "Printing initial graph" << std::endl;
    ag::printDot(dir/"supregraph_initial.dot", graph);
    logger.info() << "Converting reads" << std::endl;
    spg::PathStorage paths = converter.convertLibs({&storage}, graph);
    logger.info() << "Constructing path index" << std::endl;
    spg::PathIndex path_index(graph, paths);
    logger.info() << "Reconstructing uniqueness" << std::endl;
    UniqueClassificator classificator(dbg, storage, 0, false, false);
    classificator.classify(logger, 40000, dir / "mult");
    UniqueVertexStorage unique_storage(graph);
    for(dbg::Edge &edge : dbg.edgesUnique()) {
        if(edge.isOuter() && classificator.isUnique(edge))
            unique_storage.add(converter.map(edge));
    }
    logger.info() << "Multiplexing" << std::endl;
    AndreyRule rule(path_index, unique_storage);
    spg::Multiplexer multiplexer(graph, rule);
    multiplexer.fullMultiplex(path_index);
    logger.info() << "Printing final graph" << std::endl;
    ag::printDot(dir/"supregraph_final.dot", graph);
    return {{"supregraph_initial", dir / "supregraph_initial.dot"}, {"supregraph_final", dir / "supregraph_final.dot"}};
}

class SupreGraphPhase : public Stage {
public:
    SupreGraphPhase() : Stage(AlgorithmParameters(
            {},
            {}, ""), {"graph", "reads"}, {"supregraph_initial", "supregraph_final"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        return RunMultiplexing(logger, threads, dir, input.at("graph"), input.at("reads")[0], debug);
    }
};
int main(int argc, char **argv) {
    SupreGraphPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram multiplexing("Multiplexing", std::move(phase), std::move(parser), "Starting multiplexing procedure", "Finished multiplexing procedure");
    multiplexing.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}