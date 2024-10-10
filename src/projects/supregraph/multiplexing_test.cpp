#include "decision_rules.hpp"
#include "supregraph.hpp"
#include "multiplexer.hpp"
#include "read_storage.hpp"
#include "converter.hpp"
#include "unique_vertex_storage.hpp"
#include <dbg/dbg_construction.hpp>
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
RunMultiplexing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, size_t k, size_t w,
                  const io::Library &graph_gfa, const io::Library &reads_file, bool debug) {
    hashing::RollingHash hasher(k, 239);
    logger.info() << "Loading graph" << std::endl;
    dbg::SparseDBG dbg = graph_gfa.empty() ? DBGPipeline(logger, hasher, w, reads_file, dir, threads) :
            dbg::LoadDBGFromEdgeSequences(logger, threads, graph_gfa, hasher);
//    IdIndex<dbg::Vertex> index(dbg.vertices().begin(), dbg.vertices().end());
    logger.info() << "Loading reads" << std::endl;
    dbg::ReadAlignmentStorage storage(dbg, 0, 1000000, true, false, true);
//    dbg::ReadAlignmentStorage extra_reads(dbg, 0, 1000000, false, false, false);
    io::SeqReader reader(reads_file);
    dbg::KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, 500);
    storage.FillAlignments(logger, threads, reader.begin(), reader.end(), dbg, index);
//    ag::LoadAllReads<dbg::DBGTraits>(reads_file, {&storage, &extra_reads}, index);
    logger.info() << "Converting graph" << std::endl;
    spg::SPGConverter<dbg::DBGTraits> converter;
    spg::SupreGraph graph = converter.convert(dbg);
    ag::LoggingListener<SPGTraits> modificationLogger(graph, logger.getLoggerStream(logging::LogLevel::trace));
//    TODO: remove after listener system is normalized
//    for(Vertex &v : graph.vertices()) graph.fireAddVertex(v);
    for(Edge &edge : graph.edges()) graph.fireAddEdge(edge);
    logger.info() << "Printing initial graph" << std::endl;
    ag::printDot(dir/"supregraph1.dot", graph);
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
//    ChainRule rule(path_index, 4000);
    AndreyRule rule(path_index, unique_storage);
    spg::Multiplexer multiplexer(graph, rule, 200000);
    multiplexer.fullMultiplex(logger, threads, path_index);
    ag::printDot(dir/"supregraph2.dot", graph);
    std::vector<GraphPath> unbranching;
    for(Vertex &v : graph.vertices()) {
        if(!v.isJunction())
            continue;
        for(Edge &start : v) {
            GraphPath path = GraphPath::WalkForward(start);
            if(path.start() > path.finish().rc() || (path.start() == path.finish().rc() && path.Seq() > path.RC().Seq()))
                continue;
            if(path.size() > 2 || (path.size() == 2 && path.frontEdge().isPrefix() == path.backEdge().isPrefix()))
                unbranching.emplace_back(std::move(path));
        }
    }
    for(GraphPath &path : unbranching) {
        graph.mergePath(path);
    }
    unbranching.clear();
    for(Vertex &start : graph.vertices()) {
        if(start.isJunction() || start.front().getFinish().isJunction())
            continue;
        GraphPath forward = GraphPath ::WalkForward(start.front());
        if(forward.finish() != start)
            continue;
        bool ok = true;
        for (Vertex &v: forward.vertices()) {
            if (start.isPalindrome() == v.isPalindrome() && (v < start || v.rc() < start)) {
                ok = false;
                break;
            }
        }
        if(!ok)
            continue;
        unbranching.emplace_back(std::move(forward));
    }
    for(GraphPath &path : unbranching) {
        graph.mergePath(path);
    }
    ag::printDot(dir/"supregraph3.dot", graph);
    spg::Multiplexer multiplexer2(graph, rule, 200000000);
    multiplexer2.fullMultiplex(logger, threads, path_index);
    logger.info() << "Printing final graph" << std::endl;
    ag::printDot(dir/"supregraph4.dot", graph);
    return {{"supregraph_initial", dir / "supregraph1.dot"}, {"supregraph_final", dir / "supregraph4.dot"}};
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