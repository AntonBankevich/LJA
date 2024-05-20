#include <supregraph/multiplexer.hpp>
#include "supregraph/read_storage.hpp"
#include <supregraph/supregraph.hpp>
#include <dbg/multi_graph.hpp>
#include <dbg/visualization.hpp>
#include <dbg/dbg_construction.hpp>
#include "common/pipeline_tools.hpp"
#include "converter.hpp"
#include "assembly_graph.hpp"

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

class SPTransformPhase : public Stage {
public:
    SPTransformPhase() : Stage(AlgorithmParameters({}, {}, ""), {"graph", "reads"}, {"graph"}) {
    }

private:

protected:
    std::unordered_map<std::string, std::experimental::filesystem::path>
    innerRun(logging::Logger &logger, size_t threads,
             const std::experimental::filesystem::path &dir, bool debug,
             const AlgorithmParameterValues &parameterValues,
             const std::unordered_map<std::string, io::Library> &input) override {
        if(input.find("graph")->second.empty()) {
            logger.info() << "Constructing graph from reads" << std::endl;
            io::Library lib = oneline::initialize<std::experimental::filesystem::path>(input.find("reads")->second);
            logger.info() << lib << std::endl;
            hashing::RollingHash hasher(5001, 239);
            size_t w = 500;
            dbg::SparseDBG initial = DBGPipeline(logger, hasher, w, lib, dir, threads);
            logger.info() << "Initial graph parameters: " << initial.size() << " vertices, " << initial.edgeCount() << " edges"
                          << std::endl;
            ag::ReadLogger readLogger(threads, dir/"read_log.txt");
            dbg::ReadAlignmentStorage readStorage(initial, 0, 0, false, false, false);
            io::SeqReader reader(lib);
            dbg::KmerIndex index(initial);
            index.fillAnchors(logger, threads, initial, 500);
            readStorage.FillAlignments(logger, threads, reader.begin(), reader.end(), initial, index);
            SPGConverter<dbg::DBGTraits> converter;
            spg::SupreGraph graph = converter.Convert(initial);
            std::vector<spg::EdgeId> outerEdges;
            for(spg::Edge &edge : graph.edgesUnique()) {
                if(edge.isOuter())
                    outerEdges.emplace_back(edge.getId());
            }
            for(spg::EdgeId eid: outerEdges) {
                graph.outerEdgeToVertex(*eid);
            }
            ag::printDot<dbg::DBGTraits>(dir / "initial.dot", initial);
            ag::printDot<spg::SPGTraits>(dir / "converted.dot", graph);
            logger.info() << "Aligning reads" << std::endl;
            spg::PathStorage pathStorage(graph);
            for(ag::AlignedRead<dbg::DBGTraits> &read : readStorage) {
                ag::GraphPath<dbg::DBGTraits> opath = read.path.unpack();
                spg::GraphPath path = converter.convertPath(opath);
                if(path.size() >= 2) {
                    std::cout << read.id << std::endl;
                    for(spg::Edge &edge : path.edges()) {
                        std::cout << edge.getId() << " ";
                    }
                    std::cout << std::endl;
                }
                pathStorage.addRead(read.id, {path});
            }
            spg::ChainRule rule(pathStorage, 5001);
            spg::Multiplexer multiplexer(graph, rule);
            recreate_dir(dir/"figs");
            logger.info() << "Multiplexing graph" << std::endl;
            for(size_t cnt = 0; cnt < 100 && !multiplexer.finished(); ++cnt) {
                logger << "Step " << cnt << std::endl;
                multiplexer.multiplex();
                ag::printDot<spg::SPGTraits>(dir / "figs"/ ("graph" + itos(cnt) + ".dot"), graph);
            }
            ag::printDot<spg::SPGTraits>(dir / "resolved.dot", graph);
        } else {
            const std::experimental::filesystem::path &graph_path = input.find("graph")->second.front();
            multigraph::MultiGraph initial = multigraph::MultiGraphHelper::LoadGFA(graph_path, false);
            initial = multigraph::MultiGraphHelper::TransformToEdgeGraph(initial, 501);
            logger.info() << "Initial graph parameters: " << initial.size() << " vertices, " << initial.edgeCount() << " edges"
                          << std::endl;
            SPGConverter<multigraph::MGTraits> converter;
            spg::SupreGraph graph = converter.Convert(initial);
            ag::printDot<multigraph::MGTraits>(dir / "initial.dot", initial);
            ag::printDot<spg::SPGTraits>(dir / "converted.dot", graph);
        }
        return {{"graph",    dir / "converted.dot"}};
    }
};


int main(int argc, char **argv) {
    SPTransformPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram program("conversion", std::move(phase), std::move(parser), "Starting graph conversion", "Finished graph conversion");
    program.run(oneline::initialize<std::string, char*>(argv, argv + argc));
    return 0;
}
