//
// Created by Tatiana Dvorkina on 11.04.2022.
//


#include <dbg/graph_algorithms.hpp>
#include <dbg/minimizer_selection.hpp>
#include "dbg/sparse_dbg.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <vector>
#include <dbg/graph_stats.hpp>
#include <nano/ReadsAligner.h>
#include <nano/SGraphBuilder.h>
#include <lja/multi_graph.hpp>

using namespace dbg;

int runNano(const std::experimental::filesystem::path &graph,
            const  std::experimental::filesystem::path &ont_reads,
            const  std::experimental::filesystem::path &unique_edges,
            const size_t &threads,
            const std::experimental::filesystem::path &dir){
    const unsigned int BATCH_SIZE = 5;
    io::SeqReader reader(ont_reads);
    std::cerr << "Seq reader loaded" << std::endl;
    multigraph::MultiGraph mmg;
    mmg.LoadGFA(graph, true);
    std::cerr << "GFA loaded" << std::endl;
    multigraph::MultiGraph mg = mmg.DBG();
    std::cerr << "Graph loaded" << std::endl;
    std::unordered_map<std::string, Contig> batch;
    nano::ReadsAlignerGA reads_aligner(mg);
    nano::SGraphBuilder sgraph_builder(mg, unique_edges, true);
    int cnt = 0;
    for(StringContig scontig : reader) {
        Contig contig = scontig.makeContig();
        batch[contig.id] = contig;
        if (batch.size() > BATCH_SIZE) {
            std::unordered_map<std::string, nano::GraphContig> alignments = reads_aligner.Align(batch, graph, dir, ++cnt);
            for (auto const &[key, val]: alignments) {
                val.PrintContig();
                for (auto const edge: val.path) {
                    cout << edge << " " << edge.substr(0, edge.size() - 1) << " "
                        << mg.edges.count(atoi(edge.substr(0, edge.size() - 1).c_str())) << std::endl;
                }
            }
            sgraph_builder.LoadAlignments(alignments);
            batch.clear();
        }
    }
    std::unordered_map<std::string, nano::GraphContig> alignments = reads_aligner.Align(batch, graph, dir, ++cnt);
    for (auto const &[key, val]: alignments) {
        val.PrintContig();
    }
    sgraph_builder.LoadAlignments(alignments);
    batch.clear();

    sgraph_builder.PrintSgraph();
//    GraphSimplificator graph_simplificator(sgraph_builder.GetSGraph());
//    multigraph::MultiGraph smg = graph_simplificator.Simplify(mg);
    return 0;
}

int main(int argc, char **argv) {
    CLParser parser({"output-dir=", "unique-edges=", "graph=", "reads=", "threads=8"}, {},{"o=output-dir", "u=unique-edges", "g=graph", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    StringContig::homopolymer_compressing = true;
    const size_t threads = std::stoi(parser.getValue("threads"));
    const std::experimental::filesystem::path unique_edges(parser.getValue("unique-edges"));
    const std::experimental::filesystem::path ont_reads(parser.getValue("reads"));
    const std::experimental::filesystem::path graph(parser.getValue("graph"));
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    std::cerr << "Run Nano" << std::endl;
    runNano(graph, ont_reads, unique_edges, threads, dir);
    return 0;
}