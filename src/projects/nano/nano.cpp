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
#include <nano/GraphSimplificator.h>
#include <nano/VertexExtender.h>
#include <lja/multi_graph.hpp>

using namespace dbg;

int runNano(const std::experimental::filesystem::path &graph,
            const  std::experimental::filesystem::path &ont_reads,
            const  std::experimental::filesystem::path &unique_edges,
            const size_t threads,
            const std::experimental::filesystem::path &dir){
    const unsigned int BATCH_SIZE = 10000;
    io::SeqReader reader(ont_reads);
    std::cerr << "Seq reader loaded" << std::endl;
    multigraph::MultiGraph mmg;
    mmg.LoadGFA(graph, true);
    std::cerr << "GFA loaded" << std::endl;
    multigraph::MultiGraph mg = mmg.DBG();
    const std::experimental::filesystem::path &input_gfa = dir / "input.gfa";
    mg.printEdgeGFA(input_gfa);
    std::cerr << "Graph loaded" << std::endl;

    nano::VertexExtender initial_graph_simplificator;
    std::unordered_map<int, int> new_edges_map = initial_graph_simplificator.ExtendVertices(mg);
    //std::unordered_map<int, int> new_edges_map;
    const std::experimental::filesystem::path &input_extended_gfa = dir / "input_extended.gfa";
    mg.printEdgeGFA(input_extended_gfa);
    std::cerr << "Graph simplified " << input_extended_gfa << std::endl;

    std::unordered_map<std::string, Contig> batch;
    nano::ReadsAlignerGA reads_aligner(mg);
    nano::SGraphBuilder sgraph_builder(mg, unique_edges, new_edges_map, false);
    int cnt = 0;
    StringContig::homopolymer_compressing = true;
    StringContig::min_dimer_to_compress = 32;
    StringContig::max_dimer_size = 32;
    for(StringContig scontig : reader) {
        Contig contig = scontig.makeContig();
        batch[contig.id] = contig;
        if (batch.size() > BATCH_SIZE) {
            std::unordered_map<std::string, nano::GraphContig> alignments =
                            reads_aligner.Align(batch, input_extended_gfa, dir, threads, ++cnt);
            sgraph_builder.LoadAlignments(alignments, threads);
            batch.clear();
        }
    }
    std::unordered_map<std::string, nano::GraphContig> alignments = reads_aligner.Align(batch, input_extended_gfa, dir, threads, ++cnt);
    sgraph_builder.LoadAlignments(alignments, threads);
    batch.clear();

    sgraph_builder.SaveSGraph(dir / "scaffold_graph.tsv");
    //sgraph_builder.LoadSGraphEdges(dir / "scaffold_graph.tsv");
    nano::GraphSimplificator graph_simplificator(sgraph_builder.GetSGraph());
    graph_simplificator.ResolveWithMajor(mg);
    graph_simplificator.Simplify(mg, dir);
    const std::experimental::filesystem::path &final_gfa = dir / "final.gfa";
    mg.printEdgeGFA(final_gfa);
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