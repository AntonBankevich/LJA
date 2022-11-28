#pragma once
#include "uncompressed_output.hpp"
#include "homopolish.hpp"
#include "perfect_alignment.hpp"

std::vector<std::experimental::filesystem::path> PolishingPhase(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &output_dir,
        const std::experimental::filesystem::path &gfa_file,
        const io::Library &corrected_reads,
        const io::Library &reads, size_t dicompress, size_t min_alignment, bool skip, bool debug) {
    logger.info() << "Performing polishing and homopolymer uncompression" << std::endl;
    std::function<void()> ic_task = [&logger, threads, &output_dir, debug, &gfa_file, &corrected_reads, &reads, dicompress, min_alignment, &dir] {
        io::SeqReader reader(corrected_reads);
        multigraph::MultiGraph vertex_graph;
        vertex_graph.LoadGFA(gfa_file, true);
        multigraph::MultiGraph edge_graph = vertex_graph.DBG();
        std::vector<Contig> contigs = edge_graph.getEdges(false);
        auto res = PrintAlignments(logger, threads, contigs, reader.begin(), reader.end(), min_alignment, dir);
        std::vector<Contig> uncompressed = Polish(logger, threads, contigs, res.first, reads, dicompress);
        std::vector<Contig> assembly = printUncompressedResults(logger, threads, edge_graph, uncompressed, output_dir, debug);
        logger.info() << "Printing final assembly to " << (output_dir / "assembly.fasta") << std::endl;
        std::ofstream os_cut;
        os_cut.open(output_dir / "assembly.fasta");
        for(Contig &contig : assembly) {
            if(contig.size() > 1500)
                os_cut << ">" << contig.id << "\n" << contig.seq << "\n";
        }
        os_cut.close();
    };
    if(!skip)
        runInFork(ic_task);
    return {output_dir / "assembly.fasta", output_dir / "mdbg.gfa"};
}
