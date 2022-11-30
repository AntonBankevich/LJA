#pragma once
#include "uncompressed_output.hpp"
#include "homopolish.hpp"
#include "perfect_alignment.hpp"

std::unordered_map<std::string, std::experimental::filesystem::path> RunPolishing(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &gfa_file,
        const io::Library &corrected_reads,
        const io::Library &reads, size_t min_alignment, bool debug) {
    logger.info() << "Performing polishing and homopolymer uncompression" << std::endl;
    size_t dicompress = StringContig::max_dimer_size / 2;
    io::SeqReader reader(corrected_reads);
    multigraph::MultiGraph vertex_graph;
    vertex_graph.LoadGFA(gfa_file, true);
    multigraph::MultiGraph edge_graph = vertex_graph.DBG();
    std::vector<Contig> contigs = edge_graph.getEdges(false);
    auto res = PrintAlignments(logger, threads, contigs, reader.begin(), reader.end(), min_alignment, dir);
    std::vector<Contig> uncompressed = Polish(logger, threads, contigs, res.first, reads, dicompress);
    std::vector<Contig> assembly = printUncompressedResults(logger, threads, edge_graph, uncompressed, dir, debug);
    logger.info() << "Printing final assembly to " << (dir / "assembly.fasta") << std::endl;
    std::ofstream os_cut;
    os_cut.open(dir / "assembly.fasta");
    for(Contig &contig : assembly) {
        if(contig.size() > 1500)
            os_cut << ">" << contig.id << "\n" << contig.seq << "\n";
    }
    os_cut.close();
    return {{"assembly", dir / "assembly.fasta"}, {"graph", dir / "mdbg.gfa"}};
}

class PolishingPhase : public Stage {
public:
    PolishingPhase() : Stage(AlgorithmParameters(
            {"min_alignment=5001"},
            {}, ""), {"graph", "corrected_reads", "reads"}, {"assembly", "graph"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        logger.info() << "Started kmer counting phase for trio binning\n";
        size_t min_alignment = std::stoull(parameterValues.getValue("min_alignment"));
        return RunPolishing(logger, threads, dir, input.find("graph")->second.front(), input.find("corrected_reads")->second,
                            input.find("reads")->second, min_alignment, debug);
    }
};
