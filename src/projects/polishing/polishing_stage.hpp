#pragma once

#include "graph_polishing.hpp"
#include "uncompressed_output.hpp"
#include "homopolish.hpp"
#include "perfect_alignment.hpp"
#include <common/pipeline_tools.hpp>
#include <unordered_map>


//    Exactly one of corrected_reads and read_alignments is used: with a non-empty read_alignments (a file in
//    PrintAlignments format, e.g. from spg::PrintStorageAlignments) reads are not realigned to the graph.
std::unordered_map<std::string, std::experimental::filesystem::path> RunPolishing(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &gfa_file,
        const io::Library &corrected_reads, const std::experimental::filesystem::path &read_alignments,
        const io::Library &reads, size_t min_alignment, bool debug);

//    With realign, compressed reads (input "corrected_reads") are realigned to the graph. Otherwise ready alignments
//    to the graph vertices (input "read_alignments") are used.
class PolishingPhase : public Stage {
private:
    bool realign;
public:
    explicit PolishingPhase(bool realign = true) : Stage(AlgorithmParameters(
            {"min-alignment=500"},
            {}, ""), {"graph", realign ? "corrected_reads" : "read_alignments", "reads"}, {"assembly", "graph"}),
            realign(realign) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override;
};
