#pragma once

#include "graph_polishing.hpp"
#include "uncompressed_output.hpp"
#include "homopolish.hpp"
#include "perfect_alignment.hpp"
#include <common/pipeline_tools.hpp>
#include <unordered_map>


std::unordered_map<std::string, std::experimental::filesystem::path> RunPolishing(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &gfa_file,
        const io::Library &corrected_reads,
        const io::Library &reads, size_t min_alignment, bool debug);

class PolishingPhase : public Stage {
public:
    PolishingPhase() : Stage(AlgorithmParameters(
            {"min_alignment=5001"},
            {}, ""), {"graph", "corrected_reads", "reads"}, {"assembly", "graph"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override;
};
