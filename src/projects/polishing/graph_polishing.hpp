#pragma once

#include <sequences/sequence.hpp>
#include <unordered_map>
#include <common/disjoint_sets.hpp>
#include <common/pipeline_tools.hpp>

std::unordered_map<std::string, std::experimental::filesystem::path> RunGraphPolishing(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
        const io::Library &graph_gfa,
        const io::Library &corrected_edges, bool debug);

class GraphPolishingPhase : public Stage {
public:
    GraphPolishingPhase() : Stage(AlgorithmParameters(
            {},
            {}, ""), {"graph", "corrected_edges"}, {}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        return RunGraphPolishing(logger, threads, dir, input.at("graph"), input.at("corrected_edges"), debug);
    }
};