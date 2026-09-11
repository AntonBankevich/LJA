#pragma once

#include "multiplexer.hpp"
#include "multiplexing_stage.hpp"
#include "read_storage.hpp"
#include "decision_rules.hpp"
#include <common/pipeline_tools.hpp>
#include <dbg/dbg_read_alignment_storage.hpp>
#include <error_correction/precorrection.hpp>
#include <vector>

namespace spg {

//    Combines multiplexing and Precorrector-based error correction on a Supregraph: cores are
//    multiplexed in order of increasing vertex size (Multiplexer::core_queue is already sorted
//    this way) under a progressively loosened size threshold, with a Precorrector pass run
//    between threshold steps. Multiplexing untangles reads so later Precorrector passes see
//    better outgoing_read_count statistics, and correction cleans reads so later, larger cores
//    get better-supported resolution judgements.
//    Precorrector's isReliable/isSuspicious are both driven by the edge's outgoing_read_count
//    (suffix edges never track this field and are always treated as reliable/non-suspicious).
//    Threshold schedule: initial_core_length, +step, +step, ... up to (and including) max_core_length,
//    followed by one final unrestricted pass.
    std::unordered_map<std::string, std::experimental::filesystem::path>
    RunMultiplexingAndCorrection(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                    size_t k, size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                    const io::Library &extra_reads_files, const io::Library &paths,
                    size_t initial_core_length, size_t max_core_length, size_t core_length_step,
                    size_t reliable_read_count, size_t suspicious_read_count, bool debug);

    class MultiplexAndCorrectionPhase : public Stage {
    public:
        MultiplexAndCorrectionPhase();

    protected:
        std::unordered_map<std::string, std::experimental::filesystem::path>
        innerRun(logging::Logger &logger, size_t threads,
                 const std::experimental::filesystem::path &dir, bool debug,
                 const AlgorithmParameterValues &parameterValues,
                 const std::unordered_map<std::string, io::Library> &input) override;
    };
}
