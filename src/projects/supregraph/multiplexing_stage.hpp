#pragma once
#include "multiplexer.hpp"
#include "read_storage.hpp"
#include "converter.hpp"
#include "unique_vertex_storage.hpp"
#include "decision_rules.hpp"
#include <dbg/dbg_construction.hpp>
#include <dbg/graph_algorithms.hpp>
#include <common/pipeline_tools.hpp>
#include <dbg/dbg_read_alignment_storage.hpp>
#include <error_correction/multiplicity_estimation.hpp>
#include <dbg/aln_reads_reader.hpp>

namespace spg {

void CleanSupregraph(ag::AssemblyGraph &dbg);
UniqueClassificator ConstructUnique(logging::Logger &logger, size_t threads, const io::Library &reads_files,
                                    dbg::SparseDBG &dbg, const std::experimental::filesystem::path &dir);

std::unordered_map<std::string, std::experimental::filesystem::path>
RunMultiplexing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, size_t k,
                size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                const io::Library &extra_reads_files, const io::Library &paths, bool debug);

class SupreGraphPhase : public Stage {
public:
    SupreGraphPhase();

protected:
    std::unordered_map<std::string, std::experimental::filesystem::path>
    innerRun(logging::Logger &logger, size_t threads,
             const std::experimental::filesystem::path &dir, bool debug,
             const AlgorithmParameterValues &parameterValues,
             const std::unordered_map<std::string, io::Library> &input) override;
};
}