#pragma once
#include "supregraph.hpp"
#include "multiplexer.hpp"
#include "read_storage.hpp"
#include "converter.hpp"
#include "unique_vertex_storage.hpp"
#include "decision_rules.hpp"
#include "path_tracker.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "assembly_graph/long_vertex_multiplicity.hpp"
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

    struct SupregraphFromDBG {
        Supregraph spg;
        dbg::DBGAlignedReadStorage dbg_storage;
    };

//    Thin wrapper for every place that builds a Supregraph out of an already-loaded DBG + read storage:
//    both are taken by move (the caller is done with them in their pre-conversion form) and handed to
//    Supregraph::FromDBG. CoverageTracker and UniqueVertexStorage are deliberately NOT built here --
//    CoverageTracker's fillFromStorage must run only after every edge has been converted (constructing it
//    before that lets it react to in-loop read reroutes against sample windows that don't exist yet), and
//    UniqueVertexStorage's uniqueness threshold/criteria is a per-caller choice. Build both from the
//    returned spg/dbg_storage.
    SupregraphFromDBG
    ConstructSupregraphFromDBG(logging::Logger &logger, size_t threads, size_t k, ag::AssemblyGraph &&graph,
                               dbg::DBGAlignedReadStorage &&dbg_storage);

    struct SupregraphFromRawInput {
        Supregraph spg;
        dbg::DBGAlignedReadStorage dbg_storage;
//        Declared last so it destructs first: struct members destruct in reverse declaration order (same
//        rule as local variables), and path_storage is a listener attached to spg, so it must detach
//        before spg's own ~AssemblyGraph() (which VERIFYs its listener list is empty) runs.
        ag::DLLAlignmentStorage path_storage;
    };

//    Shared top-level entry point for every place that builds a Supregraph straight from raw pipeline
//    input (a graph/reads library, not an already-loaded graph): loads the DBG, loads reads, attaches
//    path_storage to the DBG and runs PrepareDLLPathTracker on it (debug only -- needs dbg's own k-mer
//    hashing, so this must happen before disableHashing()/conversion), then hands off to
//    ConstructSupregraphFromDBG. path_storage is attached to the DBG rather than the returned Supregraph
//    for the same reason as ConstructSupregraphFromDBG's internals: it must observe the edgeToSupreVertex
//    conversions too, and a listener's fire pointer (plus the graph's own registration of it) stays
//    correct across the moves involved (see fire_listeners.cpp).
    SupregraphFromRawInput
    ConstructSupregraphFromRawInput(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                                    size_t k, size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                                    const io::Library &extra_reads_files, const io::Library &paths, bool debug);

//    Ends with multiplicity estimation of long vertices of the final graph (see ag::LongVertexMultiplicities):
//    dir/multiplicities.tsv and Mult labels in supregraph_final.dot.
std::unordered_map<std::string, std::experimental::filesystem::path>
RunMultiplexing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, size_t k,
                size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                const io::Library &extra_reads_files, const io::Library &paths,
                const ag::MultiplicityParams &multiplicity_params, bool debug);

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