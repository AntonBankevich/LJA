#pragma once

#include "multi_graph.hpp"
#include "uncompressed_output.hpp"
#include "gap_closing.hpp"
#include "polishing/perfect_alignment.hpp"
#include "error_correction/mult_correction.hpp"
#include "error_correction/mitochondria_rescue.hpp"
#include "error_correction/manyk_correction.hpp"
#include "error_correction/precorrection.hpp"
#include "sequences/seqio.hpp"
#include "dbg/dbg_construction.hpp"
#include "common/rolling_hash.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"

#include <iostream>
#include <queue>
#include <wait.h>
#include <error_correction/dimer_correction.hpp>
#include <polishing/homopolish.hpp>
#include "common/rolling_hash.hpp"
#include <utility>
#include <ksw2/ksw_wrapper.hpp>

namespace pipeline {

struct LJAPipeline {

    std::vector<Contig> ref;
    size_t stage_num;
    LJAPipeline(io::Library ref_lib) {
        ref = io::SeqReader(ref_lib).readAllContigs();
        stage_num = 0;
    }

void PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const std::string &stage,
                                       SparseDBG &dbg, RecordStorage &readStorage, const io::Library &paths_lib, bool small);

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
AlternativeCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                             const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
                                             size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
                                             bool diploid, bool skip, bool debug, bool load);

std::vector<std::experimental::filesystem::path> NoCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                                                                           const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
                                                                                           size_t threads, size_t k, size_t w, bool skip, bool debug, bool load);
std::vector<std::experimental::filesystem::path> SecondPhase(
        logging::Logger &logger, const std::experimental::filesystem::path &dir,
        const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
        const io::Library &paths_lib, size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
        size_t unique_threshold, bool diploid, bool skip, bool debug, bool load);

    std::vector<std::experimental::filesystem::path> PolishingPhase(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &output_dir,
        const std::experimental::filesystem::path &gfa_file,
        const std::experimental::filesystem::path &corrected_reads,
        const io::Library &reads, size_t dicompress, size_t min_alignment, bool skip, bool debug);

std::vector<std::experimental::filesystem::path> MDBGPhase(
        logging::Logger &logger, size_t threads, size_t k, size_t kmdbg, size_t w, size_t unique_threshold, bool diploid,
        const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &graph_fasta,
        const std::experimental::filesystem::path &read_paths, bool skip, bool debug);

//std::vector<std::experimental::filesystem::path>
void TrioPreprocessingPhase(
            logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
            const io::Library &p_lib, const io::Library &m_lib,
            bool skip, bool debug);

std::experimental::filesystem::path TrioBinningPhase(
            logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
            std::experimental::filesystem::path pat_kmers, std::experimental::filesystem::path mat_kmers,
            std::experimental::filesystem::path contigs, bool skip, bool debug);

std::vector<std::experimental::filesystem::path> TrioSimplificationPhase(
            logging::Logger &logger, size_t threads,
            const std::experimental::filesystem::path &graph, const std::experimental::filesystem::path &binned_contigs,
            const std::experimental::filesystem::path &corrected_reads, const io::Library &reads_lib,
            const std::experimental::filesystem::path &dir, size_t saved_bridge_cutoff, bool skip, bool debug);

std::experimental::filesystem::path simplifyHaplo(logging::Logger &logger, size_t threads,
                                                      const std::experimental::filesystem::path &output_file,
                                                      const std::experimental::filesystem::path &diplo_graph,
                                                      const std::experimental::filesystem::path &haployak,
                                                      const char haplotype,
                                                      const std::experimental::filesystem::path &corrected_reads,
                                                      const io::Library &reads,
                                                      const std::experimental::filesystem::path &dir,
                                                      const size_t saved_bridge_cutoff);

};


}
