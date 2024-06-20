#pragma once

#include "diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
#include <experimental/filesystem>
ag::AlignedReadStorage<dbg::DBGTraits> MultCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, const std::experimental::filesystem::path &dir,
                                            dbg::DBGAlignedReadStorage &reads_storage, size_t unique_threshold, double initial_rel_coverage, bool diploid,
                                            bool debug);