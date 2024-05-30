#pragma once

#include "diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"
#include "dbg/sparse_dbg.hpp"
#include "assembly_graph/compact_path.hpp"
#include "dbg/graph_alignment_storage.hpp"
#include <experimental/filesystem>
dbg::ReadAlignmentStorage MultCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, const std::experimental::filesystem::path &dir,
                          dbg::ReadAlignmentStorage &reads_storage, size_t unique_threshold, double initial_rel_coverage, bool diploid,
                          bool debug);