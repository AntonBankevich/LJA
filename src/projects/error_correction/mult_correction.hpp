#pragma once

#include "diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/compact_path.hpp"
#include "dbg/graph_alignment_storage.hpp"
#include <experimental/filesystem>
RecordStorage MultCorrect(dbg::SparseDBG &dbg, logging::Logger &logger,
                 const std::experimental::filesystem::path &dir,
                 RecordStorage &reads_storage, size_t unique_threshold,
                 size_t threads, bool diploid, bool debug);