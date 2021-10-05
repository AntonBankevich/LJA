//
// Created by anton on 8/3/20.
//
#pragma once

#include "minimizer_selection.hpp"
#include "dbg_disjointigs.hpp"
#include "sparse_dbg.hpp"
#include "common/rolling_hash.hpp"
#include "sequences/sequence.hpp"
#include "common/bloom_filter.hpp"
#include "common/output_utils.hpp"
#include "common/logging.hpp"
#include "common/simple_computation.hpp"
#include "common/omp_utils.hpp"
#include <wait.h>

std::vector<hashing::htype> findJunctions(logging::Logger & logger, const std::vector<Sequence>& disjointigs,
                                 const hashing::RollingHash &hasher, size_t threads);
SparseDBG constructDBG(logging::Logger & logger, const std::vector<hashing::htype> &vertices,
                       const std::vector<Sequence> &disjointigs, const hashing::RollingHash &hasher, size_t threads);
SparseDBG DBGPipeline(logging::Logger & logger, const hashing::RollingHash &hasher, size_t w, const io::Library &lib,
                                const std::experimental::filesystem::path &dir, size_t threads,
                                const std::string& disjointigs_file = "none", const std::string &vertices_file = "none");