//
// Created by anton on 8/3/20.
//
#pragma once

#include "graph_algorithms.hpp"
#include "sparse_dbg.hpp"
#include <common/bloom_filter.hpp>
#include <common/cl_parser.hpp>
#include <common/simple_computation.hpp>
#include <common/logging.hpp>

Sequence buildDisjointig(dbg::Path &path);
void processVertex(dbg::Vertex &rec, ParallelRecordCollector<Sequence> &res);
void prepareVertex(dbg::Vertex &vertex);
void extractLinearDisjointigs(dbg::SparseDBG &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger & logger, size_t threads);
void extractCircularDisjointigs(dbg::SparseDBG &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger & logger, size_t threads);
std::vector<Sequence> extractDisjointigs(logging::Logger & logger, dbg::SparseDBG &sdbg, size_t threads);
std::vector<Sequence> constructDisjointigs(const hashing::RollingHash &hasher, size_t w, const io::Library &reads_file,
                                           const std::vector<hashing::htype> & hash_list, size_t threads,
                                           logging::Logger & logger);