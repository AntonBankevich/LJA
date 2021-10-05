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

Sequence buildDisjointig(Path &path);
void processVertex(Vertex &rec, ParallelRecordCollector<Sequence> &res);
void prepareVertex(Vertex &vertex);
void extractLinearDisjointigs(SparseDBG &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger & logger, size_t threads);
void extractCircularDisjointigs(SparseDBG &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger & logger, size_t threads);
std::vector<Sequence> extractDisjointigs(logging::Logger & logger, SparseDBG &sdbg, size_t threads);
std::vector<Sequence> constructDisjointigs(const hashing::RollingHash &hasher, size_t w, const io::Library &reads_file,
                                           const std::vector<hashing::htype> & hash_list, size_t threads,
                                           logging::Logger & logger);