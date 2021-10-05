#pragma once
#include "paths.hpp"
#include "sparse_dbg.hpp"
using namespace dbg;

template<class Iterator>
void FillSparseDBGEdges(SparseDBG &sdbg, Iterator begin, Iterator end, logging::Logger &logger, size_t threads,
                        const size_t min_read_size) {
    typedef typename Iterator::value_type ContigType;
    logger.trace() << "Starting to fill edges" << std::endl;
    std::function<void(size_t, ContigType &)> task = [&sdbg, min_read_size](size_t pos, ContigType &contig) {
        Sequence seq = contig.makeSequence();
        if (seq.size() >= min_read_size)
            sdbg.processRead(seq);
    };
    processRecords(begin, end, logger, threads, task);
    logger.trace() << "Sparse graph edges filled." << std::endl;
}

template<class Iterator>
void RefillSparseDBGEdges(SparseDBG &sdbg, Iterator begin, Iterator end, logging::Logger &logger, size_t threads) {
    logger.trace() << "Starting to fill edges" << std::endl;
    std::function<void(size_t, std::pair<Vertex *, Sequence> &)> task = [&sdbg](size_t pos, std::pair<Vertex *, Sequence> &contig) {
        sdbg.processEdge(*contig.first, contig.second);
    };
    processObjects(begin, end, logger, threads, task);
    logger.trace() << "Sparse graph edges filled." << std::endl;
}

SparseDBG LoadDBGFromFasta(const io::Library &lib, hashing::RollingHash &hasher, logging::Logger &logger, size_t threads);

template<class Iterator>
void fillCoverage(SparseDBG &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                  const hashing::RollingHash &hasher, size_t min_read_size);

SparseDBG constructSparseDBGFromReads(logging::Logger & logger, const io::Library &reads_file, size_t threads,
                                      const hashing::RollingHash &hasher, const std::vector<hashing::htype> &hash_list, size_t w);

void tieTips(logging::Logger &logger, SparseDBG &sdbg, size_t w, size_t threads);

void UpdateVertexTips(Vertex &rec, ParallelRecordCollector<Vertex *> &queue);

void findTips(logging::Logger &logger, SparseDBG &sdbg, size_t threads);

void mergeLoop(Path path);

void MergeEdge(SparseDBG &sdbg, Vertex &start, Edge &edge);

void mergeLinearPaths(logging::Logger & logger, SparseDBG &sdbg, size_t threads);

void mergeCyclicPaths(logging::Logger & logger, SparseDBG &sdbg, size_t threads);

void mergeAll(logging::Logger & logger, SparseDBG &sdbg, size_t threads);

void CalculateCoverage(const std::experimental::filesystem::path &dir, const hashing::RollingHash &hasher, const size_t w,
                       const io::Library &lib, size_t threads, logging::Logger &logger, SparseDBG &dbg);

std::experimental::filesystem::path alignLib(logging::Logger &logger, SparseDBG &dbg, const io::Library &align_lib,
        const hashing::RollingHash &hasher, const size_t w, const std::experimental::filesystem::path &dir, size_t threads);