#pragma once

#include "dbg/graph_modification.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/compact_path.hpp"
#include "error_correction/tip_correction.hpp"
#include <sequences/sequence.hpp>
#include <random>

class GapCloser {
private:
    size_t min_overlap;
    size_t max_overlap;
    size_t smallK;
    double allowed_divergence;

    struct OverlapRecord {
        OverlapRecord(size_t from, size_t to, size_t matchSizeFrom, size_t matchSizeTo) : from(from), to(to),
                                                                                          match_size_from(
                                                                                                  matchSizeFrom),
                                                                                          match_size_to(matchSizeTo) {}

        size_t from;
        size_t to;
        size_t match_size_from;
        size_t match_size_to;
    };
public:
    GapCloser(size_t min_overlap, size_t max_overlap, size_t smallK, double allowed_divergence) :
            min_overlap(min_overlap), max_overlap(max_overlap), smallK(smallK), allowed_divergence(allowed_divergence){}
    bool HasInnerDuplications(const Sequence &seq, const hashing::RollingHash &hasher);
    std::vector<Connection> GapPatches(logging::Logger &logger, dbg::SparseDBG &dbg, size_t threads);
};

void MarkUnreliableTips(dbg::SparseDBG &dbg, const std::vector<Connection> &patches);

void GapCloserPipeline(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                       const std::vector<RecordStorage *> &storges);
