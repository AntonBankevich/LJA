#include <common/logging.hpp>
#include <sequences/sequence.hpp>
#include "dbg_read_alignment_storage.hpp"

using namespace dbg;

void DBGAlignedReadStorage::checkCoverage(const SparseDBG &dbg) const {
    VERIFY_MSG(coverageTracker != nullptr, "Error: checking coverage using read storage that does not contribute to coverage");
    std::unordered_map<ConstEdgeId, size_t> map;
    for (const Edge &edge: dbg.edges()) {
        map[edge.getId()] = 0;
    }
    for (const ag::AlignedRead<DBGTraits> &read: *this) {
        if (!read.valid())
            continue;
        for (ag::PathPosition<DBGTraits> pp = read.getPath().firstPosition(); pp != read.getPath().lastPosition(); ++pp) {
            map[pp.nextEdge().getId()] += read.getPath().getSegment(pp).size();
            map[pp.nextEdge().rc().getId()] += read.getPath().getSegment(pp).size();
        }
    }
    for (const Edge &edge: dbg.edges()) {
        VERIFY_MSG(edge.intCov() == map[edge.getId()], "Coverage check failed");
    }
}

void DBGAlignedReadStorage::logReads(size_t threads, const std::experimental::filesystem::path &path) {
    readLogger = new ag::ReadLogger<DBGTraits>(*this, threads, path);
}

void DBGAlignedReadStorage::stopTrackSuffixes() {
    delete suffixes;
    suffixes = nullptr;
}

void DBGAlignedReadStorage::trackSuffixes(logging::Logger &logger, size_t threads, SparseDBG &dbg, size_t _min_len,
                                          size_t _max_len) {
    stopTrackSuffixes();
    logger.info() << "Storing suffixes of read paths of length up to " << _max_len << std::endl;
    suffixes = new ag::SuffixTracker<DBGTraits>(*this, dbg, _min_len, _max_len);
    suffixes->fillFromStorage(logger, threads);
}

DBGAlignedReadStorage::~DBGAlignedReadStorage() {
    delete suffixes;
    delete readLogger;
    delete coverageTracker;
}

DBGAlignedReadStorage
DBGAlignedReadStorage::Load(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path,
                            SparseDBG &dbg, const IdIndex<Vertex> &index, bool track_cov) {
    return {logger, threads, dbg,
            ag::AlignedReadStorage<DBGTraits>::Load(logger, threads, path, dbg, index), track_cov};
}

DBGAlignedReadStorage
DBGAlignedReadStorage::Load(logging::Logger &logger, size_t threads, std::istream &is, SparseDBG &dbg,
                            const IdIndex<Vertex> &index, bool track_cov) {
    return {logger, threads, dbg,
            ag::AlignedReadStorage<DBGTraits>::Load(logger, threads, is, dbg, index), track_cov};
}
