#include <common/logging.hpp>
#include <sequences/sequence.hpp>
#include "dbg_read_alignment_storage.hpp"

using namespace dbg;

void DBGAlignedReadStorage::checkCoverage(const SparseDBG &dbg) const {
    VERIFY_MSG(coverageTracker != nullptr, "Error: checking coverage using read storage that does not contribute to coverage");
    std::unordered_map<ag::ConstEdgeId, size_t> map;
    for (const Edge &edge: dbg.edges()) {
        map[edge.getId()] = 0;
    }
    for (const ag::AlignedRead &read: *this) {
        if (!read.valid())
            continue;
        for (ag::PathPosition pp = read.getPath().firstPosition(); pp != read.getPath().lastPosition(); ++pp) {
            map[pp.nextEdge().getId()] += read.getPath().getSegment(pp).size();
            map[pp.nextEdge().rc().getId()] += read.getPath().getSegment(pp).size();
        }
    }
    for (const Edge &edge: dbg.edges()) {
        VERIFY_MSG(edge.intCov() == map[edge.getId()], "Coverage check failed");
    }
}

void DBGAlignedReadStorage::logReads(size_t threads, const std::experimental::filesystem::path &path) {
    readLogger = new ag::ReadLogger(*this, threads, path);
}

void DBGAlignedReadStorage::logGraph(SparseDBG &dbg, std::ostream &os) {
    graphLogger = new ag::LoggingListener(dbg, os);
}

void DBGAlignedReadStorage::stopTrackSuffixes() {
    delete suffixes;
    suffixes = nullptr;
}

void DBGAlignedReadStorage::stopTrackCoverage() {
    delete coverageTracker;
    coverageTracker = nullptr;
}

void DBGAlignedReadStorage::trackSuffixes(logging::Logger &logger, size_t threads, SparseDBG &dbg, size_t _min_len,
                                          size_t _max_len) {
    stopTrackSuffixes();
    logger.info() << "Storing suffixes of read paths of length up to " << _max_len << std::endl;
    suffixes = new ag::SuffixTracker(*this, dbg, _min_len, _max_len);
    suffixes->fillFromStorage(logger, threads);
}

DBGAlignedReadStorage::~DBGAlignedReadStorage() {
    delete suffixes;
    delete readLogger;
    delete graphLogger;
    delete coverageTracker;
}

DBGAlignedReadStorage
DBGAlignedReadStorage::Load(logging::Logger &logger, size_t threads, const io::Library &lib,
                            SparseDBG &dbg, bool track_cov) {
    std::vector<ag::AlignedRead> reads;
    IdIndex<Vertex> index(dbg.vertices().begin(), dbg.vertices().end());
    for(const auto &path : lib) {
        std::ifstream is(path);
        std::vector<ag::AlignedRead> tmp = ag::AlignedReadStorage::LoadReadAlignments(is, index);
        is.close();
        for(ag::AlignedRead &read : tmp) {
            reads.emplace_back(std::move(read));
        }
    }
    return {logger, threads, dbg, std::move(reads), track_cov};
}

DBGAlignedReadStorage
DBGAlignedReadStorage::Load(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path,
                            SparseDBG &dbg, bool track_cov) {
    logger.info() << "Loading read alignments from " << path << std::endl;
    DBGAlignedReadStorage res(logger, threads, dbg,
            ag::AlignedReadStorage::Load(logger, threads, path, dbg), track_cov);
    logger.info() << "Finished loading read alignments from " << path << std::endl;
    return res;
}

DBGAlignedReadStorage
DBGAlignedReadStorage::Load(logging::Logger &logger, size_t threads, std::istream &is, SparseDBG &dbg,
                            bool track_cov) {
    return {logger, threads, dbg,
            ag::AlignedReadStorage::Load(logger, threads, is, dbg), track_cov};
}
