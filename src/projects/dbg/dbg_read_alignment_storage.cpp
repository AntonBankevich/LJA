#include <common/logging.hpp>
#include <sequences/sequence.hpp>
#include "dbg_read_alignment_storage.hpp"
#include "aln_reads_reader.hpp"

using namespace dbg;

std::vector<ag::AlignedRead> dbg::AlignReads(logging::Logger &logger, size_t threads, const io::Library &lib,
                                             SparseDBG &dbg, size_t w) {
    logger.info() << "Loading and aligning reads." << std::endl;
    KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, w);
    std::vector<ag::AlignedRead> res;
    ReadLoader(logger, threads).process(lib, [&index, &res](size_t pos, Contig &contig) {
        const std::string &id = contig.getInnerId();
        res[pos] = contig.truncSize() < index.minReadLen() ? ag::AlignedRead(id, GraphPath()) :
                   ag::AlignedRead(id, index.align(contig.getSeq(), id));
    }, [&res](size_t total) {res.resize(total);});
    logger.info() << "Read alignment completed." << std::endl;
    return res;
}

bool CoverageTracker::addReadsStartingWith(Edge &e) {
    bool leaves = false;
    for (ag::AlignedReadDirection dir: storage->getOutgoingReads(e.getStart())) {
//        e.getStart() may branch, so its bucket may hold reads leaving on other edges. The edges e replaces are still
//        corporeal, so dir.frontEdge() would resolve to them; match by code as the graph will afterwards.
        if (e.getCode().empty() || dir.getFSplits()[0] == e.getCode()[0]) {
            int64_t len(e.truncSize() - dir.leftCut());
            if (dir.getFSplits().size() <= e.getCode().size()) {
                len -= dir.rightCut();
            } else {
                leaves = true;
                e.rc().incCov(len);
            }
            e.incCov(len);
        }
    }
    return leaves;
}

void CoverageTracker::fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) {
    for (Edge &e : split.edges()) {
        bool leaves = addReadsStartingWith(e);
        VERIFY(!leaves || e == split.backEdge());
    }
}

bool CoverageTracker::checkConsistency(std::ostream &os) {
    std::unordered_map<ag::ConstEdgeId, size_t> expected;
    for (const ag::AlignedRead &read: *storage) {
        if (!read.valid() || read.getPath().isLegacy())
            continue;
        for (ag::PathPosition pp = read.getPath().firstPosition(); pp != read.getPath().lastPosition(); ++pp) {
            size_t len = read.getPath().getSegment(pp).size();
            expected[pp.nextEdge().getId()] += len;
            expected[pp.nextEdge().rc().getId()] += len;
        }
    }
    bool ok = true;
    for (const Edge &edge: ag::ResolutionListener::getGraph().edges()) {
        auto it = expected.find(edge.getId());
        size_t exp = it == expected.end() ? 0 : it->second;
        if (edge.intCov() != exp) {
            os << "CoverageTracker: edge " << edge.getId() << " has coverage " << edge.intCov()
               << " but its reads give " << exp << std::endl;
            ok = false;
        }
    }
    return ok;
}

void DBGAlignedReadStorage::checkCoverage() const {
    VERIFY_MSG(coverageTracker != nullptr, "Error: checking coverage using read storage that does not contribute to coverage");
    VERIFY_MSG(coverageTracker->checkConsistency(std::cerr), "Coverage check failed");
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

void DBGAlignedReadStorage::trackSuffixes(logging::Logger &logger, size_t threads, ag::AssemblyGraph &dbg, size_t _min_len,
                                          size_t _max_len) {
    stopTrackSuffixes();
    logger.info() << "Storing suffixes of read paths of length up to " << _max_len << std::endl;
    suffixes = new ag::SuffixTracker(*this, dbg, _min_len, _max_len);
    suffixes->fillFromStorage(logger, threads);
}

DBGAlignedReadStorage::DBGAlignedReadStorage(DBGAlignedReadStorage &&other) noexcept :
        ag::AlignedReadStorage(std::move(other)), suffixes(other.suffixes), coverageTracker(other.coverageTracker),
        readLogger(other.readLogger), graphLogger(other.graphLogger) {
    other.suffixes = nullptr;
    other.coverageTracker = nullptr;
    other.readLogger = nullptr;
    other.graphLogger = nullptr;
    if (suffixes != nullptr)
        suffixes->storage = this;
    if (coverageTracker != nullptr)
        coverageTracker->storage = this;
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
