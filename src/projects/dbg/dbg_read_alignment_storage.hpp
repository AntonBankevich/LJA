#pragma once
#include "sparse_dbg.hpp"
#include "dbg_graph_aligner.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"

namespace dbg {
    using ag::GraphPath;
//    DBG only. intCov(e) is the total size of the segments of all valid stored reads on e plus their segments on
//    e.rc(): every read contributes its forward and its backward (k+1)-mers. So a self-rc edge (a match between the
//    forward and the rc genome) gets each read's contribution twice.
    class DBGAlignedReadStorage;
    class CoverageTracker : public ag::AlignedReadStorageListener, public DBGListener {
        friend class DBGAlignedReadStorage;
    private:
        ag::AlignedReadStorage *storage;
        static void addPath(const GraphPath &path, __int64_t mult = 1) {
//            Legacy paths are empty subreads of deleted vertices: no edges, and their vertex must not be touched.
            if (!path.valid() || path.isLegacy())
                return;
            for(ag::PathPosition pp = path.firstPosition(); pp != path.lastPosition(); ++pp) {
                __int64_t val = __int64_t(path.getSegment(pp).size()) * mult;
                Edge &edge = pp.nextEdge();
                edge.incCov(val);
                edge.rc().incCov(val);
            }
        }
        void fillFromStorage(logging::Logger &logger, size_t threads) {
            omp_set_num_threads(threads);
            logger.info() << "Filling edge coverages" << std::endl;
#pragma omp parallel for default(none) schedule(dynamic, 100)
            for(size_t i = 0; i < storage->size(); i++) {
                fireAddRead((*storage)[i]);
            }
            logger.info() << "Finished filling edge coverages" << std::endl;
        }
    public:
        explicit CoverageTracker(logging::Logger &logger, size_t threads, ag::AlignedReadStorage &storage, ag::ResolutionFire &graph) :
            ag::AlignedReadStorageListener(storage, "CoverageTracker"), DBGListener(graph, "CoverageTracker"), storage(&storage) {
            fillFromStorage(logger, threads);
        }
        void fireAddRead(const ag::AlignedRead &read) override {addPath(read.getPath());}
        void fireRerouteRead(ag::AlignedRead &read) override {addPath(read.getPath(), -1); addPath(read.getCorrected());}
        void fireInvalidateRead(ag::AlignedRead &read) override {addPath(read.getPath(), -1);}
        void fireAddEdge(Edge &e) override {e.setCov(0);}

        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) override {
            for(Edge &e : path.edges())
                new_edge.incCov(int64_t(e.intCov()));
        }
//        Adds the coverage of the reads whose path starts with the new edge e; returns whether some of them leave e.
//        Reads are already moved by AlignedReadStorageMaintenance. A read ending inside e is seen again from e.rc()
//        (by the mirrored call, or by this one for a self-rc e), so only e gets its part. A read leaving e is not, so
//        both strands get it. Reads entering e from elsewhere are the other direction of such reads.
        bool addReadsStartingWith(Edge &e);
        void fireMergeTipsToEdge(Edge &new_edge, Edge &, Edge &, const AlignmentForm &, const AlignmentForm &) override {
            addReadsStartingWith(new_edge);
        }
//        No read may cross a split point, so only reads on the last piece can leave it.
        void fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) override;
        bool checkConsistency(std::ostream &os) override;
        bool fireCheckConsistency(std::ostream &os) override {return checkConsistency(os);}
    };

//    Create a proper structure that can work with multiple libraries
    class DBGAlignedReadStorage : public ag::AlignedReadStorage {
    private:
        ag::SuffixTracker * suffixes = nullptr;
        dbg::CoverageTracker * coverageTracker = nullptr;
        ag::ReadLogger *readLogger = nullptr;
        ag::LoggingListener *graphLogger = nullptr;
    public:
        ag::AlignedReadStorage &getReads() {return *this;}
        const ag::AlignedReadStorage &getReads() const {return *this;}
        ag::SuffixTracker &getSuffixes() const {return *suffixes;}
        bool tracksSuffixes() const {return suffixes != nullptr;}


        DBGAlignedReadStorage(logging::Logger &logger, size_t threads, SparseDBG &dbg, ag::AlignedReadStorage reads, bool _track_cov = false) :
                ag::AlignedReadStorage(std::move(reads)) {
            if(_track_cov) {
                coverageTracker = new CoverageTracker(logger, threads, *this, dbg);
            }
        }
        DBGAlignedReadStorage(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                                       std::vector<ag::AlignedRead> read_list, bool _track_cov = false) :
                                       ag::AlignedReadStorage(logger, threads, dbg, std::move(read_list)) {
            if(_track_cov) {
                coverageTracker = new CoverageTracker(logger, threads, *this, dbg);
            }
        }

//        The trackers move with the storage: the source no longer owns them, and they are pointed at the new storage.
        DBGAlignedReadStorage(DBGAlignedReadStorage &&other) noexcept;
        DBGAlignedReadStorage& operator=(DBGAlignedReadStorage&&) = delete;

        void logReads(size_t threads, const std::experimental::filesystem::path& path);
        void logGraph(SparseDBG &dbg, std::ostream &os);
        void trackSuffixes(logging::Logger &logger, size_t threads, ag::AssemblyGraph &dbg, size_t _min_len, size_t _max_len);
        void stopTrackSuffixes();
        void stopTrackCoverage();
        void checkCoverage() const;
        static DBGAlignedReadStorage Load(logging::Logger &logger, size_t threads, const io::Library &lib, SparseDBG &dbg, bool track_cov);
        static DBGAlignedReadStorage Load(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path, SparseDBG &dbg,
                                          bool track_cov);
        static DBGAlignedReadStorage Load(logging::Logger &logger, size_t threads, std::istream &is, SparseDBG &dbg,
                                          bool track_cov);
        ~DBGAlignedReadStorage() override;
    };

    template<class I>
    std::vector<ag::AlignedRead> AlignReads(logging::Logger &logger, size_t threads, I begin, I end, SparseDBG &dbg, size_t w) {
        logger.info() << "Loading and aligning reads."  << std::endl;
        KmerIndex index(dbg);
        index.fillAnchors(logger, threads, dbg, w);
        OrderedRecordCollector<ag::AlignedRead> tmpReads(threads);
        std::function<void(size_t, StringContig &)> read_task = [&index, &tmpReads](size_t pos, StringContig &scontig) {
            Contig contig = scontig.makeContig();
            if (contig.truncSize() < index.minReadLen()) {
                tmpReads.emplace_back(pos, ag::AlignedRead(contig.getInnerId(), GraphPath()));
            } else {
                tmpReads.emplace_back(pos, ag::AlignedRead(contig.getInnerId(), index.align(contig.getSeq(), contig.getInnerId())));
            }
        };
        processRecords(begin, end, logger, threads, read_task);
        std::vector<ag::AlignedRead> res = tmpReads.collectOrdered(logger, threads);
        logger.info() << "Read alignment completed."  << std::endl;
        return std::move(res);
    }
}
