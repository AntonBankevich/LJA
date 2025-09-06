#pragma once
#include "sparse_dbg.hpp"
#include "dbg_graph_aligner.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"

namespace dbg {

    class CoverageTracker : public ag::AlignedReadStorageListener<DBGTraits>, public ag::ResolutionListener<DBGTraits> {
    private:
        ag::AlignedReadStorage<DBGTraits> *storage;
        static void addPath(const GraphPath &path, __int64_t mult = 1) {
            for(ag::PathPosition<DBGTraits> pp = path.firstPosition(); pp != path.lastPosition(); ++pp) {
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
        explicit CoverageTracker(logging::Logger &logger, size_t threads, ag::AlignedReadStorage<DBGTraits> &storage, ag::ResolutionFire<DBGTraits> &graph) :
            ag::AlignedReadStorageListener<DBGTraits>(storage, "CoverageTracker"), ag::ResolutionListener<DBGTraits>(graph, "CoverageTracker"), storage(&storage) {
            fillFromStorage(logger, threads);
        }
        void fireAddRead(const ag::AlignedRead<DBGTraits> &read) override {addPath(read.getPath());}
        void fireRerouteRead(ag::AlignedRead<DBGTraits> &read) override {addPath(read.getPath(), -1); addPath(read.getCorrected());}
        void fireInvalidateRead(ag::AlignedRead<DBGTraits> &read) override {addPath(read.getPath(), -1);}
        void fireAddEdge(Edge &e) override {e.setCov(0);}

        void fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) override {
            for(EdgeId eid : path)
                new_edge.incCov(int64_t(eid->intCov()));
        }
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &, const AlignmentForm &) override {
    //        The same call ran for rc event will take care of coverage from left
            for(ag::AlignedReadDirection<DBGTraits> dir : storage->getOutgoingReads(right)) {
                int64_t len(new_edge.truncSize() - dir.leftCut());
                if(dir.getFSplits().size() <= left.getCode().size() + right.getCode().size())
                    len -= dir.rightCut();
                new_edge.incCov(len);
                VERIFY(len > 0);
                new_edge.rc().incCov(len);
            }
        }
        void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) override {
            for(EdgeId eid : split) {
                for (ag::AlignedReadDirection<DBGTraits> dir: storage->getOutgoingReads(*eid)) {
                    int64_t len(eid->truncSize() - dir.leftCut());
                    if(dir.getFSplits().size() <= eid->getCode().size())
                        len -= dir.rightCut();
                    else {
                        VERIFY(eid == split.back());
                        eid->rc().incCov(len);
                    }
                    eid->incCov(len);
                }
            }
        }
        void fireAddSupreVertex(Vertex &v, Edge &e) override {
            VERIFY(false);
        }
    };

//    Create a proper structure that can work with multiple libraries
    class DBGAlignedReadStorage : public ag::AlignedReadStorage<DBGTraits> {
    private:
        ag::SuffixTracker<DBGTraits> * suffixes = nullptr;
        dbg::CoverageTracker * coverageTracker = nullptr;
        ag::ReadLogger<DBGTraits> *readLogger = nullptr;
        ag::LoggingListener<DBGTraits> *graphLogger = nullptr;
    public:
        ag::AlignedReadStorage<DBGTraits> &getReads() {return *this;}
        const ag::AlignedReadStorage<DBGTraits> &getReads() const {return *this;}
        ag::SuffixTracker<DBGTraits> &getSuffixes() const {return *suffixes;}
        bool tracksSuffixes() const {return suffixes != nullptr;}


        DBGAlignedReadStorage(logging::Logger &logger, size_t threads, SparseDBG &dbg, ag::AlignedReadStorage<DBGTraits> reads, bool _track_cov = false) :
                ag::AlignedReadStorage<DBGTraits>(std::move(reads)) {
            if(_track_cov) {
                coverageTracker = new CoverageTracker(logger, threads, *this, dbg);
            }
        }
        DBGAlignedReadStorage(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                                       std::vector<ag::AlignedRead<DBGTraits>> read_list, bool _track_cov = false) :
                                       ag::AlignedReadStorage<DBGTraits>(logger, threads, dbg, std::move(read_list)) {
            if(_track_cov) {
                coverageTracker = new CoverageTracker(logger, threads, *this, dbg);
            }
        }

        DBGAlignedReadStorage(DBGAlignedReadStorage&&) = default;
        DBGAlignedReadStorage& operator=(DBGAlignedReadStorage&&) = delete;

        void logReads(size_t threads, const std::experimental::filesystem::path& path);
        void logGraph(SparseDBG &dbg, std::ostream &os);
        void trackSuffixes(logging::Logger &logger, size_t threads, SparseDBG &dbg, size_t _min_len, size_t _max_len);
        void stopTrackSuffixes();
        void checkCoverage(const SparseDBG &dbg) const;
        static DBGAlignedReadStorage Load(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path, SparseDBG &dbg,
                                          const IdIndex<Vertex> &index, bool track_cov);
        static DBGAlignedReadStorage Load(logging::Logger &logger, size_t threads, std::istream &is, SparseDBG &dbg,
                                          const IdIndex<Vertex> &index, bool track_cov);
        ~DBGAlignedReadStorage() override;
    };

    template<class I>
    std::vector<ag::AlignedRead<DBGTraits>> AlignReads(logging::Logger &logger, size_t threads, I begin, I end, SparseDBG &dbg, size_t w) {
        logger.info() << "Loading and aligning reads."  << std::endl;
        KmerIndex index(dbg);
        index.fillAnchors(logger, threads, dbg, w);
        OrderedRecordCollector<ag::AlignedRead<DBGTraits>> tmpReads(threads);
        std::function<void(size_t, StringContig &)> read_task = [&index, &tmpReads](size_t pos, StringContig &scontig) {
            Contig contig = scontig.makeContig();
            if (contig.truncSize() < index.minReadLen()) {
                tmpReads.emplace_back(pos, ag::AlignedRead<DBGTraits>(contig.getInnerId(), GraphPath()));
            } else {
                tmpReads.emplace_back(pos, ag::AlignedRead<DBGTraits>(contig.getInnerId(), index.align(contig.getSeq(), contig.getInnerId())));
            }
        };
        processRecords(begin, end, logger, threads, read_task);
        std::vector<ag::AlignedRead<DBGTraits>> res = tmpReads.collectOrdered(logger, threads);
        logger.info() << "Read alignment completed."  << std::endl;
        return std::move(res);
    }
}
