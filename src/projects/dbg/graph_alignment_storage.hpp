#pragma once
#include "sparse_dbg.hpp"
#include "dbg_graph_aligner.hpp"
#include "assembly_graph/graph_alignment_storage.hpp"

namespace dbg {
class ReadAlignmentStorage : public ag::RecordStorage<DBGTraits> {
public:
    ReadAlignmentStorage(SparseDBG &dbg, size_t _min_len, size_t _max_len,
            ag::ReadLogger &readLogger, bool _track_cov = false, bool log_changes = false,
    bool track_suffixes = true) : ag::RecordStorage<DBGTraits>(dbg, _min_len, _max_len, readLogger, _track_cov, log_changes, track_suffixes) {
    }

    ReadAlignmentStorage(ReadAlignmentStorage&&) = default;
    ReadAlignmentStorage& operator=(ReadAlignmentStorage&&) = default;

    template<class I>
    void FillAlignments(logging::Logger &logger, size_t threads, I begin, I end, SparseDBG &dbg, KmerIndex &index) {
        VERIFY(index.alignmentReady());
        if (track_cov) {
            logger.info() << "Cleaning edge coverages" << std::endl;
            for (Edge &edge: dbg.edges()) {
                edge.incCov(-edge.intCov());
            }
        }
        logger.info() << "Collecting alignments of sequences to the graph" << std::endl;
        if (track_suffixes) {
            logger.info() << "Storing suffixes of read paths of length up to " << this->max_len << std::endl;
        }
        ParallelRecordCollector<std::tuple<size_t, std::string, CompactPath>> tmpReads(threads);
        ParallelCounter cnt(threads);
        std::function<void(size_t, StringContig &)> read_task = [this, &tmpReads, &cnt, &index](size_t pos,
                                                                                                StringContig &scontig) {
            Contig contig = scontig.makeContig();
            if (contig.truncSize() < index.minReadLen()) {
                tmpReads.emplace_back(pos, contig.getInnerId(), CompactPath());
                return;
            }
            GraphPath path = index.align(contig.getSeq());
            CompactPath cpath(path);
            GraphPath rcPath = path.RC();
            CompactPath crcPath(rcPath);
            addSubpath(cpath);
            addSubpath(crcPath);
            cnt += cpath.size();
            tmpReads.emplace_back(pos, contig.getInnerId(), cpath);
        };
        processRecords(begin, end, logger, threads, read_task);
        reads.resize(tmpReads.size());
        for (auto &rec: tmpReads) {
            VERIFY(std::get<0>(rec) < reads.size());
            reads[std::get<0>(rec)] = {std::get<1>(rec), std::move(std::get<2>(rec))};
        }
        logger.info() << "Alignment collection finished. Total length of alignments is " << cnt.get() << std::endl;
    }

    void checkCoverage(const SparseDBG &dbg) const;

};
}