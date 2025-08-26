#pragma once
#include "dbg/sparse_dbg.hpp"
namespace dbg {
    inline void InvalidateBad(logging::Logger &logger, size_t threads, ag::AlignedReadStorage &reads, size_t min_read_size,
                              const std::function<bool(const Edge &)> &is_bad, const std::string &message) {
        omp_set_num_threads(threads);
        ParallelCounter cnt(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(cnt, reads, message, is_bad, min_read_size, std::cout)
        for (size_t i = 0; i < reads.size(); i++) {
            ag::AlignedRead &alignedRead = reads[i];
            if (!alignedRead.valid())
                continue;
            GraphPath &al = alignedRead.getPath();
            ag::PathPosition left = al.firstPosition();
            ag::PathPosition right = al.lastPosition();
            while (left != right && is_bad(left.nextEdge())) {
                ++left;
            }
            while (left != right && is_bad(right.prevEdge())) {
                --right;
            }
            bool middle_bad = false;
            size_t len = 0;
            for (ag::PathPosition pp = left; pp != right; ++pp) {
                if (is_bad(pp.nextEdge())) {
                    middle_bad = true;
                    break;
                }
                len += al.getSegment(pp).size();
            }
            if (middle_bad || left == right) {
                reads.delayedInvalidateRead(alignedRead, message);
                cnt += 1;
            } else if (left != al.firstPosition() || right != al.lastPosition()) {
                ag::GraphPath sub = al.subPath(left, right);
                if (sub.truncLen() < min_read_size) {
                    reads.delayedInvalidateRead(alignedRead, message);
                } else {
                    std::string extra_message =
                            message + "_EndsClipped_" + itos(al.subPath(al.firstPosition(), left).truncLen()) + "_" +
                                        itos(al.subPath(right).truncLen());
                    reads.rerouteRead(alignedRead, sub, extra_message);
                }
                cnt += 1;
            }
        }
        reads.applyCorrections(logger, threads);
        if (cnt.get() > 0)
            logger.info() << "Could not correct " << cnt.get() << " reads. They were removed or truncated."
                          << std::endl;
    }

    inline void InvalidateLowCovered(logging::Logger &logger, size_t threads, ag::AlignedReadStorage &reads, double threshold,
                              size_t min_read_size, const std::string &message) {
        const std::function<bool(const Edge &)> &is_bad = [threshold](const Edge &edge) {
            return edge.getCoverage() < threshold;
        };
        InvalidateBad(logger, threads, reads, min_read_size, is_bad, message);
    }

    inline void InvalidateSubreads(logging::Logger &logger, size_t threads, ag::AlignedReadStorage &storage, SparseDBG &graph) {
        ag::SuffixTracker tracker(storage, graph, 0, 1000000000);
        tracker.fillFromStorage(logger, threads);
        for (ag::AlignedRead &alignedRead: storage) {
            const ag::SuffixRecord &rec = tracker.getSuffixRecord(alignedRead.getPath().frontEdge());
            size_t cnt = rec.countStartsWith(alignedRead.getPath().subPath(alignedRead.getPath().firstPosition() + 1));
            VERIFY_MSG(cnt >= 1, "This function assumes that suffixes are stored for complete reads")
            if (cnt >= 2) {
                storage.delayedInvalidateRead(alignedRead, "Subread");
            }
        }
        storage.applyCorrections(logger, threads);
        logger.info() << "Uncorrected reads were removed." << std::endl;
    }

}