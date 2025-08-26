#pragma once
#include "dbg/dbg_read_alignment_storage.hpp"
#include "dbg/sparse_dbg.hpp"
#include "error_correction.hpp"
#include "correction_utils.hpp"
#include "bulge_path_marker.hpp"

namespace dbg {
    class ManyKCorrector : public AbstractCorrectionAlgorithm {
    private:
        struct Bulge {
            ag::GraphPath left;
            ag::GraphPath right;
            ag::GraphPath bulge;

            Bulge(ag::GraphPath &&left, ag::GraphPath &&right, ag::GraphPath &&bulge) :
                    left(left), right(right), bulge(bulge) {}
        };

        struct Tip {
            ag::GraphPath left;
            ag::GraphPath tip;

            Tip(ag::GraphPath &&left, ag::GraphPath &&tip) : left(left), tip(tip) {}
        };

        struct PathSegment {
            ag::PathPosition from;
            ag::PathPosition to;
            PathSegment(ag::PathPosition from, ag::PathPosition to) : from(from), to(to) {}
        };

        class ReadRecord {
        private:
            const ag::GraphPath &read;
        public:
            std::vector<PathSegment> goodRegions;

            ReadRecord(const ag::GraphPath &read, std::vector<PathSegment> goodRegions) :
                    read(read), goodRegions(std::move(goodRegions)) {}

            bool isPerfect() const { return blockNum() == 1 && !hasIncomingTip() && !hasOutgoingTip(); }

            bool isBad() const { return goodRegions.size() == 0; }

            size_t blockNum() const { return goodRegions.size(); }

            ag::GraphPath getBlock(size_t num) const;

            size_t bulgeNum() const { return goodRegions.size() - 1; }

            Bulge getBulge(size_t num);

            bool hasIncomingTip() const { return !goodRegions.empty() && goodRegions.front().from != read.firstPosition(); }

            bool hasOutgoingTip() const { return !goodRegions.empty() && goodRegions.back().to != read.lastPosition(); }

            Tip getOutgoingTip();

            Tip getIncomingTip();
        };

        void calculateReliable(const ag::GraphPath &read_path, std::vector<ag::PathPosition> &last_reliable,
                               std::vector<ag::PathPosition> &next_reliable) const;

        std::vector<PathSegment> calculateLowRegions(const std::vector<ag::PathPosition> &last_reliable,
                                                     const std::vector<ag::PathPosition> &next_reliable,
                                                     const ag::GraphPath &read_path) const;

        void mergeLow(const ag::GraphPath &read_path, std::vector<PathSegment> &positions, size_t bad_length) const;

        dbg::SparseDBG &dbg;
        DBGAlignedReadStorage &reads;
        size_t K;
        size_t expected_coverage;
        double reliable_threshold;
        double bad_threshold;
        bool diploid;
    public:
        ManyKCorrector(logging::Logger &logger, dbg::SparseDBG &dbg, DBGAlignedReadStorage &reads, size_t K,
                       size_t expectedCoverage,
                       double reliable_threshold, double bad_threshold, bool diploid) : AbstractCorrectionAlgorithm(
                "ManyKCorrector"),
                                                                                        dbg(dbg), reads(reads), K(K),
                                                                                        expected_coverage(
                                                                                                expectedCoverage),
                                                                                        reliable_threshold(
                                                                                                reliable_threshold),
                                                                                        bad_threshold(bad_threshold),
                                                                                        diploid(diploid) {
//            VERIFY(reads.getMaxLen() >= K);
        }

        void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, DBGAlignedReadStorage &reads) override;

        ReadRecord splitRead(const ag::GraphPath &read_path) const;

        ag::GraphPath uniqueExtension(const ag::GraphPath &base, size_t max_len) const;

        ag::GraphPath correctBulgeByBridging(const Bulge &bulge) const;

        ag::GraphPath correctBulgeAsDoubleTip(const Bulge &bulge) const;

        ag::GraphPath correctBulgeWithReliable(const Bulge &bulge) const;

        ag::GraphPath correctBulge(const Bulge &bulge, std::string &message) const;

        ag::GraphPath correctTipWithExtension(const Tip &tip) const;

        ag::GraphPath correctTipWithReliable(const Tip &tip) const;

        ag::GraphPath correctTip(const Tip &tip, std::string &message) const;

        std::string correctRead(const std::string &name, ag::GraphPath &read_path) override;
    };

    size_t ManyKCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, DBGAlignedReadStorage &reads_storage,
                        double threshold,
                        double reliable_threshold, size_t K, size_t expectedCoverage, bool diploid);
}