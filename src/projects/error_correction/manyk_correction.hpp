#pragma once
#include "dbg/graph_alignment_storage.hpp"
#include "dbg/sparse_dbg.hpp"
#include "error_correction.hpp"
#include "correction_utils.hpp"
#include "bulge_path_marker.hpp"

class ManyKCorrector : public AbstractCorrectionAlgorithm {
private:
    struct Bulge {
        DBGGraphPath left;
        DBGGraphPath right;
        DBGGraphPath bulge;
        Bulge(DBGGraphPath &&left, DBGGraphPath &&right, DBGGraphPath &&bulge) :
                left(left), right(right), bulge(bulge) {}
    };
    struct Tip {
        DBGGraphPath left;
        DBGGraphPath tip;
        Tip(DBGGraphPath &&left, DBGGraphPath &&tip) : left(left), tip(tip) {}
    };
    class ReadRecord {
    private:
        const DBGGraphPath &read;
    public:
        std::vector<size_t> switch_positions;
        ReadRecord(const DBGGraphPath &read, std::vector<size_t> &&switchPositions) :
                        read(read), switch_positions(switchPositions) {}
        bool isPerfect() const {return blockNum() == 1 && !hasIncomingTip() && !hasOutgoingTip();}
        bool isBad() const {return switch_positions.size() == 0;}
        size_t blockNum() const {return switch_positions.size() / 2;}
        DBGGraphPath getBlock(size_t num) const;
        size_t bulgeNum() const{return (switch_positions.size() - 2) / 2;}
        Bulge getBulge(size_t num);
        bool hasIncomingTip() const {return !switch_positions.empty() && switch_positions[0] > 0;}
        bool hasOutgoingTip() const {return !switch_positions.empty() && switch_positions.back() < read.size();}
        Tip getOutgoingTip();
        Tip getIncomingTip();
    };

    void calculateReliable(const DBGGraphPath &read_path, std::vector<size_t> &last_reliable,
                           std::vector<size_t> &next_reliable) const;
    std::vector<size_t> calculateLowRegions(const std::vector<size_t> &last_reliable,
                                            const std::vector<size_t> &next_reliable,
                                            DBGGraphPath &read_path) const;
    void mergeLow(DBGGraphPath &read_path, std::vector<size_t> &positions, size_t bad_length) const;

    dbg::SparseDBG &dbg;
    RecordStorage &reads;
    size_t K;
    size_t expected_coverage;
    double reliable_threshold;
    double bad_threshold;
    bool diploid;
public:
    ManyKCorrector(logging::Logger &logger, dbg::SparseDBG &dbg, RecordStorage &reads, size_t K, size_t expectedCoverage,
                   double reliable_threshold, double bad_threshold, bool diploid) : AbstractCorrectionAlgorithm("ManyKCorrector"),
                dbg(dbg), reads(reads), K(K), expected_coverage(expectedCoverage),
                reliable_threshold(reliable_threshold), bad_threshold(bad_threshold), diploid(diploid) {
        VERIFY(reads.getMaxLen() >= K);
    }

    void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads) override;

    ReadRecord splitRead(DBGGraphPath &read_path) const;

    DBGGraphPath uniqueExtension(const DBGGraphPath &base, size_t max_len) const;
    DBGGraphPath correctBulgeByBridging(const Bulge &bulge) const;
    DBGGraphPath correctBulgeAsDoubleTip(const Bulge &bulge) const;
    DBGGraphPath correctBulgeWithReliable(const Bulge &bulge) const;
    DBGGraphPath correctBulge(const Bulge &bulge, std::string &message) const;

    DBGGraphPath correctTipWithExtension(const Tip &tip) const;
    DBGGraphPath correctTipWithReliable(const Tip &tip) const;
    DBGGraphPath correctTip(const Tip &tip, std::string &message) const;

    std::string correctRead(DBGGraphPath &read_path) override;
};

size_t ManyKCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,RecordStorage &reads_storage, double threshold,
                double reliable_threshold, size_t K, size_t expectedCoverage, bool diploid);