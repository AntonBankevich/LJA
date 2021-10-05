#pragma once
#include "dbg/graph_alignment_storage.hpp"
#include "dbg/sparse_dbg.hpp"

class ManyKCorrector {
private:
    struct Bulge {
        GraphAlignment left;
        GraphAlignment right;
        GraphAlignment bulge;
        Bulge(GraphAlignment &&left, GraphAlignment &&right, GraphAlignment &&bulge) :
                left(left), right(right), bulge(bulge) {}
    };
    struct Tip {
        GraphAlignment left;
        GraphAlignment tip;
        Tip(GraphAlignment &&left, GraphAlignment &&tip) : left(left), tip(tip) {}
    };
    class ReadRecord {
    public:
        GraphAlignment read;
        std::vector<size_t> switch_positions;
        ReadRecord(GraphAlignment &&read, std::vector<size_t> &&switchPositions) :
                        read(read), switch_positions(switchPositions) {}
        bool isPerfect() const {return blockNum() == 1 && !hasIncomingTip() && !hasOutgoingTip();}
        bool isBad() const {return switch_positions.size() == 0;}
        size_t blockNum() const {return switch_positions.size() / 2;}
        GraphAlignment getBlock(size_t num) const;
        size_t bulgeNum() const{return (switch_positions.size() - 2) / 2;}
        Bulge getBulge(size_t num);
        bool hasIncomingTip() const {return !switch_positions.empty() && switch_positions[0] > 0;}
        bool hasOutgoingTip() const {
            return !switch_positions.empty() && switch_positions.back() < read.size();
        }
        Tip getOutgoingTip();
        Tip getIncomingTip();
    };

    void calculateReliable(const GraphAlignment &read_path, std::vector<size_t> &last_reliable,
                           std::vector<size_t> &next_reliable) const;
    std::vector<size_t> calculateLowRegions(const std::vector<size_t> &last_reliable,
                                            const std::vector<size_t> &next_reliable,
                                            GraphAlignment &read_path) const;
    void mergeLow(GraphAlignment &read_path, std::vector<size_t> &positions, size_t bad_length) const;

    dbg::SparseDBG &dbg;
    RecordStorage &reads;
    size_t K;
    size_t expected_coverage;
    double reliable_threshold;
    double bad_threshold;
public:
    ManyKCorrector(SparseDBG &dbg, RecordStorage &reads, size_t K, size_t expectedCoverage,
                   double reliable_threshold, double bad_threshold) :
                dbg(dbg), reads(reads), K(K), expected_coverage(expectedCoverage),
                reliable_threshold(reliable_threshold), bad_threshold(bad_threshold) {
        VERIFY(reads.getMaxLen() >= K);
    }

    ReadRecord splitRead(GraphAlignment &&read_path) const;

    GraphAlignment uniqueExtension(const GraphAlignment &base, size_t max_len) const;
    GraphAlignment correctBulgeByBridging(const Bulge &bulge) const;
    GraphAlignment correctBulgeAsDoubleTip(const Bulge &bulge) const;
    GraphAlignment correctBulgeWithReliable(const Bulge &bulge) const;
    GraphAlignment correctBulge(const Bulge &bulge, std::string &message) const;

    GraphAlignment correctTipWithExtension(const Tip &tip) const;
    GraphAlignment correctTipWithReliable(const Tip &tip) const;
    GraphAlignment correctTip(const Tip &tip, std::string &message) const;

    GraphAlignment correctRead(GraphAlignment &&read_path, std::string &message) const;
};

size_t ManyKCorrect(logging::Logger &logger, SparseDBG &dbg,RecordStorage &reads_storage, double threshold,
                double reliable_threshold, size_t K, size_t expectedCoverage, size_t threads);