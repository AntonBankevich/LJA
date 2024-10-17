#include "manyk_correction.hpp"
#include "correction_utils.hpp"
#include "error_correction.hpp"
#include "reliable_filler_interface.hpp"

namespace dbg {
    void ManyKCorrector::calculateReliable(const dbg::GraphPath &read_path, std::vector<size_t> &last_reliable,
                                           std::vector<size_t> &next_reliable) const {
        for (size_t i = 0; i < read_path.size(); i++) {
            const Segment<Edge> &seg = read_path[i];
            Edge &edge = seg.contig();
            if (edge.getCoverage() >= reliable_threshold || edge.is_reliable) {
                last_reliable.emplace_back(i + 1);
                next_reliable.emplace_back(i);
                size_t j = i;
                while (j > 0 && next_reliable[j - 1] == read_path.size()) {
                    next_reliable[j - 1] = i;
                    j--;
                }
            } else {
                if (i == 0)
                    last_reliable.emplace_back(0);
                else
                    last_reliable.emplace_back(last_reliable.back());
                next_reliable.emplace_back(read_path.size());
            }
        }
    }

    ManyKCorrector::ReadRecord ManyKCorrector::splitRead(dbg::GraphPath &read_path) const {
        std::vector<size_t> last_reliable;
        std::vector<size_t> next_reliable;
        calculateReliable(read_path, last_reliable, next_reliable);
        std::vector<size_t> positions = calculateLowRegions(last_reliable, next_reliable, read_path);
        mergeLow(read_path, positions, 200);
        if (positions.size() > 0 && positions[0] == 0) {
            positions = std::vector<size_t>(positions.begin() + 1, positions.end());
        } else {
            positions.insert(positions.begin(), 0);
        }
        if (positions.size() > 0 && positions.back() == read_path.size()) {
            positions.pop_back();
        } else {
            positions.push_back(read_path.size());
        }
        return {std::move(read_path), std::move(positions)};
    }

    void ManyKCorrector::mergeLow(dbg::GraphPath &read_path, std::vector<size_t> &positions, size_t bad_length) const {
        size_t new_size = 2;
        for (size_t cur = 2; cur < positions.size(); cur += 2) {
            size_t good_len = 0;
            for (size_t i = positions[cur - 1]; i < positions[cur]; i++) {
                const Segment<Edge> &seg = read_path[i];
                good_len += seg.size();
            }
            if (good_len < bad_length) {
                positions[new_size - 1] = positions[cur + 1];
            } else {
                positions[new_size] = positions[cur];
                positions[new_size + 1] = positions[cur + 1];
                new_size += 2;
            }
        }
        if (!positions.empty()) {
            positions.resize(new_size);
        }
    }

    std::vector<size_t>
    ManyKCorrector::calculateLowRegions(const std::vector<size_t> &last_reliable,
                                        const std::vector<size_t> &next_reliable,
                                        dbg::GraphPath &read_path) const {
        std::vector<size_t> positions;
        for (size_t i = 0; i < read_path.size(); i++) {
            const Segment<Edge> &seg = read_path[i];
            Edge &edge = seg.contig();
            if (edge.getCoverage() < reliable_threshold && !edge.is_reliable &&
                (edge.getStart().inDeg() == 0 || edge.getFinish().outDeg() == 0 ||
                 edge.getCoverage() <= bad_threshold)) {
                size_t left = last_reliable[i];
                size_t right = next_reliable[i];
                if (positions.empty() || left > positions.back()) {
                    positions.emplace_back(left);
                    positions.emplace_back(right);
                } else {
                    VERIFY(left == positions[positions.size() - 2])
                    VERIFY(right == positions.back());
                }
            }
        }
        return std::move(positions);
    }

    std::string ManyKCorrector::correctRead(dbg::GraphPath &read_path) {
        ReadRecord rr = splitRead(read_path);
        std::string message = "";
        if (rr.isPerfect() || rr.isBad()) {
            return "";
        }
        std::vector<std::string> messages;
        dbg::GraphPath corrected;
        if (rr.hasIncomingTip()) {
            Tip tip = rr.getIncomingTip();
            std::string tip_message;
            dbg::GraphPath tc = correctTip(tip, tip_message);
            VERIFY(tc.start() == tip.tip.start());
            VERIFY(tc.front().left == 0);
            if (!tip_message.empty()) {
                messages.emplace_back("i" + tip_message + itos(K));
                messages.emplace_back(itos(tip.tip.truncLen()));
                messages.emplace_back(itos(tc.truncLen()));
            }
            corrected += tc.RC();
        }
        corrected += rr.getBlock(0);
        for (size_t i = 0; i < rr.bulgeNum(); i++) {
            Bulge bulge = rr.getBulge(i);
            std::string bulge_message;
            dbg::GraphPath bc = correctBulge(bulge, bulge_message);
            if (!bulge_message.empty()) {
                messages.emplace_back(bulge_message + itos(K));
                messages.emplace_back(itos(bulge.bulge.truncLen()));
                messages.emplace_back(itos(bc.truncLen()));
            }
            VERIFY(!corrected.valid() || corrected.finish() == bc.start());
            corrected += bc;
            corrected += bulge.right;
        }
        if (rr.hasOutgoingTip()) {
            Tip tip = rr.getOutgoingTip();
            std::string tip_message;
            dbg::GraphPath tc = correctTip(tip, tip_message);
            VERIFY(tc.start() == tip.tip.start());
            VERIFY(tc.front().left == 0);
            if (!tip_message.empty()) {
                messages.emplace_back("o" + tip_message + itos(K));
                messages.emplace_back(itos(tip.tip.truncLen()));
                messages.emplace_back(itos(tc.truncLen()));
            }
            VERIFY(!corrected.valid() || corrected.finish() == tc.start());
            corrected += tc;
        }
        if (messages.empty())
            return "";
        message = join("_", messages);
        if (corrected.truncLen() < 100) {
#pragma omp critical
            std::cout << corrected.truncLen() << " " << message << "\n " << read_path.covStr(true) << "\n "
                      << corrected.covStr(
                              true) << std::endl;
        }
        read_path = corrected;
        return message;
    }

    dbg::GraphPath ManyKCorrector::correctTipWithExtension(const ManyKCorrector::Tip &tip) const {
        const dbg::GraphPath &left = tip.left;
        size_t tlen = tip.tip.truncLen();
        dbg::GraphPath al = uniqueExtension(tip.left, tlen);
        size_t elen = al.truncLen() - tip.left.truncLen();
        if (elen > 0 && elen + 10 >= tlen) {
            if (elen > tlen) {
                al.cutBack(elen - tlen);
            }
            return al.subPath(tip.left.size(), al.size());
        } else {
            return tip.tip;
        }
    }

    dbg::GraphPath ManyKCorrector::correctTipWithReliable(const ManyKCorrector::Tip &tip) const {
        size_t tlen = tip.tip.truncLen();
//    std::vector<dbg::GraphAlignment> alternatives = FindPlausibleTipAlternatives(tip.tip, std::max<size_t>(tlen / 100, 20), 3);
//    if(alternatives.size() == 1) {
//        if(alternatives[0].len() > tip.tip.len())
//            alternatives[0].cutBack(alternatives[0].len() - tip.tip.len());
//        return alternatives[0];
//    } else
//        return tip.tip;
        dbg::GraphPath alternative = FindReliableExtension(tip.tip.start(), tip.tip.truncLen(), 3);
        if (!alternative.valid())
            return tip.tip;
        if (alternative.truncLen() > tip.tip.truncLen()) {
            alternative.cutBack(alternative.truncLen() - tip.tip.truncLen());
        }
        return std::move(alternative);
    }

    dbg::GraphPath ManyKCorrector::correctTip(const ManyKCorrector::Tip &tip, std::string &message) const {
        dbg::GraphPath correction = correctTipWithReliable(tip);
        VERIFY(tip.tip.start() == correction.start());
        if (correction != tip.tip) {
            message = "tr";
            return std::move(correction);
        }
        correction = correctTipWithExtension(tip);
        VERIFY(tip.tip.start() == correction.start());
        if (correction != tip.tip) {
            message = "te";
            return std::move(correction);
        }
        message = "";
        return tip.tip;
    }

    dbg::GraphPath ManyKCorrector::uniqueExtension(const dbg::GraphPath &base, size_t max_len) const {
        dbg::GraphPath al = base;
        size_t start = 0;
        size_t extra_len = 0;
        size_t cur_len = base.truncLen();
        while (extra_len < max_len) {
            while (cur_len - al[start].size() >= K) {
                cur_len -= al[start].size();
                start++;
            }
            VERIFY(start < al.size());
            CompactPath cpath = CompactPath::Subpath(al, start, al.size());
            unsigned char next = reads.getRecord(cpath.start()).getUniqueExtension(cpath.cpath(), 4, 1);
            if (next == (unsigned char) -1) {
                break;
            }
            al += al.finish().getOutgoing(next);
            extra_len += al.back().size();
            cur_len += al.back().size();
        }
        return al;
    }

    dbg::GraphPath ManyKCorrector::correctBulge(const ManyKCorrector::Bulge &bulge, string &message) const {
        dbg::GraphPath corrected;
        corrected = correctBulgeWithReliable(bulge);
        if (corrected != bulge.bulge) {
            message = "br";
            return corrected;
        }
        if (bulge.bulge.truncLen() + 100 < K) {
            corrected = correctBulgeByBridging(bulge);
            if (corrected != bulge.bulge) {
                message = "bb";
                return corrected;
            }
        }
        corrected = correctBulgeAsDoubleTip(bulge);
        if (corrected != bulge.bulge) {
            message = "bd";
            return corrected;
        }
        message = "";
        return bulge.bulge;
    }

    dbg::GraphPath ManyKCorrector::correctBulgeByBridging(const ManyKCorrector::Bulge &bulge) const {
        VERIFY(bulge.bulge.truncLen() < K);
        std::vector<dbg::GraphPath> alternatives1 = reads.getRecord(bulge.bulge.start()).
                getBulgeAlternatives(bulge.bulge.finish(), 4);
        std::vector<dbg::GraphPath> alternatives;
        for (dbg::GraphPath &al: alternatives1) {
            if (al.truncLen() + 100 < bulge.bulge.truncLen() && bulge.bulge.truncLen() < al.truncLen() + 100)
                alternatives.emplace_back(std::move(al));
        }
        if (alternatives.empty())
            return bulge.bulge;
        if (alternatives.size() == 1)
            return std::move(alternatives[0]);
        size_t left_supp = 0;
        size_t right_supp = 0;
        size_t left_best = 0;
        size_t right_best = 0;
        dbg::GraphPath rc_left = bulge.left.RC();
        if (rc_left.truncLen() > K - bulge.bulge.truncLen())
            rc_left.cutBack(rc_left.truncLen() - (K - bulge.bulge.truncLen()));
        CompactPath crc_left(rc_left);
        dbg::GraphPath right = bulge.right;
        if (right.size() > K - bulge.bulge.truncLen())
            right.cutBack(right.truncLen() - (K - bulge.bulge.truncLen()));
        CompactPath cright(right);
        for (size_t i = 0; i < alternatives.size(); i++) {
            dbg::GraphPath &al = alternatives[i];
            Sequence right_ext = CompactPath(al).cpath() + cright.cpath();
            if (reads.getRecord(bulge.bulge.start()).countStartsWith(right_ext) > 0) {
                right_supp++;
                right_best = i;
            }
            Sequence left_ext = CompactPath(al.RC()).cpath() + crc_left.cpath();
            if (reads.getRecord(bulge.bulge.finish().rc()).countStartsWith(left_ext) > 0) {
                left_supp++;
                left_best = i;
            }
        }
        if ((left_supp == 1 && right_supp == 1 && left_best == right_best) || left_supp + right_supp == 1) {
            size_t best = std::max(left_best, right_best);
            return alternatives[best];
        } else {
            return bulge.bulge;
        }
    }

    dbg::GraphPath ManyKCorrector::correctBulgeAsDoubleTip(const ManyKCorrector::Bulge &bulge) const {
        size_t blen = bulge.bulge.truncLen();
        dbg::GraphPath left_ext = uniqueExtension(bulge.left, blen + 100).subPath(bulge.left.size());
        dbg::GraphPath right_ext = uniqueExtension(bulge.right.RC(), blen + 100).subPath(bulge.right.size()).RC();
        if (left_ext.truncLen() + right_ext.truncLen() < blen + std::min<size_t>(blen, 100) &&
            std::max(left_ext.truncLen(),
                     right_ext.truncLen()) + 100 > blen)
            return bulge.bulge;
        dbg::GraphPath candidate;
        for (int shift = -int(right_ext.size()) + 1; shift < int(left_ext.size()); shift++) {
            bool overlap = true;
            for (int i = 0; i < left_ext.size(); i++) {
                if (i - shift >= 0 && i - shift < right_ext.size() &&
                    left_ext[i].contig() != right_ext[i - shift].contig()) {
                    overlap = false;
                    break;
                }
            }
            if (overlap && left_ext.size() > 0 && right_ext.size() > 0) {
                dbg::GraphPath over_al = left_ext.subPath(0, std::max(0, shift)) +
                                         right_ext.subPath(std::max(0, shift) - shift);
                size_t over_len = over_al.truncLen();
                if (over_len < blen + 100 && blen < over_len + 100) {
                    if (candidate.valid())
                        return bulge.bulge;
                    candidate = std::move(over_al);
                }
            }
        }
        if (candidate.valid())
            return candidate;
        else
            return bulge.bulge;
    }

    dbg::GraphPath ManyKCorrector::correctBulgeWithReliable(const ManyKCorrector::Bulge &bulge) const {
        size_t blen = bulge.bulge.truncLen();
        std::vector<dbg::GraphPath> alternatives = FindPlausibleBulgeAlternatives(bulge.bulge,
                                                                                  std::max<size_t>(blen / 100, 20), 3);
        if (blen > bulge.bulge.finish().size() && alternatives.empty()) {
            alternatives = FindPlausibleBulgeAlternatives(bulge.bulge, blen / 10 + 32, 3);
            if (alternatives.empty()) {
                alternatives = FindPlausibleBulgeAlternatives(bulge.bulge, blen / 5 + 32, 3);
            }
        }
        if (alternatives.size() == 1)
            return alternatives[0];
        else
            return bulge.bulge;
    }

    void ManyKCorrector::initialize(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                                    dbg::ReadAlignmentStorage &reads) {
        reliableFiller->LoggedReFill(logger, dbg);
//        CoverageReliableFiller cov(reliable_threshold);
//        LengthReliableFiller len(20000, 3, 1);
//        BridgeReliableFiller bridge(40000);
//        ConnectionReliableFiller connect(reliable_threshold);
//        BulgePathMarker bulge(dbg, reads, 60000);
//        std::vector<AbstractReliableFillingAlgorithm *> algs = {&len, &cov, &bridge, &connect};
//        if (diploid)
//            algs.emplace_back(&bulge);
//        CompositeReliableFiller(std::move(algs)).LoggedReFill(logger, dbg);
    }

    ManyKCorrector::Bulge ManyKCorrector::ReadRecord::getBulge(size_t num) {
        return {read.subPath(switch_positions[num * 2], switch_positions[num * 2 + 1]),
                read.subPath(switch_positions[num * 2 + 2], switch_positions[num * 2 + 3]),
                read.subPath(switch_positions[num * 2 + 1], switch_positions[num * 2 + 2])};
    }

    ManyKCorrector::Tip ManyKCorrector::ReadRecord::getOutgoingTip() {
        return {read.subPath(switch_positions[switch_positions.size() - 2], switch_positions.back()),
                read.subPath(switch_positions.back(), read.size())};
    }

    ManyKCorrector::Tip ManyKCorrector::ReadRecord::getIncomingTip() {
        return {read.subPath(switch_positions[0], switch_positions[1]).RC(),
                read.subPath(0, switch_positions[0]).RC()};
    }

    dbg::GraphPath ManyKCorrector::ReadRecord::getBlock(size_t num) const {
        return read.subPath(switch_positions[num * 2], switch_positions[num * 2 + 1]);
    }

    size_t
    ManyKCorrect(logging::Logger &logger, size_t threads, SparseDBG &dbg, dbg::ReadAlignmentStorage &reads_storage,
                 AbstractReliableFillingAlgorithm &reliableFiller,
                 double threshold,
                 double reliable_threshold, size_t K, size_t expectedCoverage, bool diploid) {
        logger.info() << "Using K = " << K << " for error correction" << std::endl;
        ManyKCorrector algorithm(logger, dbg, reads_storage, reliableFiller, K, expectedCoverage, reliable_threshold, threshold,
                                 diploid);
        return ErrorCorrectionEngine(algorithm).run(logger, threads, dbg, reads_storage);
    }
}
