#include "manyk_correction.hpp"
#include "correction_utils.hpp"
#include "error_correction.hpp"

namespace dbg {
    void ManyKCorrector::calculateReliable(const dbg::GraphPath &read_path, std::vector<PathPosition> &last_reliable,
                                           std::vector<PathPosition> &next_reliable) const {
        for(PathPosition cur = read_path.firstPosition(); cur != read_path.lastPosition(); ++cur) {
            if(cur.nextEdge().getCoverage() >= reliable_threshold || cur.nextEdge().is_reliable) {
                last_reliable.emplace_back(cur + 1);
            } else
                last_reliable.emplace_back(last_reliable.empty() ? read_path.firstPosition() : last_reliable.back());
        }
        for(PathPosition cur = read_path.lastPosition(); cur != read_path.firstPosition(); --cur) {
            if(cur.prevEdge().getCoverage() >= reliable_threshold || cur.prevEdge().is_reliable) {
                next_reliable.emplace_back(cur - 1);
            } else
                next_reliable.emplace_back(next_reliable.empty() ? read_path.lastPosition() : next_reliable.back());
        }
        std::reverse(next_reliable.begin(), next_reliable.end());
    }

    ManyKCorrector::ReadRecord ManyKCorrector::splitRead(const dbg::GraphPath &read_path) const {
        std::vector<PathPosition> last_reliable;
        std::vector<PathPosition> next_reliable;
        calculateReliable(read_path, last_reliable, next_reliable);
        std::vector<PathSegment> bad_regions = calculateLowRegions(last_reliable, next_reliable, read_path);
        mergeLow(read_path, bad_regions, 200);
        if(bad_regions.empty())
            return {std::move(read_path), {{read_path.firstPosition(), read_path.lastPosition()}}};
        std::vector<PathSegment> good_regions;
        if (bad_regions[0].from != read_path.firstPosition())
            good_regions.emplace_back(read_path.firstPosition(), bad_regions.front().from);
        for(size_t i = 0; i + 1 < bad_regions.size(); i++)
            good_regions.emplace_back(bad_regions[i].to, bad_regions[i+1].from);
        if (bad_regions.back().to != read_path.lastPosition())
            good_regions.emplace_back(bad_regions.back().to, read_path.lastPosition());
        return {read_path, std::move(good_regions)};
    }

    void ManyKCorrector::mergeLow(const dbg::GraphPath &read_path, std::vector<PathSegment> &positions, size_t bad_length) const {
        if(positions.empty())
            return;
        size_t new_size = 1;
        for (size_t cur = 0; cur + 1 < positions.size(); cur++) {
            size_t good_len = 0;
            for (PathPosition pp = positions[cur].to; pp != positions[cur+1].from; ++pp) {
                good_len += pp.nextEdge().truncSize();
            }
            if (good_len < bad_length) {
                positions[new_size - 1].to = positions[cur + 1].to;
            } else {
                positions[new_size] = positions[cur + 1];
                new_size++;
            }
        }
        positions.erase(positions.begin() + new_size, positions.end());
    }

    std::vector<ManyKCorrector::PathSegment>
    ManyKCorrector::calculateLowRegions(const std::vector<PathPosition> &last_reliable,
                                        const std::vector<PathPosition> &next_reliable,
                                        const dbg::GraphPath &read_path) const {
        std::vector<PathSegment> positions;
        size_t i = 0;
        for(PathPosition cur = read_path.firstPosition(); cur != read_path.lastPosition(); ++cur, i++) {
            Edge &edge = cur.nextEdge();
            if (edge.getCoverage() < reliable_threshold && !edge.is_reliable &&
                (edge.getStart().inDeg() == 0 || edge.getFinish().outDeg() == 0 ||
                 edge.getCoverage() <= bad_threshold)) {
                PathPosition left = last_reliable[i];
                PathPosition right = next_reliable[i];
                positions.emplace_back(left, right);
                while(cur + 1 != right) {
                    i++;
                    ++cur;
                }
            }
        }
        return std::move(positions);
    }

    std::string ManyKCorrector::correctRead(const std::string &name, dbg::GraphPath &read_path) {
        ReadRecord rr = splitRead(read_path);
        std::string message;
        if (rr.isPerfect() || rr.isBad()) {
            return "";
        }
        std::vector<std::string> messages;
        dbg::GraphPath corrected;
        if (rr.hasIncomingTip()) {
            Tip tip = rr.getIncomingTip();
            std::string tip_message;
            dbg::GraphPath tc = correctTip(tip, tip_message);
            VERIFY(tc.getStart() == tip.tip.getStart());
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
            VERIFY(!corrected.valid() || corrected.getFinish() == bc.getStart());
            corrected += bc;
            corrected += bulge.right;
        }
        if (rr.hasOutgoingTip()) {
            Tip tip = rr.getOutgoingTip();
            std::string tip_message;
            dbg::GraphPath tc = correctTip(tip, tip_message);
            VERIFY(tc.getStart() == tip.tip.getStart());
            VERIFY(tc.front().left == 0);
            if (!tip_message.empty()) {
                messages.emplace_back("o" + tip_message + itos(K));
                messages.emplace_back(itos(tip.tip.truncLen()));
                messages.emplace_back(itos(tc.truncLen()));
            }
            VERIFY(!corrected.valid() || corrected.getFinish() == tc.getStart());
            corrected += tc;
        }
        if (messages.empty())
            return "";
        message = join("_", messages);
        read_path = std::move(corrected);
        return message;
    }

    dbg::GraphPath ManyKCorrector::correctTipWithExtension(const ManyKCorrector::Tip &tip) const {
        const dbg::GraphPath &left = tip.left;
        size_t tlen = tip.tip.truncLen();
        dbg::GraphPath al = uniqueExtension(tip.left, tlen);
        size_t elen = al.truncLen();
        if (elen > 0 && elen + 10 >= tlen) {
            if (elen > tlen) {
                al.cutBack(elen - tlen);
            }
            return std::move(al);
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
        dbg::GraphPath alternative = FindReliableExtension(tip.tip.getStart(), tip.tip.truncLen(), 3);
        if (!alternative.valid())
            return tip.tip;
        if (alternative.truncLen() > tip.tip.truncLen()) {
            alternative.cutBack(alternative.truncLen() - tip.tip.truncLen());
        }
        return std::move(alternative);
    }

    dbg::GraphPath ManyKCorrector::correctTip(const ManyKCorrector::Tip &tip, std::string &message) const {
        dbg::GraphPath correction = correctTipWithExtension(tip);
        VERIFY(tip.tip.getStart() == correction.getStart());
        if (correction != tip.tip) {
            message = "te";
            return std::move(correction);
        }
        correction = correctTipWithReliable(tip);
        VERIFY(tip.tip.getStart() == correction.getStart());
        if (correction != tip.tip) {
            message = "tr";
            return std::move(correction);
        }
        message = "";
        return tip.tip;
    }

    dbg::GraphPath ManyKCorrector::uniqueExtension(const dbg::GraphPath &base, size_t max_len) const {
        dbg::GraphPath al = base;
        al.setCutLeft(0);
        PathPosition cut_pos = al.lastPosition();
        PathPosition start = al.firstPosition();
        size_t extra_len = 0;
        size_t cur_len = base.truncLen() + base.leftCut();
        while (extra_len < max_len) {
            while (cur_len - start.nextEdge().truncSize() >= K) {
                cur_len -= start.nextEdge().truncSize();
                ++start;
            }
            GraphPath cpath = al.subPath(start + 1, al.lastPosition());
            EdgeId next = SuffixSupportedExtension(reads.getSuffixes().getSuffixRecord(start.nextEdge()), cpath, 4, 1);
            if (!next.valid()) {
                break;
            }
            al += *next;
            extra_len += al.back().size();
            cur_len += al.back().size();
        }
        return al.subPath(cut_pos);
    }

    dbg::GraphPath ManyKCorrector::correctBulge(const ManyKCorrector::Bulge &bulge, string &message) const {
        dbg::GraphPath corrected;
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
        corrected = correctBulgeWithReliable(bulge);
        if (corrected != bulge.bulge) {
            message = "br";
            return corrected;
        }
        message = "";
        return bulge.bulge;
    }

    dbg::GraphPath ManyKCorrector::correctBulgeByBridging(const ManyKCorrector::Bulge &bulge) const {
        VERIFY(bulge.bulge.truncLen() < K);
        std::vector<dbg::GraphPath> alternatives1 =
                SuffixSupportedBulgeAlternatives(reads.getSuffixes(), bulge.bulge, 4);
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
        if (rc_left.truncLen() + bulge.bulge.truncLen() > K)
            rc_left.cutBack(rc_left.truncLen() - (K - bulge.bulge.truncLen()));
        dbg::GraphPath right = bulge.right;
        if (right.truncLen() + bulge.bulge.truncLen() > K)
            right.cutBack(right.truncLen() - (K - bulge.bulge.truncLen()));
        for (size_t i = 0; i < alternatives.size(); i++) {
            dbg::GraphPath &al = alternatives[i];
            dbg::GraphPath right_ext = al + right;
            if (reads.getSuffixes().getSuffixRecord(bulge.left.backEdge()).countStartsWith(right_ext) > 0) {
                right_supp++;
                right_best = i;
            }
            dbg::GraphPath left_ext = al.RC() + rc_left;
            if (reads.getSuffixes().getSuffixRecord(right.frontEdge().rc()).countStartsWith(left_ext) > 0) {
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
        dbg::GraphPath left_ext = uniqueExtension(bulge.left, blen + 100);
        dbg::GraphPath right_ext = uniqueExtension(bulge.right.RC(), blen + 100).RC();
        if (left_ext.truncLen() + right_ext.truncLen() < blen + std::min<size_t>(blen, 100) &&
            std::max(left_ext.truncLen(),
                     right_ext.truncLen()) + 100 > blen)
            return bulge.bulge;
        dbg::GraphPath candidate;
        std::vector<EdgeId> left = oneline::map(left_ext.edges().begin(), left_ext.edges().end(), Edge::IdTransformer());
        std::vector<EdgeId> right = oneline::map(right_ext.edges().begin(), right_ext.edges().end(), Edge::IdTransformer());
        for (int shift = -int(right.size()) + 1; shift < int(left.size()); shift++) {
            bool overlap = true;
            for (int i = 0; i < left.size(); i++) {
                if (i - shift >= 0 && i - shift < right.size() &&
                    left[i] != right[i - shift]) {
                    overlap = false;
                    break;
                }
            }
            if (overlap && left.size() > 0 && right.size() > 0) {
                dbg::GraphPath over_al = GraphPath(left.begin(), left.begin() + std::max(0, shift)) +
                        GraphPath(right.begin() + std::max(0, shift) - shift, right.end());;
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
        if (blen > bulge.bulge.getFinish().size() && alternatives.empty()) {
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
                                    dbg::DBGAlignedReadStorage &reads) {
        CoverageReliableFiller cov(reliable_threshold);
        LengthReliableFiller len(20000, 3, 1);
        BridgeReliableFiller bridge(40000);
        ConnectionReliableFiller connect(reliable_threshold);
        BulgePathMarker bulge(dbg, reads, 60000);
        std::vector<AbstractReliableFillingAlgorithm *> algs = {&len, &cov, &bridge, &connect};
        if (diploid)
            algs.emplace_back(&bulge);
        CompositeReliableFiller(std::move(algs)).LoggedReFill(logger, dbg);
    }

    ManyKCorrector::Bulge ManyKCorrector::ReadRecord::getBulge(size_t num) {
        return {getBlock(num), getBlock(num + 1),
                read.subPath(goodRegions[num].to, goodRegions[num + 1].from)};
    }

    ManyKCorrector::Tip ManyKCorrector::ReadRecord::getOutgoingTip() {
        return {getBlock(goodRegions.size() - 1),
                read.subPath(goodRegions.back().to)};
    }

    ManyKCorrector::Tip ManyKCorrector::ReadRecord::getIncomingTip() {
        return {getBlock(0).RC(),
                read.subPath(read.firstPosition(), goodRegions.front().from).RC()};
    }

    dbg::GraphPath ManyKCorrector::ReadRecord::getBlock(size_t num) const {
        return read.subPath(goodRegions[num].from, goodRegions[num].to);
    }

    size_t
    ManyKCorrect(logging::Logger &logger, size_t threads, SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads_storage,
                 double threshold,
                 double reliable_threshold, size_t K, size_t expectedCoverage, bool diploid) {
        logger.info() << "Using K = " << K << " for error correction" << std::endl;
        ManyKCorrector algorithm(logger, dbg, reads_storage, K, expectedCoverage, reliable_threshold, threshold,
                                 diploid);
        return ErrorCorrectionEngine(algorithm).run(logger, threads, dbg, reads_storage);
    }
}