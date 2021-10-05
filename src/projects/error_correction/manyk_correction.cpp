#include "manyk_correction.hpp"
#include "correction_utils.hpp"

void ManyKCorrector::calculateReliable(const GraphAlignment &read_path, std::vector<size_t> &last_reliable,
                                       std::vector<size_t> &next_reliable) const {
    for(size_t i = 0; i < read_path.size(); i++) {
        const Segment<Edge> &seg= read_path[i];
        Edge & edge = seg.contig();
        if (edge.getCoverage() >= reliable_threshold || edge.is_reliable) {
            last_reliable.emplace_back(i + 1);
            next_reliable.emplace_back(i);
            size_t j = i;
            while(j > 0 && next_reliable[j - 1] == read_path.size()) {
                next_reliable[j - 1] = i;
                j--;
            }
        } else {
            if(i == 0)
                last_reliable.emplace_back(0);
            else
                last_reliable.emplace_back(last_reliable.back());
            next_reliable.emplace_back(read_path.size());
        }
    }
}

ManyKCorrector::ReadRecord ManyKCorrector::splitRead(GraphAlignment &&read_path) const {
    std::vector<size_t> last_reliable;
    std::vector<size_t> next_reliable;
    calculateReliable(read_path, last_reliable, next_reliable);
    std::vector<size_t> positions = calculateLowRegions(last_reliable, next_reliable, read_path);
    mergeLow(read_path, positions, 200);
    if(positions.size() > 0 && positions[0] == 0) {
        positions = std::vector<size_t>(positions.begin() + 1, positions.end());
    } else {
        positions.insert(positions.begin(), 0);
    }
    if(positions.size() > 0 && positions.back() == read_path.size()) {
        positions.pop_back();
    } else {
        positions.push_back(read_path.size());
    }
    return {std::move(read_path), std::move(positions)};
}

void ManyKCorrector::mergeLow(GraphAlignment &read_path, std::vector<size_t> &positions, size_t bad_length) const {
    size_t new_size = 2;
    for(size_t cur = 2; cur < positions.size(); cur+= 2) {
        size_t good_len = 0;
        for(size_t i = positions[cur - 1]; i < positions[cur]; i++) {
            const Segment<Edge> &seg = read_path[i];
            good_len += seg.size();
        }
        if(good_len < bad_length) {
            positions[new_size - 1] = positions[cur + 1];
        } else {
            positions[new_size] = positions[cur];
            positions[new_size + 1] = positions[cur + 1];
            new_size += 2;
        }
    }
    if(!positions.empty()) {
        positions.resize(new_size);
    }
}

std::vector<size_t>
ManyKCorrector::calculateLowRegions(const std::vector<size_t> &last_reliable, const std::vector<size_t> &next_reliable,
                                    GraphAlignment &read_path) const {
    std::vector<size_t> positions;
    for(size_t i = 0; i < read_path.size(); i++) {
        const Segment<Edge> &seg = read_path[i];
        Edge &edge = seg.contig();
        if (edge.getCoverage() < reliable_threshold && !edge.is_reliable &&
            (edge.start()->inDeg() == 0 || edge.end()->outDeg() == 0 || edge.getCoverage() <= bad_threshold)) {
            size_t left = last_reliable[i];
            size_t right = next_reliable[i];
            if(positions.empty() || left > positions.back()) {
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

GraphAlignment ManyKCorrector::correctRead(GraphAlignment &&read_path, string &message) const {
    ReadRecord rr = splitRead(std::move(read_path));
    message = "";
    if(rr.isPerfect() || rr.isBad()) {
        return std::move(rr.read);
    }
    std::vector<std::string> messages;
    GraphAlignment corrected;
    if(rr.hasIncomingTip()) {
        Tip tip = rr.getIncomingTip();
        std::string tip_message;
        GraphAlignment tc = correctTip(tip, tip_message);
        VERIFY(tc.start() == tip.tip.start());
        VERIFY(tc.front().left == 0);
        if(!tip_message.empty()) {
            messages.emplace_back("i" + tip_message + itos(K));
            messages.emplace_back(itos(tip.tip.len()));
            messages.emplace_back(itos(tc.len()));
        }
        corrected += tc.RC();
    }
    corrected += rr.getBlock(0);
    for(size_t i = 0; i < rr.bulgeNum(); i++) {
        Bulge bulge = rr.getBulge(i);
        std::string bulge_message;
        GraphAlignment bc = correctBulge(bulge, bulge_message);
        if(!bulge_message.empty()) {
            messages.emplace_back(bulge_message + itos(K));
            messages.emplace_back(itos(bulge.bulge.len()));
            messages.emplace_back(itos(bc.len()));
        }
        VERIFY(!corrected.valid() || corrected.finish() == bc.start());
        corrected += bc;
        corrected += bulge.right;
    }
    if(rr.hasOutgoingTip()) {
        Tip tip = rr.getOutgoingTip();
        std::string tip_message;
        GraphAlignment tc = correctTip(tip, tip_message);
        VERIFY(tc.start() == tip.tip.start());
        VERIFY(tc.front().left == 0);
        if(!tip_message.empty()) {
            messages.emplace_back("o" + tip_message + itos(K));
            messages.emplace_back(itos(tip.tip.len()));
            messages.emplace_back(itos(tc.len()));
        }
        VERIFY(!corrected.valid() || corrected.finish() == tc.start());
        corrected += tc;
    }
    message = join("_", messages);
    if(corrected.len() < 100) {
#pragma omp critical
            std::cout << corrected.len() << " " << message << "\noppa " << rr.read.str(true) << "\noppa " << corrected.str(true) << std::endl;
    }
    return std::move(corrected);
}

GraphAlignment ManyKCorrector::correctTipWithExtension(const ManyKCorrector::Tip &tip) const {
    const GraphAlignment &left = tip.left;
    size_t tlen = tip.tip.len();
    GraphAlignment al = uniqueExtension(tip.left, tlen);
    size_t elen = al.len() - tip.left.len();
    if(elen > 0 && elen + 10 >= tlen) {
        if(elen > tlen) {
            al.cutBack(elen - tlen);
        }
        return al.subalignment(tip.left.size(), al.size());
    } else {
        return tip.tip;
    }
}

GraphAlignment ManyKCorrector::correctTipWithReliable(const ManyKCorrector::Tip &tip) const {
    size_t tlen = tip.tip.len();
//    std::vector<dbg::GraphAlignment> alternatives = FindPlausibleTipAlternatives(tip.tip, std::max<size_t>(tlen / 100, 20), 3);
//    if(alternatives.size() == 1) {
//        if(alternatives[0].len() > tip.tip.len())
//            alternatives[0].cutBack(alternatives[0].len() - tip.tip.len());
//        return alternatives[0];
//    } else
//        return tip.tip;
    GraphAlignment alternative = FindReliableExtension(tip.tip.start(), tip.tip.len(), 4);
    if(!alternative.valid())
        return tip.tip;
    if(alternative.len() > tip.tip.len()) {
        alternative.cutBack(alternative.len() - tip.tip.len());
    }
    return std::move(alternative);
}

GraphAlignment ManyKCorrector::correctTip(const ManyKCorrector::Tip &tip, std::string &message) const {
    GraphAlignment correction = correctTipWithExtension(tip);
    VERIFY(tip.tip.start() == correction.start());
    if(correction != tip.tip) {
        message = "te";
        return std::move(correction);
    }
    correction = correctTipWithReliable(tip);
    VERIFY(tip.tip.start() == correction.start());
    if(correction != tip.tip) {
        message = "tr";
        return std::move(correction);
    }
    message = "";
    return tip.tip;
}

GraphAlignment ManyKCorrector::uniqueExtension(const GraphAlignment &base, size_t max_len) const {
    GraphAlignment al = base;
    size_t start = 0;
    size_t extra_len = 0;
    size_t cur_len = base.len();
    while(extra_len < max_len){
        while(cur_len - al[start].size() >= K) {
            cur_len -= al[start].size();
            start++;
        }
        VERIFY(start < al.size());
        CompactPath cpath(al, start, al.size());
        unsigned char next = reads.getRecord(cpath.start()).getUniqueExtension(cpath.cpath(), 4, 1);
        if(next == (unsigned char)-1) {
            break;
        }
        al += al.finish().getOutgoing(next);
        extra_len += al.back().size();
        cur_len += al.back().size();
    }
    return al;
}

GraphAlignment ManyKCorrector::correctBulge(const ManyKCorrector::Bulge &bulge, string &message) const {
    GraphAlignment corrected;
    if(bulge.bulge.len() + 100 < K) {
        corrected = correctBulgeByBridging(bulge);
        if(corrected != bulge.bulge) {
            message = "bb";
            return corrected;
        }
    }
    corrected = correctBulgeAsDoubleTip(bulge);
    if(corrected != bulge.bulge) {
        message = "bd";
        return corrected;
    }
    corrected = correctBulgeWithReliable(bulge);
    if(corrected != bulge.bulge) {
        message = "br";
        return corrected;
    }
    message = "";
    return bulge.bulge;
}

GraphAlignment ManyKCorrector::correctBulgeByBridging(const ManyKCorrector::Bulge &bulge) const {
    VERIFY(bulge.bulge.len() < K);
    std::vector<GraphAlignment> alternatives = reads.getRecord(bulge.bulge.start()).
            getBulgeAlternatives(bulge.bulge.finish(), 4);
    if(alternatives.empty())
        return bulge.bulge;
    if(alternatives.size() == 1)
        return std::move(alternatives[0]);
    size_t left_supp = 0;
    size_t right_supp = 0;
    size_t left_best = 0;
    size_t right_best = 0;
    GraphAlignment rc_left = bulge.left.RC();
    if(rc_left.len() > K - bulge.bulge.len())
        rc_left.cutBack(rc_left.len() - (K - bulge.bulge.len()));
    CompactPath crc_left(rc_left);
    GraphAlignment right = bulge.right;
    if(right.size() > K - bulge.bulge.len())
        right.cutBack(right.len() - (K - bulge.bulge.len()));
    CompactPath cright(right);
    for(size_t i = 0; i < alternatives.size(); i++) {
        GraphAlignment &al = alternatives[i];
        Sequence right_ext = CompactPath(al).cpath() + cright.cpath();
        if(reads.getRecord(bulge.bulge.start()).countStartsWith(right_ext) > 0) {
            right_supp++;
            right_best = i;
        }
        Sequence left_ext = CompactPath(al.RC()).cpath() + crc_left.cpath();
        if(reads.getRecord(bulge.bulge.finish().rc()).countStartsWith(left_ext) > 0) {
            left_supp++;
            left_best = i;
        }
    }
    if((left_supp == 1 && right_supp == 1 && left_best == right_best) || left_supp + right_supp == 1) {
        size_t best = std::max(left_best, right_best);
        return alternatives[best];
    } else {
        return bulge.bulge;
    }
}

GraphAlignment ManyKCorrector::correctBulgeAsDoubleTip(const ManyKCorrector::Bulge &bulge) const {
    size_t blen = bulge.bulge.len();
    GraphAlignment left_ext = uniqueExtension(bulge.left, blen + 100).subalignment(bulge.left.size());
    GraphAlignment right_ext = uniqueExtension(bulge.right.RC(), blen + 100).subalignment(bulge.right.size()).RC();
    if(left_ext.len() + right_ext.len() < blen + std::min<size_t>(blen, 500))
        return bulge.bulge;
    GraphAlignment candidate;
    for(int shift = -int(right_ext.size()) + 1; shift < int(left_ext.size()); shift++) {
        bool overlap = true;
        for(int i = 0; i < left_ext.size(); i++) {
            if(i - shift >= 0 && i - shift < right_ext.size() &&
               left_ext[i].contig() != right_ext[i - shift].contig()) {
                overlap = false;
                break;
            }
        }
        if(overlap && left_ext.size() > 0 && right_ext.size() > 0) {
            GraphAlignment over_al =left_ext.subalignment(0, std::max(0, shift)) +
                                    right_ext.subalignment(std::max(0, shift) - shift);
            size_t over_len = over_al.len();
            if(over_len < blen + 100 && blen < over_len + 100){
                if(candidate.valid())
                    return bulge.bulge;
                candidate = std::move(over_al);
            }
        }
    }
    if(candidate.valid())
        return candidate;
    else
        return bulge.bulge;
}

GraphAlignment ManyKCorrector::correctBulgeWithReliable(const ManyKCorrector::Bulge &bulge) const {
    size_t blen = bulge.bulge.len();
    std::vector<dbg::GraphAlignment> alternatives = FindPlausibleBulgeAlternatives(bulge.bulge, std::max<size_t>(blen / 100, 20), 3);
    if(alternatives.size() == 1)
        return alternatives[0];
    else
        return bulge.bulge;
}

ManyKCorrector::Bulge ManyKCorrector::ReadRecord::getBulge(size_t num) {
    return  {read.subalignment(switch_positions[num * 2], switch_positions[num * 2 + 1]),
             read.subalignment(switch_positions[num * 2 + 2], switch_positions[num * 2 + 3]),
             read.subalignment(switch_positions[num * 2 + 1], switch_positions[num * 2 + 2])};
}

ManyKCorrector::Tip ManyKCorrector::ReadRecord::getOutgoingTip() {
    return {read.subalignment(switch_positions[switch_positions.size() - 2], switch_positions.back()),
            read.subalignment(switch_positions.back(), read.size())};
}

ManyKCorrector::Tip ManyKCorrector::ReadRecord::getIncomingTip() {
    return {read.subalignment(switch_positions[0], switch_positions[1]).RC(),
            read.subalignment(0, switch_positions[0]).RC()};
}

GraphAlignment ManyKCorrector::ReadRecord::getBlock(size_t num) const {
    return read.subalignment(switch_positions[num * 2], switch_positions[num * 2 + 1]);
}

size_t ManyKCorrect(logging::Logger &logger, SparseDBG &dbg, RecordStorage &reads_storage, double threshold,
                    double reliable_threshold, size_t K, size_t expectedCoverage, size_t threads) {
    FillReliableWithConnections(logger, dbg, reliable_threshold);
    logger.info() << "Correcting low covered regions in reads with K = " << K << std::endl;
    ManyKCorrector corrector(dbg, reads_storage, K, expectedCoverage, reliable_threshold, threshold);
    ParallelRecordCollector<std::string> results(threads);
    ParallelCounter cnt(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(std::cout, corrector, reads_storage, results, threshold, logger, reliable_threshold, cnt)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        std::stringstream ss;
        std::vector<std::string> messages;
        AlignedRead &alignedRead = reads_storage[read_ind];
        if (!alignedRead.valid())
            continue;
        CompactPath &initial_cpath = alignedRead.path;
        std::string message;
        GraphAlignment corrected = corrector.correctRead(initial_cpath.getAlignment(), message);
        if(!message.empty()) {
            reads_storage.reroute(alignedRead, corrected, message);
            cnt += 1;
        }
    }
    reads_storage.applyCorrections(logger, threads);
    logger.info() << "Corrected low covered regions in " << cnt.get() << " reads with K = " << K << std::endl;
    return cnt.get();
}
