#include "dimer_correction.hpp"
using namespace dbg;
using namespace ag;

Sequence truncSubseq(const GraphPath &path, PathPosition pp, size_t sz) {
    SequenceBuilder sb;
    while (pp != path.lastPosition()) {
        Segment<Edge> seg = path.getSegment(pp);
        if (seg.size() >= sz) {
            sb.append(seg.shrinkRightToLen(sz).truncSeq());
            sz = 0;
            break;
        } else {
            sb.append(seg.truncSeq());
            sz -= seg.size();
        }
        ++pp;
    }
    return sb.BuildSequence();
}

std::string DimerCorrector::correctRead(const std::string &name, ag::GraphPath &path) {
    size_t corrected = 0;
    size_t k = path.getStart().size();
    Sequence remaining_seq = path.truncSeq();
    std::vector<std::string> message;
    for (PathPosition pp = path.firstPosition(); pp != path.lastPosition(); ++pp) {
        Edge &edge = pp.nextEdge();
        if ((pp == path.firstPosition() && !path.startClosed()) || (pp + 1 == path.lastPosition() && !path.endClosed()))
            continue;
        Sequence seq = edge.getStart().getSeq();
        size_t at_cnt1 = 2;
        while (at_cnt1 < k && seq[k - at_cnt1 - 1] == seq[k - at_cnt1 + 1])
            at_cnt1 += 1;
        VERIFY_MSG(at_cnt1 <= max_at,
                   "at_cnt1 < max_at failed " + itos(at_cnt1) + " " + itos(max_at)); // ATAT should be compressed in reads to some length < k
        if (at_cnt1 < 4) //Tandem repeat should be at least 4 nucleotides long
            continue;
        Sequence unit = seq.Subseq(k - 2);
        ag::GraphPath atPrefix(edge.getStart());
        VERIFY(atPrefix.getStart() == edge.getStart());
        atPrefix.extend(unit);
        if (!atPrefix.valid())
            continue;
        Sequence extension = truncSubseq(path, pp, k + max_at);
        size_t at_cnt2 = 0;
        while (at_cnt2 < extension.size() && extension[at_cnt2] == unit[at_cnt2 % 2])
            at_cnt2 += 1;
        VERIFY_MSG(at_cnt2 <= max_at,
                   "at_cnt2 < max_at failed");// ATAT should be compressed in reads to some length max_at < k
        if (at_cnt2 % 2 != 0 || extension.size() < at_cnt2 + k - at_cnt1)
            continue;
        extension = extension.Subseq(0, at_cnt2 + k - at_cnt1);
        ag::GraphPath bulgeSide(edge.getStart());
        bulgeSide.extend(extension);
        VERIFY_MSG(bulgeSide.valid(), "Extension along an existing path failed");
        if (!bulgeSide.endClosed())
            continue;
        Sequence end_seq = extension.Subseq(at_cnt2);
        std::vector<GraphPath> candidates = {bulgeSide};
        if (at_cnt2 > 0) {
            ag::GraphPath candidate(edge.getStart());
            candidate.extend(end_seq);
            if (!candidate.valid())
                continue;
            VERIFY_OMP(candidate.endClosed(), "Candidate alignment end is not closed in case 1");
            candidates.emplace_back(candidate);
        } else {
            size_t max_variation = std::max<size_t>(6, (at_cnt1 + at_cnt2) / 3);
            max_variation = std::min(max_variation, at_cnt1 / 2);
            size_t len = 2;
            while (len <= max_variation && atPrefix.valid()) {
                ag::GraphPath candidate = atPrefix;
                candidate.extend(end_seq);
                if (candidate.valid()) {
                    VERIFY_OMP(candidate.endClosed(), "Candidate alignment end is not closed in case 2");
                    candidates.emplace_back(candidate);
                }
                atPrefix.extend(unit);
                len += 2;
            }
        }
        if (candidates.size() == 1)
            continue;
        size_t best_val = 0;
        size_t best = 0;
        for (size_t i = 0; i < candidates.size(); i++) {
            const ag::SuffixRecord &rec = reads_storage.getSuffixes().getSuffixRecord(candidates[i].frontEdge());
            size_t support = rec.countStartsWith(candidates[i].subPath(candidates[i].firstPosition() + 1, candidates[i].lastPosition()));
            if (support > best_val) {
                best_val = support;
                best = i;
            }
        }
        if (best_val == 0) {
#pragma omp critical
            logger.trace() << "Unsupported path during dinucleotide correction" << std::endl;
        }
        if (best == 0)
            continue;
        message.emplace_back(itos(at_cnt1) + "_" + itos(at_cnt2) + "_" + itos(candidates[best].truncLen()));
        GraphPath tail = path.subPath(pp + candidates[0].calculateSize(), path.lastPosition());
        path.shorten(path.firstPosition(), pp);
        pp = path.lastPosition();
        path += candidates[best];
        path += tail;
        corrected++;
    }
    return join("_", message);
}
