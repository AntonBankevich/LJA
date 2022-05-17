#include "dimer_correction.hpp"
using namespace dbg;

std::string DimerCorrector::correctRead(GraphAlignment &path) {
    size_t corrected = 0;
    size_t k = path.start().seq.size();
    Sequence remaining_seq = path.truncSeq();
    std::vector<std::string> message;
    for (size_t path_pos = 0; path_pos < path.size(); path_pos++) {
        if (path[path_pos].left > 0 || path[path_pos].right < path[path_pos].contig().size())
            continue;
        Sequence seq = path.getVertex(path_pos).seq;
        size_t at_cnt1 = 2;
        while (at_cnt1 < k && seq[k - at_cnt1 - 1] == seq[k - at_cnt1 + 1])
            at_cnt1 += 1;
        VERIFY_MSG(at_cnt1 <= max_at,
                   "at_cnt1 < max_at failed " + itos(at_cnt1) + " " + itos(max_at)); // ATAT should be compressed in reads to some length < k
        if (at_cnt1 < 4) //Tandem repeat should be at least 4 nucleotides long
            continue;
        Sequence unit = seq.Subseq(k - 2);
        GraphAlignment atPrefix(path.getVertex(path_pos));
        atPrefix.extend(unit);
        if (!atPrefix.valid())
            continue;
        Sequence extension = path.truncSeq(path_pos, k + max_at);
        size_t at_cnt2 = 0;
        while (at_cnt2 < extension.size() && extension[at_cnt2] == unit[at_cnt2 % 2])
            at_cnt2 += 1;
        VERIFY_MSG(at_cnt2 <= max_at,
                   "at_cnt2 < max_at failed");// ATAT should be compressed in reads to some length max_at < k
        if (at_cnt2 % 2 != 0 || extension.size() < at_cnt2 + k - at_cnt1)
            continue;
        extension = extension.Subseq(0, at_cnt2 + k - at_cnt1);
        GraphAlignment bulgeSide(path.getVertex(path_pos));
        bulgeSide.extend(extension);
        VERIFY_MSG(bulgeSide.valid(), "Extension along an existing path failed");
        if (!bulgeSide.endClosed())
            continue;
        Sequence end_seq = extension.Subseq(at_cnt2);
        std::vector<CompactPath> candidates = {CompactPath(bulgeSide)};
        if (at_cnt2 > 0) {
            GraphAlignment candidate(path.getVertex(path_pos));
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
                GraphAlignment candidate = atPrefix;
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
        const VertexRecord &rec = reads_storage.getRecord(path.getVertex(path_pos));
        for (size_t i = 0; i < candidates.size(); i++) {
            size_t support = rec.countStartsWith(candidates[i].cpath());
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
        message.emplace_back(itos(at_cnt1) + "_" + itos(at_cnt2) + "_" + itos(candidates[best].getPath().len()));
        path = path.reroute(path_pos, path_pos + candidates[0].size(), candidates[best].getPath());
        corrected++;
    }
    return join("_", message);
}
