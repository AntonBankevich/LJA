//
// Created by anton on 07.08.2023.
//
#include "ksw_aligner.hpp"

size_t KSWAligner::window = 1000;

__int64_t KSWAligner::cost(const char *tseq, const char *qseq, const AlignmentForm &cigar) const {
    unsigned long from_pos = 0;
    unsigned long to_pos = 0;
    __int64_t res = 0;
    for(const CigarPair &cp: cigar) {
        if(cp.type == CigarEvent::M) {
            for(unsigned long i = 0; i < cp.length; i++) {
                if(tseq[to_pos + i] == qseq[from_pos + i])
                    res += sc_mch;
                else
                    res -= sc_mis;
            }
        }
        if(cp.type != CigarEvent::I) {
            to_pos += cp.length;
            res -= gapo + gape * (cp.length - 1);
        }
        if(cp.type != CigarEvent::D) {
            from_pos += cp.length;
            res -= gapo + gape * (cp.length - 1);
        }
    }
    return res;
}

AlignmentForm KSWAligner::iterativeBandAlign(const char *tseq, const char *qseq) const {
    unsigned long l1 = strlen(tseq);
    unsigned long l2 = strlen(qseq);
    if(max_width < std::max(l1, l2) - std::min(l1, l2)) {
        return {};
    }
    size_t cur_width = min_width;
    cur_width += std::min<int>(std::max(l1, l2) - std::min(l1, l2) + min_width, max_width);
    __int64_t prev_cost = -1000000;
    while(true) {
        auto res = runAlignment(tseq, qseq, cur_width);
        __int64_t new_cost = cost(tseq, qseq, res);
        if(new_cost == prev_cost || cur_width == max_width)
            return std::move(res);
        prev_cost = new_cost;
        cur_width = std::min(cur_width * 2, max_width);
    }
    return runAlignment(tseq, qseq, cur_width);
}

AlignmentForm KSWAligner::iterativeBandExtend(const char *tseq, const char *qseq) const {
    __int64_t prev_cost = -1000000;
    size_t cur_width = min_width;
    while(true) {
        auto res = runAlignment(tseq, qseq, cur_width, 100000);
        __int64_t new_cost = cost(tseq, qseq, res);
        if(new_cost == prev_cost)
            return std::move(res);
//          if(MaxAlignmentShift(res) < min_width && Divergence(tseq, qseq, res) < max_divergence) {
//              return std::move(res);
//          }
        prev_cost = new_cost;
        cur_width = std::min(cur_width * 2, max_width);
    }
    return runAlignment(tseq, qseq, cur_width, 100000);
}

AlignmentForm KSWAligner::extendAlignment(const char *tseq, const char *qseq) const {
    size_t tlen = strlen(tseq);
    size_t qlen = strlen(qseq);
    AlignmentForm res;
    char subquery[window + 1];
    while(true) {
        size_t subqlen = std::min(window, qlen  - res.queryLength());
        memcpy(subquery, qseq + res.queryLength(), subqlen);
        subquery[subqlen] = 0;
        const char *newTarget = tseq + res.targetLength();
        const char *newQuery = subquery;
        AlignmentForm extension = iterativeBandExtend(newTarget, newQuery);
        if(extension.queryLength() == 0 || extension.targetLength() == 0)
            break;
        size_t match = 10;
        AlignmentForm::AlignmentColumnIterator iter = AlignmentHelper::LastLongMatch(newTarget, newQuery, extension, match);
        if(iter == extension.columns().begin())
            break;
        extension = extension.Prefix(iter);
        res += extension;
    }
    return std::move(res);
}

AlignmentForm KSWAligner::runAlignment(const char *tseq, const char *qseq, int width, int end_bonus) const {
    return {align_ksw(tseq, qseq, sc_mch, sc_mis, gapo, gape, width, end_bonus)};
}

AlignmentForm KSWAligner::globalAlignment(const std::string &tseq, const std::string &qseq) const {
    return alignByExtension(tseq.c_str(), qseq.c_str());
}

AlignmentForm KSWAligner::extendAlignment(const std::string &tseq, const std::string &qseq) const {
    return extendAlignment(tseq.c_str(), qseq.c_str());
}

AlignmentForm KSWAligner::alignByExtension(const char *tseq, const char *qseq) const {
    size_t tlen = strlen(tseq);
    size_t qlen = strlen(qseq);
    if(tlen < 1000)
        return iterativeBandAlign(tseq, qseq);
    AlignmentForm initial = extendAlignment(tseq, qseq);
    if(initial.empty() || initial.queryLength() + 50 < qlen || initial.targetLength() + 50 < tlen) {
        return {};
    }
    while(!initial.empty() && initial.targetLength() + 1000 > tlen)
        initial.pop_back();
    initial = initial.Prefix(AlignmentHelper::LastLongMatch(tseq, qseq, initial, 10));
    AlignmentForm extension = iterativeBandAlign(tseq + initial.targetLength(), qseq + initial.queryLength());
    initial += extension;
    return std::move(initial);
}
