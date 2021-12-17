#pragma once

#include "sequences/sequence.hpp"

inline size_t edit_distance(Sequence s1, Sequence s2) {
    size_t left_skip = 0;
    while(left_skip < s1.size() && left_skip < s2.size() && s1[left_skip] == s2[left_skip]) {
        left_skip++;
    }
    s1 = s1.Subseq(left_skip, s1.size());
    s2 = s2.Subseq(left_skip, s2.size());
    size_t right_skip = 0;
    while(right_skip < s1.size() && right_skip < s2.size() && s1[s1.size() - 1 - right_skip] == s2[s2.size() - 1 - right_skip]) {
        right_skip++;
    }
    s1 = s1.Subseq(0, s1.size() - right_skip);
    s2 = s2.Subseq(0, s2.size() - right_skip);
    std::vector<std::vector<size_t>> d(s1.size() + 1, std::vector<size_t>(s2.size() + 1));
    d[0][0] = 0;
    for(unsigned int i = 1; i <= s1.size(); ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= s2.size(); ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= s1.size(); ++i)
        for(unsigned int j = 1; j <= s2.size(); ++j)
            d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    return d[s1.size()][s2.size()];
}

inline std::pair<size_t, size_t> bestPrefix(const Sequence &s1, const Sequence &_s2) {
    if(_s2.startsWith(s1))
        return {s1.size(), s1.size()};
    Sequence s2 = _s2.Subseq(0, std::min(_s2.size(), s1.size() * 2));
    std::vector<size_t> prev(s2.size() + 1);
    std::vector<size_t> cur(s2.size() + 1);
    for(unsigned int j = 0; j <= s2.size(); ++j) cur[j] = j;
    for(unsigned int i = 1; i <= s1.size(); ++i) {
        std::swap(prev, cur);
        cur[0] = i;
        for(unsigned int j = 1; j <= s2.size(); ++j)
            cur[j] = std::min({ prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    }
    size_t res = s2.size();
    for(size_t j = 0; j <= s2.size(); j++)
        if(cur[j] < cur[res])
            res = j;
    return {res, cur[res]};
}

inline std::pair<size_t, size_t> CheckOverlap(const Sequence &s1, const Sequence &s2, size_t min_overlap, size_t max_overlap, double allowed_divergence) {
    Sequence a = s1.Subseq(s1.size() - std::min(s1.size(), max_overlap));
    Sequence b = s2.Subseq(0, std::min(s2.size(), max_overlap));
    int64_t mult = a.size() + 1;
    int64_t match = 1 * mult;
    int64_t mismatch = 10 * mult;
    int64_t indel = 10 * mult;
    std::vector<int64_t> res(a.size() + 1);
    for(size_t i = 0; i <= a.size(); i++) {
        res[i] = i;
    }
    std::vector<int64_t> prev(a.size() + 1);
    size_t best = 0;
    int64_t best_val = res[a.size()];
    for(size_t j = 1; j <= b.size(); j++) {
        std::swap(prev, res);
        res[0] = prev[0] - indel;
        for(size_t i = 1; i <= a.size(); i++) {
            if(a[i - 1] == b[j - 1]) {
                res[i] = prev[i - 1] + match;
            } else {
                res[i] = std::max(res[i - 1] - indel, std::max(prev[i] - indel, prev[i - 1] - mismatch));
            }
        }
        if(best_val < res[a.size()]) {
            best = j;
            best_val = res[a.size()];
        }
    }
    size_t l1 = a.size() - (best_val % mult);
    size_t l2 = best;
    best_val = best_val / mult * mult;
    double min_val = match * (1 - allowed_divergence) - allowed_divergence * std::max(indel, mismatch);
    if(l1 < min_overlap || l2 < min_overlap || best_val < std::max(l1, l2) * min_val)
        return {0, 0};
    return {l1, l2};
}