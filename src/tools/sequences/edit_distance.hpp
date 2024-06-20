#pragma once

#include "sequences/sequence.hpp"

inline size_t edit_distance(Sequence s1, Sequence s2, size_t max_diff) {
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
    std::vector<size_t> prev(s2.size() + 1);
    std::vector<size_t> cur(s2.size() + 1);
    size_t from = 0;
    size_t to = s2.size();
    for(unsigned int j = 0; j <= s2.size(); ++j) cur[j] = j;
    for(unsigned int i = 1; i <= s1.size(); ++i) {
        if(from > to)
            return max_diff;
        std::swap(prev, cur);
        cur[from] = prev[from] + 1;
        for(unsigned int j = from + 1; j <= to; ++j)
            cur[j] = std::min({ prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
        if(to + 1 <= s2.size()) {
            cur[to + 1] = std::min({ cur[to] + 1, prev[to] + (s1[i - 1] == s2[to] ? 0 : 1) });
            to++;
        }
        while(from <= to && cur[from] > max_diff)
            from++;
        while(from <= to && cur[to] > max_diff)
            to--;
    }
    if(from <= s2.size() && to >= s2.size())
        return cur[s2.size()];
    else
        return max_diff;
}

class ArrayWithDefault {
    size_t min_index;
    size_t max_index;
    std::vector<size_t> vals;
    size_t default_value = size_t(-1) / 2;
public:
    size_t &operator[](size_t ind) {
        VERIFY(ind >= min_index);
        VERIFY(ind <= max_index);
        return vals[ind - min_index];
    }
    size_t get(size_t ind) const {
        if(ind <= max_index && ind >= min_index) {
            return vals[ind - min_index];
        } else {
            return default_value;
        }
    }
    size_t minIndex() const {return min_index;}
    size_t maxIndex() const {return max_index;}
    ArrayWithDefault(size_t min_index, size_t max_index) : min_index(min_index), max_index(max_index),
                    vals(max_index - min_index + 1, default_value) {}
};
//inline std::pair<size_t, size_t> bestPrefix(const Sequence &s1, const Sequence &_s2) {
//    if(_s2.startsWith(s1))
//        return {s1.size(), s1.size()};
//    Sequence s2 = _s2.Subseq(0, std::min(_s2.size(), s1.size() * 2));
//    size_t d = std::max<size_t>(std::max(s1.size(), s2.size()) / 50, 20);
//    if(s2.size() < s1.size() - d)
//        return {s2.size(), size_t(-1) / 2};
//    ArrayWithDefault prev(0, d);
//    for(size_t j = 0; j <= d; ++j) prev[j] = j;
//    for(size_t i = 1; i <= s1.size(); ++i) {
//        ArrayWithDefault cur(i - std::min(i, d), std::min(i + d, s2.size()));
//        if(cur.minIndex() == 0)
//            cur[0] = prev.get(0) + 1;
//        for(size_t j = std::max<size_t>(cur.minIndex(), 1); j <= cur.maxIndex(); ++j)
//            cur[j] = std::min({ prev.get(j) + 1, cur.get(j - 1) + 1, prev.get(j - 1) + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
//        std::swap(cur, prev);
//    }
//    size_t res = s2.size();
//    for(size_t j = prev.minIndex(); j <= prev.maxIndex(); j++)
//        if(prev.get(j) < prev.get(res))
//            res = j;
//    return {res, prev.get(res)};
//}

inline std::pair<size_t, size_t> bestPrefix(const Sequence &s1, const Sequence &_s2, size_t max_diff = -1) {
    if(max_diff == -1)
        max_diff = std::max(s1.size(), _s2.size());
    if(_s2.startsWith(s1))
        return {s1.size(), s1.size()};
    Sequence s2 = _s2.Subseq(0, std::min(_s2.size(), s1.size() * 2));
    std::vector<size_t> prev(s2.size() + 1);
    std::vector<size_t> cur(s2.size() + 1);
    size_t from = 0;
    size_t to = s2.size();
    for(unsigned int j = 0; j <= s2.size(); ++j) cur[j] = j;
    for(unsigned int i = 1; i <= s1.size(); ++i) {
        if(from > to)
            return {std::min(s1.size(), s2.size()), std::min(s1.size(), s2.size())};
        std::swap(prev, cur);
        cur[from] = prev[from] + 1;
        for(unsigned int j = from + 1; j <= to; ++j)
            cur[j] = std::min({ prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
        if(to + 1 <= s2.size()) {
            cur[to + 1] = std::min({ cur[to] + 1, prev[to] + (s1[i - 1] == s2[to] ? 0 : 1) });
            to++;
        }
        while(from <= to && cur[from] > max_diff)
            from++;
        while(from <= to && cur[to] > max_diff)
            to--;
    }
    if(from > to)
        return {std::min(s1.size(), s2.size()), std::min(s1.size(), s2.size())};
    size_t res = s2.size();
    for(size_t j = from; j <= to; j++)
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