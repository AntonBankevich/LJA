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

inline std::pair<size_t, size_t> bestPrefix(const Sequence &s1, const Sequence &s2) {
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
