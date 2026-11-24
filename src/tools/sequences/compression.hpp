#pragma once
#include "common/verify.hpp"
#include <cstddef>
#include <string>

//    Homopolymer compression, then shortening of dimer repeats longer than min_dimer_to_compress to at most
//    max_dimer_size.
struct CompressionParameters {
    bool homopolymer_compressing = false;
    size_t min_dimer_to_compress = 1000000000;
    size_t max_dimer_size = 1000000000;
    size_t dimer_step = 1;

//    The parameters set by a program from its command line, taken by readers unless they are given others.
    static CompressionParameters global;

//    No compression at all: reads as they are.
    static CompressionParameters None() {return {};}

//    Sets the dimer parameters from "<min_dimer_to_compress>,<max_dimer_size>,<dimer_step>".
    void setDimers(const std::string &s);
};

//    Compresses the size symbols of s in place and returns the compressed size. S gives the symbols by get(i) and
//    takes them by set(i, v); it only ever writes a position at or before the one it reads, so packed and character
//    storage compress alike.
template<class S>
size_t CompressInPlace(S &s, size_t size, const CompressionParameters &p) {
    if (!p.homopolymer_compressing)
        return size;
    size_t unique = 0;
    for (size_t i = 0; i < size; i++) {
        if (unique == 0 || s.get(i) != s.get(unique - 1)) {
            s.set(unique, s.get(i));
            unique++;
        }
    }
    size = unique;
    VERIFY(p.min_dimer_to_compress <= p.max_dimer_size);
    VERIFY(p.min_dimer_to_compress >= 4);
    VERIFY(p.dimer_step == 1);
    if (p.min_dimer_to_compress >= size)
        return size;
    size_t cur = 2;
    size_t at_len = 2;
    for (size_t i = 2; i <= size; i++) {
        if (i < size && s.get(i) == s.get(cur - 2)) {
            s.set(cur, s.get(i));
            cur++;
            at_len += 1;
        } else {
            if (at_len > p.min_dimer_to_compress) {
                if (at_len > p.max_dimer_size) {
                    size_t d = 0;
                    if (i == size || at_len == i)
                        d = 2;
                    cur -= (at_len - p.max_dimer_size + 1) / 2 * 2 + d;
                    at_len -= (at_len - p.max_dimer_size + 1) / 2 * 2 + d;
                }
                size_t corr_at_len = at_len / p.dimer_step * p.dimer_step;
                VERIFY(corr_at_len == at_len);
                cur -= (at_len - corr_at_len) / 2 * 2;
                at_len -= (at_len - corr_at_len) / 2 * 2;
            }
            at_len = 2;
            if (i < size) {
                s.set(cur, s.get(i));
                cur++;
            }
        }
    }
    return cur;
}
