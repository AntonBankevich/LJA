#pragma once
#include "ksw2.h"
#include "common/verify.hpp"
#include <vector>
#include <cstring>

enum CigarEvent : char {
    I = 'I', D = 'D', M = 'M',
};

inline CigarEvent CigarEventFromChar(char c) {
    if(c == 'I')
        return I;
    else if (c == 'D')
        return D;
    else
        return M;
}

struct CigarPair {
    CigarEvent type;
    size_t length;
    CigarPair(CigarEvent type, size_t len) : type(type), length(len) {}
    CigarPair(char type, size_t len) : type(CigarEventFromChar(type)), length(len) {}
    CigarPair Reverse() const {
        switch(type) {
            case M:
                return *this;
            case I:
                return {D, length};
            case D:
                return {I, length};
            default:
                VERIFY(false);
        }
    }
};

inline std::vector<CigarPair> RcCigar(const std::vector<CigarPair> &cigar) {
    return {cigar.rbegin(), cigar.rend()};
}

inline size_t MaxAlignmentShift(std::vector<CigarPair> &cigars) {
    int shift = 0;
    int min_shift = 0;
    int max_shift = 0;
    for(size_t i = 0; i + 1 < cigars.size(); i++) {
        if (cigars[i].type == 'D') {
            shift -= cigars[i].length; // NOLINT(cppcoreguidelines-narrowing-conversions)
        } else if (cigars[i].type == 'I') {
            shift += cigars[i].length; // NOLINT(cppcoreguidelines-narrowing-conversions)
        }
        min_shift = std::min(min_shift, shift);
        max_shift = std::max(max_shift, shift);
    }
    min_shift = std::min(min_shift, min_shift - shift);
    max_shift = std::max(max_shift, max_shift - shift);
    return std::max(std::abs(min_shift), std::abs(max_shift));
}

inline double Divergence(const char *tseq, const char *qseq, std::vector<CigarPair> &cigar) {
    size_t from_pos = 0;
    size_t to_pos = 0;
    size_t match_size = 0;
    for(CigarPair &cp: cigar) {
        if(cp.type == 'M') {
            for(size_t i = 0; i < cp.length; i++) {
                if(tseq[to_pos + i] == qseq[from_pos + i])
                    match_size++;
            }
        }
        if(cp.type != 'I') {
            to_pos += cp.length;
        }
        if(cp.type != 'D') {
            from_pos += cp.length;
        }
    }
    size_t len = std::max(from_pos, to_pos);
    return double(len - match_size) / double(len);
}


inline std::vector<CigarPair> align_ksw(const char *tseq, const char *qseq, int8_t sc_mch, int8_t sc_mis, int gapo, int gape, int width, int end_bonus = 0){
    int i; int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int flag = 0;
    if(end_bonus != 0)
        flag = KSW_EZ_EXTZ_ONLY;
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, width, -1, end_bonus, flag, &ez);
    std::vector<CigarPair> res;
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
//        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
//        ss << ez.cigar[i]>>4 << "MID"[ez.cigar[i]&0xf]);
        res.emplace_back("MID"[ez.cigar[i]&0xf], ez.cigar[i]>>4);
//    putchar('\n');
    free(ez.cigar); free(ts); free(qs);
    return res;
}


