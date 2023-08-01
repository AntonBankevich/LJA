#pragma once
#include "ksw2.h"
#include <vector>
#include <cstring>

struct cigar_pair {
    char type;
    size_t length;
    cigar_pair(char type, size_t len):type(type), length(len) {}
};

inline std::vector<cigar_pair> RcCigar(const std::vector<cigar_pair> &cigar) {
    return {cigar.rbegin(), cigar.rend()};
}

inline size_t MaxAlignmentShift(std::vector<cigar_pair> &cigars) {
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

inline double Divergence(const char *tseq, const char *qseq, std::vector<cigar_pair> &cigar) {
    size_t from_pos = 0;
    size_t to_pos = 0;
    size_t match_size = 0;
    for(cigar_pair &cp: cigar) {
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


inline std::vector<cigar_pair> align_ksw(const char *tseq, const char *qseq, int8_t sc_mch, int8_t sc_mis, int gapo, int gape, int width, int end_bonus = 0){
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
    std::vector<cigar_pair> res;
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
//        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
//        ss << ez.cigar[i]>>4 << "MID"[ez.cigar[i]&0xf]);
        res.emplace_back("MID"[ez.cigar[i]&0xf], ez.cigar[i]>>4);
//    putchar('\n');
    free(ez.cigar); free(ts); free(qs);
    return res;
}


class KSWAligner {
private:
    int8_t sc_mch;
    int8_t sc_mis;
    int gapo;
    int gape;
    inline int64_t cost(const char *tseq, const char *qseq, std::vector<cigar_pair> &cigar) {
        size_t from_pos = 0;
        size_t to_pos = 0;
        int64_t res = 0;
        for(cigar_pair &cp: cigar) {
            if(cp.type == 'M') {
                for(size_t i = 0; i < cp.length; i++) {
                    if(tseq[to_pos + i] == qseq[from_pos + i])
                        res += sc_mch;
                    else
                        res -= sc_mis;
                }
            }
            if(cp.type != 'I') {
                to_pos += cp.length;
                res -= gapo + gape * (cp.length - 1);
            }
            if(cp.type != 'D') {
                from_pos += cp.length;
                res -= gapo + gape * (cp.length - 1);
            }
        }
        return res;
    }

public:
    KSWAligner(int8_t scMch, int8_t scMis, int gapo, int gape) : sc_mch(scMch), sc_mis(scMis), gapo(gapo), gape(gape) {}


    inline int64_t cost(const char *tseq, const char *qseq, std::vector<cigar_pair> &cigar) const {
        size_t from_pos = 0;
        size_t to_pos = 0;
        int64_t res = 0;
        for(cigar_pair &cp: cigar) {
            if(cp.type == 'M') {
                for(size_t i = 0; i < cp.length; i++) {
                    if(tseq[to_pos + i] == qseq[from_pos + i])
                        res += sc_mch;
                    else
                        res -= sc_mis;
                }
            }
            if(cp.type != 'I') {
                to_pos += cp.length;
                res -= gapo + gape * (cp.length - 1);
            }
            if(cp.type != 'D') {
                from_pos += cp.length;
                res -= gapo + gape * (cp.length - 1);
            }
        }
        return res;
    }

    std::vector<cigar_pair> align(const char *tseq, const char *qseq, int width, int end_bonus = 0) const {
        return align_ksw(tseq, qseq, sc_mch, sc_mis, gapo, gape, width, end_bonus);
    }

    std::vector<cigar_pair> iterativeBandAlign(const char *tseq, const char *qseq, int min_width, int max_width) const {
        size_t l1 = strlen(tseq);
        size_t l2 = strlen(qseq);
        if(max_width < std::max(l1, l2) - std::min(l1, l2)) {
            if (l1 < l2)
                return {{'M', l1}, {'D', l2 - l1}};
            else
                return {{'M', l2}, {'I', l2 - l1}};
        }
        min_width += std::min<int>(std::max(l1, l2) - std::min(l1, l2) + min_width, max_width);
        int64_t prev_cost = -1000000;
        while(true) {
            auto res = align(tseq, qseq, min_width);
            int64_t new_cost = cost(tseq, qseq, res);
            if(new_cost == prev_cost || min_width == max_width)
                return std::move(res);
            prev_cost = new_cost;
            min_width = std::min(min_width * 2, max_width);
        }
        return align(tseq, qseq, min_width);
    }

    std::vector<cigar_pair> iterativeBandExtend(const char *tseq, const char *qseq, int min_width, int max_width) const {
        int64_t prev_cost = -1000000;
        while(true) {
            auto res = align(tseq, qseq, min_width, 100000);
            int64_t new_cost = cost(tseq, qseq, res);
            if(new_cost == prev_cost)
                return std::move(res);
//          if(MaxAlignmentShift(res) < min_width && Divergence(tseq, qseq, res) < max_divergence) {
//              return std::move(res);
//          }
            prev_cost = new_cost;
            min_width = std::min(min_width * 2, max_width);
        }
        return align(tseq, qseq, min_width, 100000);
    }

    std::vector<cigar_pair> align(const std::string &tseq, const std::string &qseq, int width, int end_bonus) const {
        return align(tseq.c_str(), qseq.c_str(), width, end_bonus);
    }

    std::vector<cigar_pair> iterativeBandAlign(const std::string &tseq, const std::string &qseq, int min_width, int max_width) const {
        return iterativeBandAlign(tseq.c_str(), qseq.c_str(), min_width, max_width);
    }
    std::vector<cigar_pair> iterativeBandExtend(const std::string &tseq, const std::string &qseq, int min_width, int max_width) const {
        return iterativeBandExtend(tseq.c_str(), qseq.c_str(), min_width, max_width);
    }
};

