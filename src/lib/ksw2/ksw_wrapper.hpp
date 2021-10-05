#include "ksw2.h"
#include <vector>
#include <cstring>

struct cigar_pair {
    char type;
    size_t length;
    cigar_pair(char type, size_t len):type(type), length(len) {}
};
std::vector<cigar_pair> align_ksw(const char *tseq, const char *qseq, int8_t sc_mch, int8_t sc_mis, int gapo, int gape, int width){
    int i; int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
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
    ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, width, -1, 0, 0, &ez);
    std::vector<cigar_pair> res;
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
//        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
//        ss << ez.cigar[i]>>4 << "MID"[ez.cigar[i]&0xf]);
        res.push_back(cigar_pair("MID"[ez.cigar[i]&0xf], ez.cigar[i]>>4));
//    putchar('\n');
    free(ez.cigar); free(ts); free(qs);
    return res;
}
