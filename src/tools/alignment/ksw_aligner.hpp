#pragma once

#include <ksw2/ksw_wrapper.hpp>
#include "alignment_form.hpp"
#include <vector>
#include <cctype>
#include <string>

class KSWAligner {
private:
    __int8_t sc_mch;
    __int8_t sc_mis;
    int gapo;
    int gape;
    static unsigned long window;
    size_t max_width = 10000;
    size_t min_width = 10;

    AlignmentForm runAlignment(const char *tseq, const char *qseq, int width, int end_bonus = 0) const;

    AlignmentForm iterativeBandAlign(const char *tseq, const char *qseq) const;
    AlignmentForm iterativeBandExtend(const char *tseq, const char *qseq) const;
    AlignmentForm extendAlignment(const char *tseq, const char *qseq) const;
    AlignmentForm alignByExtension(const char *tseq, const char *qseq) const;
public:
    KSWAligner(__int8_t scMch, __int8_t scMis, int gapo, int gape) :
            sc_mch(scMch), sc_mis(scMis), gapo(gapo), gape(gape) {}

    inline __int64_t cost(const char *tseq, const char *qseq, const AlignmentForm &cigar) const;

    AlignmentForm globalAlignment(const std::string &tseq, const std::string &qseq) const;
    AlignmentForm extendAlignment(const std::string &tseq, const std::string &qseq) const;
};