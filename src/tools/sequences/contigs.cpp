//
// Created by anton on 19.12.2019.
//

#include "contigs.hpp"

bool StringContig::homopolymer_compressing = false;
size_t StringContig::min_dimer_to_compress = 1000000000;
size_t StringContig::max_dimer_size = 1000000000;
size_t StringContig::dimer_step = 1;

void StringContig::compress() {
    if(!homopolymer_compressing)
        return;
    seq.erase(std::unique(seq.begin(), seq.end()), seq.end());
    VERIFY(min_dimer_to_compress <= max_dimer_size);
    VERIFY(min_dimer_to_compress >= 4);
    VERIFY(dimer_step == 1);
    if(min_dimer_to_compress >= seq.size())
        return;
    size_t cur = 2;
    size_t at_len = 2;
    for(size_t i = 2; i <= seq.size(); i++) {
        if (i < seq.size() && seq[i] == seq[cur - 2]) {
            seq[cur] = seq[i];
            cur++;
            at_len += 1;
        } else {
            if(at_len > min_dimer_to_compress) {
                if(at_len > max_dimer_size) {
                    cur -= (at_len - max_dimer_size) / 2 * 2;
                    at_len -= (at_len - max_dimer_size) / 2 * 2;
                }
                size_t corr_at_len = at_len / dimer_step * dimer_step;
                VERIFY(corr_at_len == at_len);
                cur -= (at_len - corr_at_len) / 2 * 2;
                at_len -= (at_len - corr_at_len) / 2 * 2;
            }
            at_len = 2;
            if(i < seq.size()) {
                seq[cur] = seq[i];
                cur++;
            }
        }
    }
    seq.erase(seq.begin() + cur, seq.end());
}
