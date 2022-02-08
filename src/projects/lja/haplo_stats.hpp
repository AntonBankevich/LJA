//
// Created by lab44 on 08.02.22.
//
#pragma once

#include <unordered_set>
#include <unordered_map>

#include <utility>
#include <vector>
#include <iostream>
#include <array>
#include <algorithm>
#include "common/logging.hpp"
#include <common/string_utils.hpp>
#include <common/verify.hpp>

inline char other_haplo(char c) {
    if (c == 'm') return 'p';
    else if (c == 'p') return 'm';
    else VERIFY(false);
}

inline char is_defined_haplo(char c){
    return (c == 'm' || c == 'p');
}


struct HaplotypeStats {
    char haplotype;
    std::string label;
//order: p, m, as in file;
    std::array<int, 2> decisive_strips;
    std::array<int, 2> decisive_counts;
    int total_kmers;

//        32         m       0       273     28      390     22      24      112906  17
//#s->seq[i].name, type, s->cnt[i].sc[0], s->cnt[i].sc[1],c[0<<2|2], c[2<<2|0], c[0<<2|1], c[1<<2|0], s->cnt[i].nk, c[0])
//TODO: arr 6, 7, 9
    HaplotypeStats(std::string s) {
        std::vector<std::string> tokens = ::split(s);
        haplotype = tokens[1][0];
        label = tokens[0];
        decisive_strips = std::array<int,2>{stoi(tokens[2]), stoi(tokens[3])};
        decisive_counts = std::array<int,2>{stoi(tokens[4]), stoi(tokens[5])};
        total_kmers = stoi(tokens[8]);
    }
    HaplotypeStats (std::string label, HaplotypeStats a, HaplotypeStats b):label(label) {
        for (size_t i = 0; i < 2; i++) {
            decisive_counts[i] = a.decisive_counts[i] + b.decisive_counts[i];
            decisive_strips[i] = a.decisive_strips[i] + b.decisive_strips[i];
        }
        if (is_defined_haplo(a.haplotype) && is_defined_haplo(b.haplotype)) {
            if (a.haplotype == b.haplotype) {
                haplotype = a.haplotype;
            } else {
                std::cout << "Merging different haplotypes " << a.label << " " << b.label << endl;
                haplotype = 'a';
            }
        } else {
            haplotype = a.haplotype;
            if (!is_defined_haplo(haplotype))
                haplotype = b.haplotype;
        }

    }
    HaplotypeStats():haplotype('u'),label(""), decisive_counts{0,0}, decisive_strips{0,0}, total_kmers(0) {}
    bool is_undefined() {
        return (haplotype != 'm' and haplotype != 'p');
    }

};


inline char AssignBulge(HaplotypeStats top_h, HaplotypeStats bottom_h) {
    size_t pat_count = top_h.decisive_counts[0] * bottom_h.decisive_counts[1];
    size_t mat_count = top_h.decisive_counts[1] * bottom_h.decisive_counts[0];
    size_t top_total = top_h.decisive_counts[0] + top_h.decisive_counts[1];
    size_t bottom_total = bottom_h.decisive_counts[0] + bottom_h.decisive_counts[1];
    char decision = 'a';

    if (mat_count > pat_count || (top_total == 0 && bottom_h.decisive_counts[0] >  bottom_h.decisive_counts[1])
        || (bottom_total == 0 && top_h.decisive_counts[1] > top_h.decisive_counts[0]))
        decision = 'm';
    else if (mat_count < pat_count || (top_total == 0 && bottom_h.decisive_counts[1] >  bottom_h.decisive_counts[0])
             || (bottom_total == 0 && top_h.decisive_counts[0] > top_h.decisive_counts[1]))
        decision = 'p';
    return decision;
}

typedef  std::unordered_map<std::string, HaplotypeStats> haplo_map_type;
//primitive clean graph
