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

namespace trio {


enum Haplotype: char {
    Paternal = 'p',
    Maternal = 'm',
    Shared = 'a',
    Unknown = '0'
};
inline Haplotype other_haplo(Haplotype c) {
    if (c == Haplotype::Maternal) return Haplotype::Paternal;
    else if (c == Haplotype::Paternal) return Haplotype::Maternal;
    else
        VERIFY(false);
}

inline bool is_defined_haplo(Haplotype c) {
    return (c == Haplotype::Maternal || c == Haplotype::Paternal);
}


struct HaplotypeStats {
    Haplotype haplotype;
    std::vector<multigraph::Edge::id_type> label;
//order: p, m, as in file;
    std::array<int, 2> decisive_strips;
    std::array<int, 2> decisive_counts;
    int total_kmers;

//        32         m       0       273     28      390     22      24      112906  17
//#s->seq[i].name, type, s->cnt[i].sc[0], s->cnt[i].sc[1],c[0<<2|2], c[2<<2|0], c[0<<2|1], c[1<<2|0], s->cnt[i].nk, c[0])
//TODO: arr 6, 7, 9
    HaplotypeStats(std::string s) {
        std::vector<std::string> tokens = ::split(s);
        haplotype = Haplotype(tokens[1][0]);
        label = {Parse<multigraph::MGEdge::id_type>(tokens[0])};
//Numbers of distinctive kmers in strips longer than k - eps
        decisive_strips = std::array<int, 2>{stoi(tokens[2]), stoi(tokens[3])};
//Numbers of distinctive kmers
        decisive_counts = std::array<int, 2>{stoi(tokens[4]), stoi(tokens[5])};
        total_kmers = stoi(tokens[8]);
    }

    void appendKmerStats(HaplotypeStats other) {
        for (size_t i = 0; i < 2; i++) {
            decisive_counts[i] += other.decisive_counts[i];
            decisive_strips[i] += other.decisive_strips[i];
        }
        if (is_defined_haplo(haplotype) && is_defined_haplo(other.haplotype)) {
            if (haplotype != other.haplotype) {
                haplotype = Haplotype::Shared;
            }
        } else {
            if (!is_defined_haplo(haplotype))
                haplotype = other.haplotype;
        }

    }

    HaplotypeStats() : haplotype(Haplotype::Unknown), decisive_counts{0, 0}, decisive_strips{0, 0}, total_kmers(0) {}

    bool is_undefined() {
        return (haplotype != Haplotype::Maternal and haplotype != Haplotype::Paternal);
    }

};


inline Haplotype AssignBulge(HaplotypeStats top_h, HaplotypeStats bottom_h) {
    size_t pat_count = top_h.decisive_counts[0] * bottom_h.decisive_counts[1];
    size_t mat_count = top_h.decisive_counts[1] * bottom_h.decisive_counts[0];
    size_t top_total = top_h.decisive_counts[0] + top_h.decisive_counts[1];
    size_t bottom_total = bottom_h.decisive_counts[0] + bottom_h.decisive_counts[1];
    Haplotype decision = Haplotype::Shared;

    if (mat_count > pat_count || (top_total == 0 && bottom_h.decisive_counts[0] > bottom_h.decisive_counts[1])
        || (bottom_total == 0 && top_h.decisive_counts[1] > top_h.decisive_counts[0]))
        decision = Haplotype::Maternal;
    else if (mat_count < pat_count || (top_total == 0 && bottom_h.decisive_counts[1] > bottom_h.decisive_counts[0])
             || (bottom_total == 0 && top_h.decisive_counts[0] > top_h.decisive_counts[1]))
        decision = Haplotype::Paternal;
    return decision;
}

//primitive clean graph
}