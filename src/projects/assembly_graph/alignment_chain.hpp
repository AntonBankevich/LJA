#pragma once
#include "sequences/contigs.hpp"
namespace ag {
    template<class U, class V>
    class AlignmentChain {
    public:
        Segment <U> seg_from;
        Segment <V> seg_to;

        AlignmentChain(const Segment <U> &seg_from_, const Segment <V> &seg_to_) : seg_from(seg_from_),
                                                                                   seg_to(seg_to_) {
            VERIFY(seg_from_.size() == seg_to_.size());
        }

        AlignmentChain(U& contig_from, V& contig_to, size_t start_from, size_t start_to, size_t len) :
                seg_from(contig_from, start_from, start_from + len), seg_to(contig_to, start_to, start_to + len) {}

        size_t size() { return seg_from.size(); }

        AlignmentChain RC() const {
            return {seg_from.RC(), seg_to.RC()};
        }

        bool operator<(const ag::AlignmentChain<U, V> &other) const {
            if (seg_to != other.seg_to)
                return seg_to < other.seg_to;
            else
                return seg_from < other.seg_from;
        }
        bool operator==(const ag::AlignmentChain<U, V> &other) const {
            return seg_from == other.seg_from && seg_to == other.seg_to;
        }

        bool operator!=(const ag::AlignmentChain<U, V> &other) const {
            return !(*this == other);
        }
    };
}