#pragma once
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
    };
}