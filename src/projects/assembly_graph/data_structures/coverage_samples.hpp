#pragma once
#include "common/double_linked_list.hpp"
#include "common/iterator_utils.hpp"
#include <cstddef>

namespace ag {

//    Doubly linked list of records that represent coverage reports for substrings of a vertex.
//    Also stores accumulated weight and average of measurements,
    struct CoverageSamples {
//        Most records are about the substrings of the same size k+1, originating from initial
//        DBG coverage records. Such records are stored in segments of multiple k+1-mers. This
//        type of record is called kpomer. During multiplexing information about individual requests
//        for repeat bridging reads is recorded as resolution type of record..
        enum class SampleType { kpomer, resolution };

//        Sample is a record that stores the coverage of one (resolution) or multiple (kpomer) substrings
//        of a string. To support O(1) unio of lists from different vertices enter and exit shifts are added.
//        True coordinates of samples are derived by iterators.
        struct Sample {
            size_t start = 0;
            size_t finish = 0;
            size_t raw_support = 0;
            __int64_t enter_shift = 0;
            __int64_t exit_shift = 0;
            SampleType type = SampleType::kpomer;
            Sample() = default;
            Sample(size_t start, size_t finish, size_t raw_support, SampleType type)
                : start(start),
                  finish(finish),
                  raw_support(raw_support),
                  type(type) {
            }
        };
        ds::DoublyLinkedList<Sample> samples = {};
        double weight = 0;
//      This integral support value takes into account the bias introduced by difference in segment lengths.
        double support = 0;
//      This integral support value is a direct sum of all support of all substrings. Useful because if there are
//      n reads passing through the entire vertex, the average support will be exactly n.
        size_t raw_support = 0;

        // Temporary objects created by iterators that store true coordinates of substrings.
        struct SampleView {
            Sample &sample;
            __int64_t shift;
            size_t size() const;
            size_t start() const;
            size_t finish() const;
        };

//        Bidirectional iterator over list of records and computes true coordinates of substrings.
        class SampleIterator {
        private:
            ds::ListPosition<Sample> pos;
            __int64_t acc = 0;
            void enter();
            void leave();
        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type        = SampleView;
            using difference_type   = std::ptrdiff_t;
            using pointer           = Sample*;
            using reference         = SampleView;

            explicit SampleIterator(ds::ListPosition<Sample> position);

            bool operator==(const SampleIterator &other) const;
            bool operator!=(const SampleIterator &other) const;
            SampleIterator &operator++();
            SampleIterator operator++(int);

            SampleView operator*() const;
        };

        IterableStorage<SampleIterator> getSamples();
        IterableStorage<SampleIterator> getRSamples();

        void addKpomerChunk(size_t start, size_t finish, size_t chunk_support, size_t k, double multiplier);
        void addResolutionSample(size_t start, size_t finish, size_t sample_support, double multiplier);

//        Shift coordinates of all substrings by delta.
        CoverageSamples &shift(__int64_t delta);
        CoverageSamples &operator+=(CoverageSamples &&other);
    };
}
