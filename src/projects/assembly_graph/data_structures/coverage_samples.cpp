#include "coverage_samples.hpp"

#include <utility>

using namespace ag;

size_t CoverageSamples::SampleView::size() const {return sample.finish - sample.start;}

size_t CoverageSamples::SampleView::start() const { return size_t(__int64_t(sample.start) + shift); }
size_t CoverageSamples::SampleView::finish() const { return size_t(__int64_t(sample.finish) + shift); }

void CoverageSamples::SampleIterator::enter() { acc += pos.isBackward() ? -pos->exit_shift : pos->enter_shift; }
void CoverageSamples::SampleIterator::leave() { acc += pos.isBackward() ? -pos->enter_shift : pos->exit_shift; }

CoverageSamples::SampleIterator::SampleIterator(ds::ListPosition<Sample> position) : pos(position) { enter(); }

bool CoverageSamples::SampleIterator::operator==(const SampleIterator &other) const { return pos == other.pos; }
bool CoverageSamples::SampleIterator::operator!=(const SampleIterator &other) const { return !(*this == other); }

CoverageSamples::SampleIterator &CoverageSamples::SampleIterator::operator++() {
    leave();
    ++pos;
    enter();
    return *this;
}

CoverageSamples::SampleIterator CoverageSamples::SampleIterator::operator++(int) {
    SampleIterator t = *this;
    ++(*this);
    return t;
}

CoverageSamples::SampleView CoverageSamples::SampleIterator::operator*() const { return {*pos, acc}; }

IterableStorage<CoverageSamples::SampleIterator> CoverageSamples::getSamples() {
    return {SampleIterator(samples.forward().begin()), SampleIterator(samples.forward().end())};
}

IterableStorage<CoverageSamples::SampleIterator> CoverageSamples::getRSamples() {
    return {SampleIterator(samples.backward().begin()), SampleIterator(samples.backward().end())};
}

void CoverageSamples::addKpomerChunk(size_t start, size_t finish, size_t chunk_support, size_t k, double multiplier) {
    double scaled_support = double(chunk_support) * multiplier;
    samples.push_back(Sample(start, finish, chunk_support, SampleType::kpomer));
    weight += double(finish - start - k);
    support += scaled_support;
    raw_support += chunk_support;
}

void CoverageSamples::addResolutionSample(size_t start, size_t finish, size_t sample_support, double multiplier) {
    double scaled_support = double(sample_support) * multiplier;
    samples.push_back(Sample(start, finish, sample_support, SampleType::resolution));
    weight += 1;
    support += scaled_support;
    raw_support += sample_support;
}

CoverageSamples &CoverageSamples::shift(__int64_t delta) {
    if (!samples.empty()) {
        samples.head()->enter_shift += delta;
        samples.tail()->exit_shift -= delta;
    }
    return *this;
}

CoverageSamples &CoverageSamples::operator+=(CoverageSamples &&other) {
    samples += std::move(other.samples);
    weight += other.weight;
    support += other.support;
    raw_support += other.raw_support;
    return *this;
}
