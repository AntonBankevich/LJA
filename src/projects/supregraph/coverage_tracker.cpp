#include "coverage_tracker.hpp"
#include <cmath>

using namespace spg;

CoverageTracker::CoverageTracker(AssemblyGraph &graph, ag::AlignedReadStorage &storage, logging::Logger &logger, size_t threads) :
        ag::AlignedReadStorageListener(storage, "CoverageTracker"), SPGListener(graph, "CoverageTracker"),
        graph(graph), storage(&storage) {
    fillLengthHistogram(storage, threads);
    total_bases = windowCapacity(1);
    kpomer_multiplier = multiplier(graph.getK() + 1);
    fillFromStorage(logger, threads);
}

void CoverageTracker::fillLengthHistogram(ag::AlignedReadStorage &storage, size_t threads) {
    omp_set_num_threads(threads);
    std::vector<std::vector<size_t>> thread_histograms(threads);
#pragma omp parallel for default(none) shared(storage, thread_histograms) schedule(dynamic, 100)
    for (size_t i = 0; i < storage.size(); i++) {
        if (!storage[i].valid() || storage[i].getPath().isLegacy())
            continue;
        std::vector<size_t> &hist = thread_histograms[omp_get_thread_num()];
        size_t len = storage[i].getPath().len();
        if (hist.size() <= len)
            hist.resize(len + 1, 0);
        hist[len]++;
    }
    std::vector<size_t> length_histogram;
    for (const std::vector<size_t> &hist : thread_histograms) {
        if (length_histogram.size() < hist.size())
            length_histogram.resize(hist.size(), 0);
        for (size_t l = 0; l < hist.size(); l++)
            length_histogram[l] += hist[l];
    }
    std::vector<size_t> cnt_reads_loeq_than;
    std::vector<size_t> sum_reads_loeq_than;
    cnt_reads_loeq_than.assign(length_histogram.size() + 1, 0);
    sum_reads_loeq_than.assign(length_histogram.size() + 1, 0);
    window_capacity.assign(length_histogram.size() + 1, 0);
    for (size_t l = length_histogram.size(); l > 0; l--) {
        cnt_reads_loeq_than[l-1] = cnt_reads_loeq_than[l] + length_histogram[l - 1];
        sum_reads_loeq_than[l-1] = sum_reads_loeq_than[l] + length_histogram[l - 1] * (l - 1);
        window_capacity[l-1] = sum_reads_loeq_than[l-1] - (l - 2) * cnt_reads_loeq_than[l-1];
    }
}

double CoverageTracker::windowCapacity(size_t s) const {
    if (s >= window_capacity.size())
        return 0;
    return window_capacity[s];
}

void CoverageTracker::adjustSingleSupport(ag::CoverageSamples &info, ag::CoverageSamples::SampleView &view, size_t left, size_t right, __int64_t mult) {
    bool is_kpomer = view.sample.type == ag::CoverageSamples::SampleType::kpomer;
    size_t unit_size = is_kpomer ? graph.getK() + 1 : view.size();
    size_t overlap_left = std::max(view.start(), left);
    size_t overlap_right = std::min(view.finish(), right);
    if (overlap_right >= overlap_left + unit_size) {
        double m = is_kpomer ? kpomer_multiplier : multiplier(unit_size);
        __int64_t raw_delta = __int64_t(overlap_right - overlap_left - unit_size + 1) * mult;
        double delta = double(raw_delta) * m;
//        Reads are replayed concurrently, and several of them may update the same vertex.
#pragma omp atomic
        info.support += delta;
#pragma omp atomic
        view.sample.raw_support += raw_delta;
#pragma omp atomic
        info.raw_support += raw_delta;
    }
}

void CoverageTracker::adjustSupport(Vertex &v, size_t left, size_t right, __int64_t mult) {
    ag::CoverageSamples &info = v.coverage_info;
    const size_t k = graph.getK();
    if (right == v.size()) {
        for (ag::CoverageSamples::SampleView view : info.getRSamples()) {
            if (view.finish() <= left + k)
                break;
            adjustSingleSupport(info, view, left, right, mult);
        }
    } else {
        for (ag::CoverageSamples::SampleView view : info.getSamples()) {
            if (view.start() + k >= right)
                break;
            adjustSingleSupport(info, view, left, right, mult);
        }
    }
}

void CoverageTracker::adjustBothStrands(Vertex &v, size_t left, size_t right, __int64_t mult) {
    adjustSupport(v, left, right, mult);
    adjustSupport(v.rc(), v.size() - right, v.size() - left, mult);
}

//    The rc path visits v.rc() for every v of the path, with mirrored bounds, so both are credited in one walk.
void CoverageTracker::processPath(const GraphPath &path, __int64_t mult) {
//    A legacy path is a subread of a deleted vertex: its windows are gone with the vertex.
    if (!path.valid() || path.isLegacy())
        return;
    if (path.empty()) {
        Vertex &v = path.getStart();
        adjustBothStrands(v, path.leftCut(), v.size() - path.rightCut(), mult);
        return;
    }
    adjustBothStrands(path.getStart(), path.leftCut(), path.getStart().size(), mult);
    for (Vertex &v : path.innerVertices())
        adjustBothStrands(v, 0, v.size(), mult);
    adjustBothStrands(path.getFinish(), 0, path.getFinish().size() - path.rightCut(), mult);
}

void CoverageTracker::fillFromStorage(logging::Logger &logger, size_t threads) {
    omp_set_num_threads(threads);
    logger.info() << "Filling SPG coverage from reads" << std::endl;
#pragma omp parallel for default(none) schedule(dynamic, 100)
    for (size_t i = 0; i < storage->size(); i++) {
        fireAddRead((*storage)[i]);
    }
    logger.info() << "Finished filling SPG coverage from reads" << std::endl;
}

void CoverageTracker::fireAddRead(const ag::AlignedRead &read) {
    processPath(read.getPath(), 1);
}

void CoverageTracker::fireRerouteRead(ag::AlignedRead &read) {
    processPath(read.getPath(), -1);
    processPath(read.getCorrected(), 1);
}

void CoverageTracker::fireInvalidateRead(ag::AlignedRead &read) {
    processPath(read.getPath(), -1);
}

void CoverageTracker::fireEdgeToSupreVertex(Vertex &v, ag::Edge &e) {
//    SampleSegmentTracker fired first; a window on v could not be credited with the reads already on e.
    VERIFY(v.coverage_info.empty());
}

void CoverageTracker::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
//    SampleSegmentTracker gives a new vertex its one zero-support resolution window only when both edges of its pair
//    have several partners; the pair's support is credited to it.
    for (const auto &rec : resolution) {
        Vertex &new_vertex = *rec.first;
        if (!new_vertex.coverage_info.empty()) {
            VERIFY(new_vertex.coverage_info.samples.calculateSize() == 1);
            adjustSupport(new_vertex, 0, new_vertex.size(), __int64_t(rec.second.getSupport()));
        }
    }
}

bool CoverageTracker::checkConsistency(std::ostream &os) {
    for (AlignedRead &read : *storage)
        if (read.valid()) {
            processPath(read.getPath(), -1);
        }
    bool ok = true;
    for (Vertex &v : graph.vertices())
        if (v.hasCoverageInfo()) {
            bool windows_zero = true;
            for (ag::CoverageSamples::SampleView view : v.coverage_info.getSamples())
                windows_zero = windows_zero && view.sample.raw_support == 0;
            if (!windows_zero || v.coverage_info.raw_support != 0 || std::abs(v.coverage_info.support) > 1e-6) {
                os << "CoverageTracker: vertex " << v.getId() << " keeps raw_support=" << v.coverage_info.raw_support
                   << " support=" << v.coverage_info.support << (windows_zero ? "" : " (some window nonzero)")
                   << " after subtracting all stored reads" << std::endl;
                ok = false;
            }
        }
    for (AlignedRead &read : *storage)
        if (read.valid()) {
            processPath(read.getPath(), 1);
        }
    return ok;
}

std::function<std::string(const Vertex &)> CoverageTracker::getVertexLabeler() const {
    std::function<std::string(const Vertex &)> res = [](const Vertex &v) -> std::string {
        if (!v.hasCoverageInfo())
            return "SPGCov:-";
        return "SPGCov:" + std::to_string(v.getSPGCoverage()) + "(raw:" + std::to_string(v.getRawSPGCoverage()) + ")";
    };
    return res;
}
