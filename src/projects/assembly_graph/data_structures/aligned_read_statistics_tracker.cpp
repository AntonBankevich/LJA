#include "aligned_read_statistics_tracker.hpp"

void ag::AlignedReadStatisticsTracker::addPath(const GraphPath &path, __int64_t mult) {
    if (!path.valid())
        return;
    if (path.empty()) {
        Vertex &v = path.getStart();
#pragma omp atomic
        v.subread_length += path.len() * mult;
#pragma omp atomic
        v.subread_count += mult;
        // if (v.outDeg() == 1 && v.front().isSuffix()) {
        //     if (path.rightCut() <= v.frontVertex().size()) {
        //         VERIFY(path.leftCut() < v.size() - v.frontVertex().size());
        //         v.front().rc().outgoing_read_count += mult;
        //     }
        // }
    } else {
        for (Vertex &v : path.innerVertices()) {
#pragma omp atomic
            v.covering_read_count+=mult;
        }
        VERIFY(!path.frontEdge().isPrefix());
        VERIFY(!path.backEdge().isSuffix());
#pragma omp atomic
        path.backEdge().read_tail_length += (path.backEdge().truncSize() - path.rightCut()) * mult;
#pragma omp atomic
        path.backEdge().read_tail_count += mult;
        for (Edge &e : path.edges()) {
            if (!e.isSuffix()) {
#pragma omp atomic
                e.outgoing_read_count += mult;
            }
        }
        if (!path.frontEdge().isSuffix()) {
#pragma omp atomic
            path.frontEdge().outgoing_read_count -= mult;
        }
        // if (!path.frontEdge().isSuffix() && path.leftCut() > path.getStart().size()) {
        //     std::cout << path.frontEdge() << std::endl;
        //
        //     path.frontEdge().outgoing_read_count -= mult;
        //     VERIFY(path.frontEdge().outgoing_read_count <1000000000);
        // }
    }
}

void ag::AlignedReadStatisticsTracker::fillFromStorage(logging::Logger &logger, size_t threads) {
    omp_set_num_threads(threads);
    logger.info() << "Filling edge coverages" << std::endl;
#pragma omp parallel for default(none) schedule(dynamic, 100)
    for(size_t i = 0; i < storage->size(); i++) {
        fireAddRead((*storage)[i]);
    }
    logger.info() << "Finished filling edge coverages" << std::endl;
}

void ag::AlignedReadStatisticsTracker::processNewOuterVertex(Vertex &new_vertex) {
    for (AlignedReadDirection dir : storage->getSubstringReads(new_vertex.getId())) {
        new_vertex.subread_length += dir.getPath().len();
        new_vertex.subread_count += 1;
    }
    Edge &inc = new_vertex.incFront();
    for (AlignedReadDirection dir : storage->getOutgoingReads(inc.rc())) {
        inc.read_tail_length += inc.truncSize() - dir.leftCut();
        inc.read_tail_count += 1;
    }
}

ag::AlignedReadStatisticsTracker::AlignedReadStatisticsTracker(logging::Logger &logger, size_t threads,
                                                               ag::AssemblyGraph &graph, ag::AlignedReadStorage &storage, ag::SuffixTracker &suffix_tracker):
    ag::AlignedReadStorageListener(storage, "AlignedReadStatisticsTracker"),
    ag::ResolutionListener(graph, "AlignedReadStatisticsTracker"), storage(&storage),
    suffix_tracker(&suffix_tracker) {
    fillFromStorage(logger, threads);
    for (Edge &edge : graph.edges()) {
        fireAddEdge(edge);
    }
}

void ag::AlignedReadStatisticsTracker::fireAddRead(const ag::AlignedRead &read) {
    addPath(read.getPath().RC());
    addPath(read.getPath());
}

void ag::AlignedReadStatisticsTracker::fireRerouteRead(ag::AlignedRead &read) {
    addPath(read.getPath().RC(), -1);
    addPath(read.getPath(), -1);
    addPath(read.getCorrected().RC());
    addPath(read.getCorrected());
}

void ag::AlignedReadStatisticsTracker::fireInvalidateRead(ag::AlignedRead &read) {
    addPath(read.getPath().RC(), -1);
    addPath(read.getPath(), -1);
}

void ag::AlignedReadStatisticsTracker::fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) {
    if (new_edge.isPrefix()) {
        for (Edge &edge : path.edges()) {
            new_edge.read_tail_length += edge.read_tail_length;
            new_edge.read_tail_count += edge.read_tail_count;
        }
    }
    if (!new_edge.isSuffix()) {
        new_edge.min_equivalent_size = path.frontEdge().min_equivalent_size;
        new_edge.outgoing_read_count = path.frontEdge().outgoing_read_count;
    }
}

void ag::AlignedReadStatisticsTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    processNewOuterVertex(v);
    const SuffixRecord &rec = suffixTracker().getSuffixRecord(v.rc().front());
    v.covering_read_count = rec.getNumberOfPaths() - v.incFront().read_tail_count;
    Edge &inc = v.incFront();
    inc.min_equivalent_size = e.min_equivalent_size;
    inc.outgoing_read_count = e.outgoing_read_count;
}

void ag::AlignedReadStatisticsTracker::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    for (Vertex &v : path.innerVertices()) {
        new_vertex.subread_length += v.subread_length;
        new_vertex.subread_count += v.subread_count;
    }
    processNewOuterVertex(new_vertex);
    new_vertex.covering_read_count = path.frontEdge().getFinish().covering_read_count -
        new_vertex.incFront().read_tail_count + path.frontEdge().read_tail_count;
    Edge &inc = new_vertex.incFront();
    inc.min_equivalent_size = path.frontEdge().min_equivalent_size;
    inc.outgoing_read_count = path.frontEdge().outgoing_read_count;
}

void ag::AlignedReadStatisticsTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    for (auto &rec : resolution) {
        Vertex &new_vertex = *rec.first;
        processNewOuterVertex(new_vertex);
        // new_vertex.covering_read_count = rec.second.getSupport() - new_vertex.incFront().read_tail_count -
        //     new_vertex.front().rc().read_tail_count + new_vertex.subread_count;
        //Above is the correct formula. Its calculation is split between processing of core and core.rc() because
        //processNewOuterVertex fills read_tail_count only for incoming edges ofy new vertices.
        new_vertex.covering_read_count += rec.second.getSupport() - new_vertex.incFront().read_tail_count + new_vertex.subread_count;
        new_vertex.rc().covering_read_count -= new_vertex.incFront().read_tail_count;
    }
    if (core.inDeg() == 1 && core.outDeg() != 1) {
        core.incFrontVertex().subread_length += core.subread_length;
        core.incFrontVertex().subread_count += core.subread_count;
    }
    if (core.inDeg() > 1 && core.outDeg() == 1) {
        core.frontVertex().subread_length += core.subread_length;
        core.frontVertex().subread_count += core.subread_count;
    }
    std::unordered_map<EdgeId, size_t> in_deg;
    for (Edge &edge : core)
        in_deg[edge.getId()] = 0;
    for (std::pair<const ObjectId<Vertex>, InOutEdgePair> it: resolution) {
        in_deg[it.second.outgoing().getId()]++;
    }
    for (std::pair<const ObjectId<Vertex>, InOutEdgePair> it: resolution) {
        Edge &new_inc = it.first->incFront();
        Edge &old_out = it.second.outgoing();
        if (in_deg[it.second.outgoing().getId()] == 1) {
            new_inc.min_equivalent_size = old_out.min_equivalent_size;
            new_inc.outgoing_read_count = old_out.outgoing_read_count;
        } else {
            new_inc.min_equivalent_size = core.size() + 2;
            new_inc.outgoing_read_count = it.second.getSupport();
            VERIFY(it.second.getSupport() == suffix_tracker->getSuffixRecord(it.second.incoming()).countStartsWith(GraphPath(it.second.outgoing())));
        }
    }
}

ag::CoverageSamplingTracker::CoverageSamplingTracker(ag::AssemblyGraph &graph, ag::AlignedReadStorage &storage,
                                                       size_t k, size_t threads) :
        ag::AlignedReadStorageListener(storage, "CoverageSamplingTracker"),
        ag::ResolutionListener(graph, "CoverageSamplingTracker"), k(k) {
    fillLengthHistogram(storage, threads);
    total_bases = windowCapacity(1);
    kpomer_multiplier = multiplier(k + 1);
}

void ag::CoverageSamplingTracker::fillLengthHistogram(ag::AlignedReadStorage &storage, size_t threads) {
    omp_set_num_threads(threads);
    std::vector<std::vector<size_t>> thread_histograms(threads);
#pragma omp parallel for default(none) shared(storage, thread_histograms) schedule(dynamic, 100)
    for (size_t i = 0; i < storage.size(); i++) {
        if (!storage[i].valid())
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

double ag::CoverageSamplingTracker::windowCapacity(size_t s) const {
    if (s >= window_capacity.size())
        return 0;
    return window_capacity[s];
}

void ag::CoverageSamplingTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    VERIFY(e.getStartSize() == k);
//        Converting DBG coverage into k+1-mer samples
    v.coverage_info.addKpomerChunk(0, v.size(), e.intCov(), k, kpomer_multiplier);
}

void ag::CoverageSamplingTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    if (!core.isJunction()) {
        resolution.begin()->first->coverage_info += std::move(core.coverage_info.shift(core.rc().front().truncSize()));
    } else if (core.inDeg() == 1) {
        core.incFrontVertex().coverage_info += std::move(core.coverage_info.shift(core.rc().front().truncSize()));
    } else if (core.outDeg() == 1) {
        core.frontVertex().coverage_info += std::move(core.coverage_info);
    } else {
        for (const auto &rec : resolution) {
            Vertex &new_vertex = *rec.first;
            const InOutEdgePair &pair = rec.second;
            size_t start = rec.second.incoming().rc().truncSize() - 1;
            size_t finish = rec.second.incoming().fullSize() + 1;
            new_vertex.coverage_info.addResolutionSample(start, finish, pair.getSupport(), multiplier(finish - start));
        }
    }
}

void ag::CoverageSamplingTracker::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    std::vector<Vertex *> inner;
    for (Vertex &v : path.innerVertices())
        inner.push_back(&v);
    size_t offset = path.getStart().size();
    size_t idx = 0;
    for (Edge &e : path.edges()) {
        offset += e.truncSize();
        if (idx == inner.size())
            continue;
        Vertex &v = *inner[idx];
        VERIFY(offset >= v.size());
        v.coverage_info.shift(__int64_t(offset) - __int64_t(v.size()));
        new_vertex.coverage_info += std::move(v.coverage_info);
        idx++;
    }
}

void ag::CoverageSamplingTracker::adjustSupport(Vertex &v, size_t left, size_t right, __int64_t mult) {
    CoverageSamples &info = v.coverage_info;
//        Samples are stored in position order regardless of kind. If the affected range reaches the
//        vertex's own end (right == v.size()), scanning back-to-front lets us stop as soon as we drop
//        below left instead of walking through the whole unaffected prefix; otherwise scanning
//        front-to-back stops as soon as we reach right. Either direction is correct — this only picks the
//        one that can abort earliest.
    bool backward = right == v.size();
    IterableStorage<CoverageSamples::SampleIterator> samples = backward ? info.getRSamples() : info.getSamples();
    for (CoverageSamples::SampleView view : samples) {
        if (backward ? view.finish() <= left : view.start() >= right)
            break;
//        unit_size is the length of a single witness: k+1 nucleotides for a kpomer k+1-mer, or the whole
//        sample for a resolution witness (one atomic, all-or-nothing unit). [overlap_left, overlap_right)
//        is the honest overlap between the sample's own window and the read's covered range; the number
//        of unit_size-long witnesses that fit entirely inside it is its length minus unit_size plus 1.
//        Each such witness is scaled by multiplier(unit_size) -- see the class declaration -- so votes
//        for samples of different lengths remain comparable/combinable.
        bool is_kpomer = view.sample.type == CoverageSamples::SampleType::kpomer;
        size_t unit_size = is_kpomer ? k + 1 : view.size();
        size_t overlap_left = std::max(view.start(), left);
        size_t overlap_right = std::min(view.finish(), right);
        if (overlap_right > overlap_left + unit_size - 1) {
            double m = is_kpomer ? kpomer_multiplier : multiplier(unit_size);
            __int64_t raw_delta = __int64_t(overlap_right - overlap_left - unit_size + 1) * mult;
            info.support += double(raw_delta) * m;
            view.sample.raw_support += raw_delta;
            info.raw_support += raw_delta;
        }
    }
}

void ag::CoverageSamplingTracker::processPath(const GraphPath &path, __int64_t mult) {
    if (!path.valid())
        return;
    if (path.empty()) {
        Vertex &v = path.getStart();
        adjustSupport(v, path.leftCut(), v.size() - path.rightCut(), mult);
        return;
    }
    adjustSupport(path.getStart(), path.leftCut(), path.getStart().size(), mult);
    for (Vertex &v : path.innerVertices())
        adjustSupport(v, 0, v.size(), mult);
    adjustSupport(path.getFinish(), 0, path.getFinish().size() - path.rightCut(), mult);
}

void ag::CoverageSamplingTracker::fireAddRead(const ag::AlignedRead &read) {
    VERIFY(false);
}

void ag::CoverageSamplingTracker::fireRerouteRead(ag::AlignedRead &read) {
    processPath(read.getPath(), -1);
    processPath(read.getPath().RC(), -1);
    processPath(read.getCorrected(), 1);
    processPath(read.getCorrected().RC(), 1);
}

void ag::CoverageSamplingTracker::fireInvalidateRead(ag::AlignedRead &read) {
    processPath(read.getPath(), -1);
    processPath(read.getPath().RC(), -1);
}
