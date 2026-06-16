#include "aligned_read_statistics_tracker.hpp"

void ag::AlignedReadStatisticsTracker::addPath(const GraphPath &path, __int64_t mult) {
    if (!path.valid())
        return;
    if (path.empty()) {
        path.getStart().subread_length += path.len() * mult;
        path.getStart().subread_count += mult;
    } else {
        for (Vertex &v : path.innerVertices()) {
            v.covering_read_count+=mult;
        }
        VERIFY(!path.frontEdge().isPrefix());
        VERIFY(!path.backEdge().isSuffix());
        path.backEdge().read_tail_length += (path.backEdge().fullSize() - path.rightCut()) * mult;
        path.backEdge().read_tail_count += mult;
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

ag::AlignedReadStatisticsTracker::AlignedReadStatisticsTracker(logging::Logger &logger, size_t threads,
                                                               ag::ResolutionFire &graph, ag::AlignedReadStorage &storage, ag::SuffixTracker &suffix_tracker):
    ag::AlignedReadStorageListener(storage, "AlignedReadStatisticsTracker"),
    ag::ResolutionListener(graph, "AlignedReadStatisticsTracker"), storage(&storage),
    suffix_tracker(&suffix_tracker) {
    fillFromStorage(logger, threads);
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
}

void ag::AlignedReadStatisticsTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    processNewOuterVertex(v);
}

void ag::AlignedReadStatisticsTracker::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    for (Vertex &v : path.innerVertices()) {
        new_vertex.subread_length += v.subread_length;
        new_vertex.subread_count += v.subread_count;
    }
    processNewOuterVertex(new_vertex);
}

void ag::AlignedReadStatisticsTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    for (Vertex &new_vertex : resolution.newVertices()) {
        processNewOuterVertex(new_vertex);
    }
    if (core.inDeg() == 1 && core.outDeg() != 1) {
        core.rc().front().getFinish().rc().subread_length += core.subread_length;
        core.rc().front().getFinish().rc().subread_count += core.subread_count;
    }
    if (core.inDeg() > 1 && core.outDeg() == 1) {
        core.front().getFinish().subread_length += core.subread_length;
        core.front().getFinish().subread_count += core.subread_count;
    }
}
