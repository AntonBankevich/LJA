#include "aligned_read_statistics_tracker.hpp"

void ag::AlignedReadStatisticsTracker::addPath(const GraphPath &path, __int64_t mult) {
    if (!path.valid())
        return;
    if (path.empty()) {
        Vertex &v = path.getStart();
        v.subread_length += path.len() * mult;
        v.subread_count += mult;
        // if (v.outDeg() == 1 && v.front().isSuffix()) {
        //     if (path.rightCut() <= v.front().getFinish().size()) {
        //         VERIFY(path.leftCut() < v.size() - v.front().getFinish().size());
        //         v.front().rc().outgoing_read_count += mult;
        //     }
        // }
    } else {
        for (Vertex &v : path.innerVertices()) {
            v.covering_read_count+=mult;
        }
        VERIFY(!path.frontEdge().isPrefix());
        VERIFY(!path.backEdge().isSuffix());
        path.backEdge().read_tail_length += (path.backEdge().fullSize() - path.rightCut()) * mult;
        path.backEdge().read_tail_count += mult;
        path.frontEdge().outgoing_read_count -= mult;
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
    Edge &inc = v.rc().front().rc();
    inc.min_equivalent_size = e.min_equivalent_size;
    inc.outgoing_read_count = e.outgoing_read_count;
}

void ag::AlignedReadStatisticsTracker::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    for (Vertex &v : path.innerVertices()) {
        new_vertex.subread_length += v.subread_length;
        new_vertex.subread_count += v.subread_count;
    }
    processNewOuterVertex(new_vertex);
    Edge &inc = new_vertex.rc().front().rc();
    inc.min_equivalent_size = path.frontEdge().min_equivalent_size;
    inc.outgoing_read_count = path.frontEdge().outgoing_read_count;
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
    std::unordered_map<EdgeId, size_t> in_deg;
    for (Edge &edge : core)
        in_deg[edge.getId()] = 0;
    for (std::pair<const ObjectId<Vertex>, InOutEdgePair> it: resolution) {
        in_deg[it.second.outgoing().getId()]++;
    }
    for (std::pair<const ObjectId<Vertex>, InOutEdgePair> it: resolution) {
        Edge &new_inc = it.first->rc().front().rc();
        if (in_deg[it.second.outgoing().getId()] == 1) {
            new_inc.min_equivalent_size = it.second.outgoing().min_equivalent_size;
            new_inc.outgoing_read_count = it.second.outgoing().outgoing_read_count;
        } else {
            new_inc.min_equivalent_size = core.size() + 2;
            suffix_tracker->getSuffixRecord(it.second.incoming()).countStartsWith(GraphPath(it.second.outgoing()));
        }
    }
}
