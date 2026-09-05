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
