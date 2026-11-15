#include "sample_segment_tracker.hpp"
#include "assembly_graph.hpp"

using namespace ag;

// New windows are introduced only if this is applied during the initial conversion to supregraph
void SampleSegmentTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    size_t k = getGraph().getK();
    if (e.getStart().size() == k && e.getFinish().size() == k)
        v.coverage_info.addKpomerChunk(0, v.size(), 0, k, 0.0);
}

void SampleSegmentTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    VERIFY(core.inDeg() > 1);
    VERIFY(core.outDeg() > 1);
    for (const auto &rec : resolution) {
//        With a single partner on either side, the strings one longer than core are already in the adjacent vertices.
        if (resolution.connections(rec.second.incoming()) < 2 || resolution.connections(rec.second.outgoing()) < 2)
            continue;
        Vertex &new_vertex = *rec.first;
        VERIFY(new_vertex.size() >= core.size() + 2);
        VERIFY(rec.second.incoming().rc().truncSize() >= new_vertex.incFront().rc().truncSize() + 1)
        VERIFY(rec.second.outgoing().truncSize() >= new_vertex.front().truncSize() + 1);
        VERIFY(new_vertex.size() - (rec.second.outgoing().truncSize() - new_vertex.front().truncSize() - 1) == rec.second.incoming().rc().truncSize() - new_vertex.incFront().rc().truncSize() - 1 + core.size() + 2)
        new_vertex.coverage_info.addResolutionSample(rec.second.incoming().rc().truncSize() - new_vertex.incFront().rc().truncSize() - 1,
            new_vertex.size() - (rec.second.outgoing().truncSize() - new_vertex.front().truncSize() - 1),
            0, 0);
    }
}

void SampleSegmentTracker::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    // Only inner vertices are absorbed into the merged vertex; move their coverage samples to the
    // matching offset. The start and end of the path survive and keep their own coverage.
    __int64_t shift = 0;
    for (Edge &e : path.edges()) {
        shift += e.rc().truncSize();
        if (e.getFinish() == path.getFinish())
            break;
        e.getFinish().coverage_info.shift(shift);
        new_vertex.coverage_info += std::move(e.getFinish().coverage_info);
    }
}

void SampleSegmentTracker::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    for (Vertex &v : path.innerVertices())
        VERIFY(v.coverage_info.empty());
}

// Inner vertices of a prefix run are prefixes of its finish and their segments come before the finish's own; inner
// vertices of a suffix run (the mirrored call) are suffixes of its start and their segments come after the start's.
void SampleSegmentTracker::fireMergePrefixPath(const RAGraphPath &path, Edge &new_edge) {
    if (path.frontEdge().isPrefix()) {
        CoverageSamples inner;
        for (Vertex &v : path.innerVertices())
            inner += std::move(v.coverage_info);
        path.getFinish().coverage_info.insertFront(std::move(inner));
    } else {
        Vertex &start = path.getStart();
        for (Vertex &v : path.innerVertices()) {
            v.coverage_info.shift(__int64_t(start.size()) - __int64_t(v.size()));
            start.coverage_info += std::move(v.coverage_info);
        }
    }
}

// Forward call: core is a suffix of its single in-neighbour u, so its segments go after u's. Mirrored call: core is a
// prefix of its single out-neighbour, so its segments go before that one's.
void SampleSegmentTracker::fireFlipVertex(Edge &edge, const std::unordered_map<EdgeId, EdgeId> &extensions) {
    if (edge.isSuffix()) {
        Vertex &core = edge.getFinish();
        Vertex &u = edge.getStart();
        core.coverage_info.shift(__int64_t(u.size()) - __int64_t(core.size()));
        u.coverage_info += std::move(core.coverage_info);
    } else {
        edge.getFinish().coverage_info.insertFront(std::move(edge.getStart().coverage_info));
    }
}

bool SampleSegmentTracker::checkConsistency(std::ostream &os) {
    const size_t k = getGraph().getK();
    size_t errors = 0;
    typedef std::tuple<__int64_t, __int64_t, CoverageSamples::SampleType, size_t> Rec;
    auto records = [](Vertex &v) {
        std::vector<Rec> res;
        for (CoverageSamples::SampleView s : v.coverage_info.getSamples())
            res.emplace_back(__int64_t(s.start()), __int64_t(s.finish()), s.sample.type, s.sample.raw_support);
        return res;
    };
    for (Vertex &v : getGraph().vertices()) {
        auto report = [&]() -> std::ostream & {
            errors++;
            return os << "Sample segment inconsistency at vertex " << v.getInnerId() << ": ";
        };
        const __int64_t size = v.size();
        std::vector<Rec> recs = records(v);
        if (v.size() == k && !recs.empty())
            report() << "a k-mer has samples\n";
        double weight = 0;
        size_t raw = 0;
        for (size_t i = 0; i < recs.size(); i++) {
            __int64_t s = std::get<0>(recs[i]), f = std::get<1>(recs[i]);
            bool kpomer = std::get<2>(recs[i]) == CoverageSamples::SampleType::kpomer;
            if (s < 0 || s >= f || f > size || (kpomer && f - s < __int64_t(k) + 1))
                report() << "sample [" << s << ", " << f << ") does not fit\n";
            if (i > 0 && (s < std::get<0>(recs[i - 1]) || f < std::get<1>(recs[i - 1]) ||
                          (s == std::get<0>(recs[i - 1]) && f == std::get<1>(recs[i - 1]))))
                report() << "sample [" << s << ", " << f << ") is out of order\n";
            weight += kpomer ? double(f - s - __int64_t(k)) : 1;
            raw += std::get<3>(recs[i]);
        }
        if (raw != v.coverage_info.raw_support || std::abs(weight - v.coverage_info.weight) > 1e-6)
            report() << "weight/raw_support differ from the sums over the samples\n";
        std::vector<Rec> rc = records(v.rc());
        bool mirror = rc.size() == recs.size();
        for (size_t i = 0; mirror && i < recs.size(); i++) {
            const Rec &r = recs[recs.size() - 1 - i];
            mirror = rc[i] == Rec(size - std::get<1>(r), size - std::get<0>(r), std::get<2>(r), std::get<3>(r));
        }
        if (!mirror)
            report() << "samples of the rc are not the mirror\n";
    }
    if (errors > 0)
        os << "Sample segment inconsistency: " << errors << " violation(s) found" << std::endl;
    return errors == 0;
}
