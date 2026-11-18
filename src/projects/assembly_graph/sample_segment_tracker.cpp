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
//        With a single partner on either side, the strings one longer than core are already in the adjacent vertices: a
//        1-to-1 pair's new vertex gets no sample segment.
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

// Every removed vertex's segments go to the dominant vertex D at its offset in it. Samples are kept in path order: those
// of removed vertices before D in P go in front of D's own, the others after them. As in fireEdgeToSupreVertex, a new D
// spelling a single edge between two k-mers (only on a DBG) gets its (k+1)-mer window.
void SampleSegmentTracker::fireReplacePath(const PathReplacement &rep) {
    const RAGraphPath &old_path = rep.oldPath();
    Vertex &dominant = rep.getDominant();
    CoverageSamples before;
    for (size_t i = 0; i <= old_path.size(); i++) {
        if (rep.survives(i))
            continue;
        CoverageSamples &samples = old_path.getVertex(i).coverage_info;
        VERIFY(rep.oldPosition(i) >= rep.dominantPosition());
        samples.shift(__int64_t(rep.oldPosition(i) - rep.dominantPosition()));
        if (rep.dominantOldIndex() != size_t(-1) && i < rep.dominantOldIndex())
            before += std::move(samples);
        else
            dominant.coverage_info += std::move(samples);
    }
    dominant.coverage_info.insertFront(std::move(before));
    const size_t k = getGraph().getK();
    if (rep.dominantOldIndex() == size_t(-1) && old_path.size() == 1 && old_path.getStart().size() == k &&
            old_path.getFinish().size() == k)
        dominant.coverage_info.addKpomerChunk(0, dominant.size(), 0, k, 0.0);
}

void SampleSegmentTracker::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    for (Vertex &v : path.innerVertices())
        VERIFY(v.coverage_info.empty());
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
