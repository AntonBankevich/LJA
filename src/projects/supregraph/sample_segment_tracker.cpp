#include "sample_segment_tracker.hpp"

using namespace spg;
using namespace ag;

// New windows are introduced only if this is applied during the initial conversion to supregraph
void SampleSegmentTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    if (e.getStart().size() == k && e.getFinish().size() == k)
        v.coverage_info.addKpomerChunk(0, v.size(), 0, k, 0.0);
}

void SampleSegmentTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    VERIFY(core.inDeg() > 1);
    VERIFY(core.outDeg() > 1);
    for (const auto &rec : resolution) {
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
    new_vertex.coverage_info = std::move(path.getStart().coverage_info);
    __int64_t shift = 0;
    for (Edge &e : path.edges()) {
        shift += e.rc().truncSize();
        e.getFinish().coverage_info.shift(shift);
        new_vertex.coverage_info += std::move(e.getFinish().coverage_info);
    }
}
