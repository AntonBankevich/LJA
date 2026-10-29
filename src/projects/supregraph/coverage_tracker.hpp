#pragma once

#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/data_structures/read_alignment_storage.hpp"

namespace spg {
    using ag::AssemblyGraph;
    using ag::AlignedReadStorage;
    using ag::AlignedRead;
    using ag::Vertex;
    using ag::GraphPath;

//    Computes per-vertex Supregraph coverage purely from reads: every number here comes from walking a
//    read's path and crediting whichever sample segments (anchors) it covers (see SampleSegmentTracker,
//    which owns those segments' identity/structure but never writes their support values) -- with one
//    exception: a resolveVertex's genuinely ambiguous case (neither inDeg nor outDeg is 1) creates brand
//    new vertices whose only incoming/outgoing edges are new too, so no read can be mid-reroute across
//    them yet when fireResolveVertex fires -- there is nothing for fireRerouteRead to replay against.
//    SampleSegmentTracker already created each such vertex's one resolution-sample window at zero support;
//    fireResolveVertex here patches it to the real value directly from the resolution result's
//    InOutEdgePair::getSupport() (an exact count of reads already computed during resolution planning, not
//    a guess), via the same adjustSupport used for read replay. The other three resolveVertex cases
//    transfer core's own (already-real) coverage_info onto an existing vertex and need no further touch.
    class CoverageTracker : public ag::AlignedReadStorageListener, public ag::ResolutionListener {
    private:
        AssemblyGraph &graph;
        ag::AlignedReadStorage *storage;
        size_t k;
        std::vector<double> window_capacity;
        double total_bases = 0;
//      Distribution of read lengths is used to inform the multiplier to be used for sampling contribution to coverage
        double windowCapacity(size_t s) const;
        double multiplier(size_t s) const {double wc = windowCapacity(s); return wc == 0 ? 0.0 : total_bases / wc;}
        double kpomer_multiplier;
//        One-pass, thread-parallel scan of every read's full aligned length, building cnt_ge/sum_ge.
        void fillLengthHistogram(ag::AlignedReadStorage &storage, size_t threads);
//       Similar to DBG path processing for coverage update, but instead of k+1-mers, segments of variable
//       size are considered and processed individually.
        void processPath(const GraphPath &path, __int64_t mult);
        void adjustSingleSupport(ag::CoverageSamples &info, ag::CoverageSamples::SampleView &view, size_t left,
                                 size_t right, __int64_t mult);
//        Sample information is stored in one of two types of records: k+1-mer chunk record or a record for
//        a larger segment. This method chooses how to process new information properly.
        void adjustSupport(Vertex &v, size_t left, size_t right, __int64_t mult);
    public:
        CoverageTracker(AssemblyGraph &graph, ag::AlignedReadStorage &storage, size_t k, logging::Logger &logger, size_t threads);

//        Explicit, on-demand replay of every read currently in storage. Not called automatically from the
//        constructor: at construction time (before DBG->SPG conversion) there are no sample windows yet
//        for reads to credit. Call once after conversion is complete, then let incremental listening
//        (fireAddRead/fireRerouteRead/fireInvalidateRead) take over for everything after that.
        void fillFromStorage(logging::Logger &logger, size_t threads);

        void fireAddRead(const ag::AlignedRead &read) override;
        void fireRerouteRead(ag::AlignedRead &read) override;
        void fireInvalidateRead(ag::AlignedRead &read) override;

        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;
        void fireEdgeToSupreVertex(Vertex &v, ag::Edge &e) override {};

//        Subtracts every currently-stored read's contribution from every vertex's sample coverage (the
//        exact inverse of fireAddRead/fillFromStorage), verifies every vertex's support/raw_support
//        returns to the zero baseline SampleSegmentTracker always creates segments at, then adds the
//        reads back. A non-zero residual means the stored coverage no longer corresponds to what the
//        current read paths actually support -- either a read's contribution was double-counted/missed,
//        or coverage was seeded some other way than replaying reads. Non-destructive: restores state
//        before returning either way, like SuffixTracker::fireCheckConsistency.
        bool fireCheckConsistency() override;

        std::function<std::string(const Vertex &)> getVertexLabeler() const;
    };
}
