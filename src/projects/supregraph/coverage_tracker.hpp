#pragma once

#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/data_structures/read_alignment_storage.hpp"
#include "spg_listener.hpp"

namespace spg {
    using ag::AssemblyGraph;
    using ag::AlignedReadStorage;
    using ag::AlignedRead;
    using ag::Vertex;
    using ag::GraphPath;

//    Per-vertex supregraph coverage from reads. SampleSegmentTracker owns the sample windows; this tracker owns their
//    support. Invariant: every window's raw support is what replaying every valid stored read (both orientations)
//    gives; vertex totals are the sums over its windows.
//    Contracts:
//    - Attach after DBG -> SPG conversion: windows created by edgeToSupreVertex between two k-mers can't be credited
//      (reads are indexed by start vertex only), so such a conversion is not allowed while attached.
//    - The support of a resolved pair (in, out) is the number of read paths, in either orientation, that pass
//      in -> out and reach at least one nucleotide past core on both sides. It is credited to the resolution window
//      SampleSegmentTracker creates on the new vertex.
//    Updates are atomic, so reads can be replayed concurrently.
    class CoverageTracker : public ag::AlignedReadStorageListener, public SPGListener {
    private:
        AssemblyGraph &graph;
        ag::AlignedReadStorage *storage;
        std::vector<double> window_capacity;
        double total_bases = 0;
//      Distribution of read lengths is used to inform the multiplier to be used for sampling contribution to coverage
        double windowCapacity(size_t s) const;
        double multiplier(size_t s) const {double wc = windowCapacity(s); return wc == 0 ? 0.0 : total_bases / wc;}
        double kpomer_multiplier;
//        One-pass, thread-parallel scan of every read's full aligned length, building cnt_ge/sum_ge.
        void fillLengthHistogram(ag::AlignedReadStorage &storage, size_t threads);
//       Similar to DBG path processing for coverage update, but instead of k+1-mers, segments of variable
//       size are considered and processed individually. Credits both the path and its rc.
        void processPath(const GraphPath &path, __int64_t mult);
//        adjustSupport on [left, right) of v and on the mirrored range of v.rc().
        void adjustBothStrands(Vertex &v, size_t left, size_t right, __int64_t mult);
        void adjustSingleSupport(ag::CoverageSamples &info, ag::CoverageSamples::SampleView &view, size_t left,
                                 size_t right, __int64_t mult);
//        Sample information is stored in one of two types of records: k+1-mer chunk record or a record for
//        a larger segment. This method chooses how to process new information properly.
        void adjustSupport(Vertex &v, size_t left, size_t right, __int64_t mult);
    public:
        CoverageTracker(AssemblyGraph &graph, ag::AlignedReadStorage &storage, logging::Logger &logger, size_t threads);

//        Replays every read currently in storage; the constructor calls it.
        void fillFromStorage(logging::Logger &logger, size_t threads);

        void fireAddRead(const ag::AlignedRead &read) override;
        void fireRerouteRead(ag::AlignedRead &read) override;
        void fireInvalidateRead(ag::AlignedRead &read) override;

        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;
        void fireEdgeToSupreVertex(Vertex &v, ag::Edge &e) override;

//        Checks the invariant: subtracts every stored read, expects every vertex back at zero support, then adds the
//        reads back (non-destructive). Single-threaded only.
        bool checkConsistency(std::ostream &os) override;
        bool fireCheckConsistency(std::ostream &os) override {return checkConsistency(os);}

        std::function<std::string(const Vertex &)> getVertexLabeler() const;
    };
}
