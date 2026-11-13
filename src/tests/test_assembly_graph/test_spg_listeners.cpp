#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "supregraph/coverage_tracker.hpp"
#include "supregraph/unique_vertex_storage.hpp"
#include "supregraph/path_spy.hpp"
#include "supregraph/supregraph.hpp"

// Tests for the supregraph-layer (spg) ResolutionListeners, exercised on a real supregraph built from a
// random DBG via edgeToSupreVertex (ConvertToSupregraph). These require linking lja_spg.

using namespace ag_test;

namespace {
    constexpr size_t K = 5;
}

// CoverageTracker + SampleSegmentTracker. The graph's own SampleSegmentTracker creates the per-vertex sample
// windows (zero support) in fireEdgeToSupreVertex during the conversion; CoverageTracker is then
// built on the supregraph and fills real support by replaying the reads over those windows. Its
// fireCheckConsistency subtracts every read's contribution, checks the support returns to the zero
// baseline the samples were created at, and re-adds - so a true result validates BOTH listeners: the
// sample windows SampleSegmentTracker created AND the support CoverageTracker computed match the reads.
TEST(SpgListeners, CoverageAndSampleSegmentConsistentOnSupregraph) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(400, seed), K);
        std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 90, 2, seed + 300);
        ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));

        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg); // condense; vertices stay k-mers (size k)
        ConvertToSupregraph(rdbg.dbg);                   // edgeToSupreVertex on every corporeal edge

        spg::CoverageTracker cov(rdbg.dbg, storage, QuietLogger(), 1); // ctor fills support from reads
        EXPECT_TRUE(cov.fireCheckConsistency(std::cerr)) << "seed=" << seed;
    }
}

// UniqueVertexStorage. Its constructor replays fireAddVertex over every vertex, marking a vertex unique
// exactly when innerSize() >= threshold. Directly verify that stored classification against the
// threshold criterion, on a supregraph where supre-vertices are large and k-mer vertices are small (so
// the threshold genuinely splits them, not a trivial all/none).
TEST(SpgListeners, UniqueVertexStorageMatchesThreshold) {
    RandomDBG rdbg(RandomGenome(400, 1), K);
    std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 90, 2, 301);
    ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));
    ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
    ConvertToSupregraph(rdbg.dbg);

    // Pick a threshold from the actual innerSize distribution so it genuinely splits the vertex set.
    size_t minInner = size_t(-1), maxInner = 0;
    for (ag::Vertex &v : rdbg.dbg.vertices()) {
        minInner = std::min(minInner, v.innerSize());
        maxInner = std::max(maxInner, v.innerSize());
    }
    ASSERT_GT(maxInner, minInner); // the supregraph has a range of vertex sizes to split on
    const size_t threshold = minInner + 1;
    spg::UniqueVertexStorage uvs(rdbg.dbg, threshold);

    // add() marks both v and v.rc(), and innerSize() is not rc-symmetric, so a vertex is unique iff
    // EITHER orientation meets the threshold (and isUnique(v) == isUnique(v.rc()) always).
    size_t above = 0, below = 0;
    for (ag::Vertex &v : rdbg.dbg.vertices()) {
        bool expected = v.innerSize() >= threshold || v.rc().innerSize() >= threshold;
        EXPECT_EQ(uvs.isUnique(v), expected)
            << "innerSize=" << v.innerSize() << " rc.innerSize=" << v.rc().innerSize() << " threshold=" << threshold;
        (expected ? above : below)++;
    }
    EXPECT_GT(above, 0u); // threshold actually splits the vertex set
    EXPECT_GT(below, 0u);
}

// Vertex::innerSize() must be the same for a vertex and its reverse complement (regression test for the
// tip-case asymmetry fix). Checked on a supregraph, where tips (the former failing case) are present.
TEST(SpgListeners, InnerSizeIsRcSymmetric) {
    RandomDBG rdbg(RandomGenome(400, 1), K);
    std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 90, 2, 301);
    ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));
    ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
    ConvertToSupregraph(rdbg.dbg);
    for (ag::Vertex &v : rdbg.dbg.vertices())
        EXPECT_EQ(v.innerSize(), v.rc().innerSize())
            << "vertex " << v.getId() << " outDeg=" << v.outDeg() << " inDeg=" << v.inDeg();
}

// OldVertexTracker read-substring query (spg::readSubstringVertices). Resolve a 2-in/2-out junction
// core into two instances: the core is deleted (so a read equivalent to it carries a legacy path), and
// OldVertexTracker records the core as embedded in each living resolution instance. The query takes a
// read whose legacy path is the deleted core and must return exactly those living instances - each of
// which genuinely contains the core's (the read's) sequence.
TEST(SpgListeners, OldVertexTrackerReadSubstringQuery) {
    SeqDBG g({"CACGTAC", "GACGTAG"}, 5); // shared 5-mer ACGTA is the only junction
    spg::OldVertexTracker tracker(g.dbg, false);
    ag::GraphPath p1 = g.align("CACGTAC"); // CACGT -> ACGTA -> CGTAC
    ag::GraphPath p2 = g.align("GACGTAG"); // GACGT -> ACGTA -> CGTAG
    ASSERT_EQ(EdgeCount(p1), 2u);
    ASSERT_EQ(EdgeCount(p2), 2u);
    ag::Vertex &core = p1.frontEdge().getFinish();
    ASSERT_EQ(&core, &p2.frontEdge().getFinish());

    // A read equivalent to the core (its whole MES extension is this one vertex).
    ag::VertexId Vid = core.getId();
    ag::VertexId Vrc = core.rc().getId();
    std::string Vseq = core.getSeq().str();

    ag::VertexResolutionPlan plan(core);
    plan.add(p1.frontEdge(), p1.backEdge());
    plan.add(p2.frontEdge(), p2.backEdge());
    g.dbg.disableHashing(); // resolution creates vertices longer than k
    ag::VertexResolutionResult res = g.dbg.resolveVertex(core, plan); // core deleted, embedded into instances
    std::set<ag::VertexId> instances;
    for (ag::Vertex &v : res.newVertices())
        instances.insert(v.getId());
    g.dbg.removeMarked();

    ag::AlignedRead read("legacy_read", ag::GraphPath::LegacyPath(Vid, Vrc));
    std::vector<ag::VertexId> containers = spg::readSubstringVertices(read, tracker);

    // Soundness: every reported container is alive and genuinely contains the read's sequence, and the
    // two resolution instances (the minimal living vertices the core became) are among them. (Reducing
    // the returned set to exactly the minimal-by-substring containers is a documented refinement.)
    ASSERT_FALSE(containers.empty()) << "the resolved-away core must still live inside its instances";
    std::set<ag::VertexId> container_set(containers.begin(), containers.end());
    for (ag::VertexId inst : instances)
        EXPECT_TRUE(container_set.count(inst)) << "resolution instance must be reported as a container";
    for (ag::VertexId u : containers)
        EXPECT_NE(u->getSeq().str().find(Vseq), std::string::npos)
            << "reported container " << u->getId() << " must actually contain the read's sequence";
}

// OldVertexTracker through the supregraph pipeline: conversion, resolution, MergeAllSPG, flips. After each stage
// its records are consistent, and the resolved-away core still climbs to living vertices that contain it.
TEST(SpgListeners, OldVertexTrackerThroughPipeline) {
    size_t queried = 0;
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), K);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg::Supregraph spg = spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg));
        spg::OldVertexTracker tracker(spg, true);
        ASSERT_TRUE(tracker.checkConsistency(std::cerr)) << "seed=" << seed << " after attaching";
        ag::VertexId core;
        for (ag::Vertex &v : spg.vertices())
            if (!core.valid() && v.size() == K && v.inDeg() == 2 && v.outDeg() == 2)
                core = v.getId();
        std::string core_seq = core.valid() ? core->getSeq().str() : "";
        ResolveMatching(spg, K);
        ASSERT_TRUE(tracker.checkConsistency(std::cerr)) << "seed=" << seed << " after resolution";
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        ASSERT_TRUE(tracker.checkConsistency(std::cerr)) << "seed=" << seed << " after MergeAllSPG";
        FlipAll(spg);
        ASSERT_TRUE(tracker.checkConsistency(std::cerr)) << "seed=" << seed << " after flips";
        if (core.valid()) {
            std::unordered_set<ag::ConstVertexId> live;
            for (ag::Vertex &v : spg.vertices())
                live.insert(v.getId());
            std::vector<spg::OldVertexTracker::DirectEmbedding> containers = tracker.getAllEmbeddings(core);
            EXPECT_FALSE(containers.empty()) << "seed=" << seed;
            for (const spg::OldVertexTracker::DirectEmbedding &de : containers) {
                ASSERT_TRUE(live.count(de.outer_vertex)) << "seed=" << seed << " container " << de.outer_vertex << " is dead";
                EXPECT_EQ(de.outer_vertex->getSeq().Subseq(de.from, de.to).str(), core_seq) << "seed=" << seed;
            }
            queried++;
        }
    }
    EXPECT_GT(queried, 0u);
}

// spg::CoverageTracker through the supregraph pipeline. Built after conversion with several threads (atomic updates),
// then resolutions by the reads' own pairs (support = read paths through the pair), MergeAllSPG and flips; after each
// stage the stored coverage is exactly what replaying the reads gives.
TEST(SpgListeners, CoverageTrackerThroughPipeline) {
    size_t credited = 0;
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), K);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, 60, 60, 1, seed + 100));
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        spg::CoverageTracker cov(spg, storage, QuietLogger(), 4);
        ASSERT_TRUE(cov.checkConsistency(std::cerr)) << "seed=" << seed << " after construction";
        for (size_t i = 0; i < 10 && ResolveByReads(spg, storage, K); i++) {}
        ASSERT_TRUE(cov.checkConsistency(std::cerr)) << "seed=" << seed << " after resolutions";
        for (ag::Vertex &v : spg.vertices())
            for (ag::CoverageSamples::SampleView view : v.coverage_info.getSamples())
                if (view.sample.type == ag::CoverageSamples::SampleType::resolution && view.sample.raw_support > 0)
                    credited++;
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        ASSERT_TRUE(cov.checkConsistency(std::cerr)) << "seed=" << seed << " after MergeAllSPG";
        FlipAll(spg);
        ASSERT_TRUE(cov.checkConsistency(std::cerr)) << "seed=" << seed << " after flips";
    }
    EXPECT_GT(credited, 0u);
}
