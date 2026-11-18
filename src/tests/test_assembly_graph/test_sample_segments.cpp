#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "supregraph/supregraph.hpp"

// Sample segments by definition: every sample segment s of a vertex v is a smallest string equivalent to v, where
// equivalence is with respect to all genomes the graph allows. Checked directly on small graphs.

using namespace ag_test;

namespace {
//    A position in the graph's sequence: v->getSeq()[i] if e is null, otherwise e->truncSeq()[i].
    struct Pos {
        ag::Vertex *v;
        ag::Edge *e;
        size_t i;
        char nucl() const { return "ACGT"[e == nullptr ? v->getSeq()[i] : e->truncSeq()[i]]; }
    };

    void AppendFromEnd(ag::Vertex &x, std::vector<Pos> &res) {
        for (ag::Edge &f : x) {
            if (f.truncSize() == 0)
                AppendFromEnd(f.getFinish(), res);
            else
                res.push_back({nullptr, &f, 0});
        }
    }

//    Positions of the next nucleotide on every genome continuing from p.
    std::vector<Pos> Next(const Pos &p) {
        if (p.e == nullptr && p.i + 1 < p.v->size())
            return {{p.v, nullptr, p.i + 1}};
        if (p.e != nullptr && p.i + 1 < p.e->truncSize())
            return {{nullptr, p.e, p.i + 1}};
        std::vector<Pos> res;
        AppendFromEnd(p.e == nullptr ? *p.v : p.e->getFinish(), res);
        return res;
    }

    bool Forced(Pos p, const std::string &rest) {
        for (char c : rest) {
            std::vector<Pos> next = Next(p);
            if (next.size() != 1 || next[0].nucl() != c)
                return false;
            p = next[0];
        }
        return true;
    }

    bool OccurrencesForced(const Pos &p, const std::string &s, size_t m, const std::string &rest) {
        if (p.nucl() != s[m])
            return true;
        if (m + 1 == s.size())
            return Forced(p, rest);
        for (const Pos &q : Next(p))
            if (!OccurrencesForced(q, s, m + 1, rest))
                return false;
        return true;
    }

//    Whether every occurrence of s in the graph is continued by rest on every genome the graph allows.
    bool ForcedRight(ag::AssemblyGraph &graph, const std::string &s, const std::string &rest) {
        for (ag::Vertex &v : graph.vertices()) {
            for (size_t i = 0; i < v.size(); i++)
                if (!OccurrencesForced({&v, nullptr, i}, s, 0, rest))
                    return false;
            for (ag::Edge &e : v)
                for (size_t i = 0; i < e.truncSize(); i++)
                    if (!OccurrencesForced({nullptr, &e, i}, s, 0, rest))
                        return false;
        }
        return true;
    }

    std::string RC(const std::string &s) { return (!Sequence(s)).str(); }

//    s = v[from, to) is equivalent to v.
    bool Equivalent(ag::AssemblyGraph &graph, const std::string &v, size_t from, size_t to) {
        std::string s = v.substr(from, to - from);
        return ForcedRight(graph, s, v.substr(to)) && ForcedRight(graph, RC(s), RC(v.substr(0, from)));
    }

//    For s equivalent to its vertex: no shorter substring is, i.e. neither s without its last nor without its first
//    nucleotide is equivalent to s.
    bool Minimal(ag::AssemblyGraph &graph, const std::string &s, size_t k) {
        return s.size() == k + 1 ||
               (!ForcedRight(graph, s.substr(0, s.size() - 1), s.substr(s.size() - 1)) &&
                !ForcedRight(graph, RC(s.substr(1)), RC(s.substr(0, 1))));
    }

//    Checks every sample segment of every vertex; reports failures into err.
    bool SegmentsAreMinimalEquivalent(ag::AssemblyGraph &graph, std::string &err) {
        const size_t k = graph.getK();
        std::stringstream ss;
        for (ag::Vertex &v : graph.vertices()) {
            std::string seq = v.getSeq().str();
            for (ag::CoverageSamples::SampleView sample : v.coverage_info.getSamples()) {
                std::vector<std::pair<size_t, size_t>> segments;
                if (sample.sample.type == ag::CoverageSamples::SampleType::kpomer)
                    for (size_t i = sample.start(); i + k + 1 <= sample.finish(); i++)
                        segments.emplace_back(i, i + k + 1);
                else
                    segments.emplace_back(sample.start(), sample.finish());
                for (const auto &seg : segments) {
                    if (!Equivalent(graph, seq, seg.first, seg.second))
                        ss << "vertex " << v.getInnerId() << ": segment [" << seg.first << ", " << seg.second
                           << ") is not equivalent to the vertex\n";
                    else if (!Minimal(graph, seq.substr(seg.first, seg.second - seg.first), k))
                        ss << "vertex " << v.getInnerId() << ": segment [" << seg.first << ", " << seg.second
                           << ") is not minimal\n";
                }
            }
        }
        err = ss.str();
        return err.empty();
    }
}

TEST(SampleSegments, MinimalEquivalentAfterConversion) {
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        ConvertToSupregraph(rdbg.dbg);
        std::string err;
        EXPECT_TRUE(SegmentsAreMinimalEquivalent(rdbg.dbg, err)) << "seed=" << seed << "\n" << err;
    }
}

// The checker is not vacuous: a segment that branches is reported as not equivalent, an equivalent segment that has an
// equivalent substring as not minimal.
TEST(SampleSegments, CheckerRejectsBadSegments) {
    const size_t k = 5;
    RandomDBG rdbg(RandomGenome(300, 1), k);
    ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
    ConvertToSupregraph(rdbg.dbg);
    ag::Vertex *branching = nullptr;
    ag::Vertex *long_vertex = nullptr;
    for (ag::Vertex &v : rdbg.dbg.vertices()) {
        if (v.size() == k || v.inDeg() != 1)
            continue;
        if (branching == nullptr && v.incFront().getStart().outDeg() > 1)
            branching = &v;
        if (long_vertex == nullptr && v.size() >= k + 4)
            long_vertex = &v;
    }
    ASSERT_NE(branching, nullptr);
    ASSERT_NE(long_vertex, nullptr);
    std::string err;

    branching->coverage_info.addResolutionSample(0, k, 0, 0);
    EXPECT_FALSE(SegmentsAreMinimalEquivalent(rdbg.dbg, err));
    EXPECT_NE(err.find("[0, 5) is not equivalent"), std::string::npos) << err;
    branching->coverage_info = ag::CoverageSamples();

    long_vertex->coverage_info.addResolutionSample(1, k + 3, 0, 0);
    EXPECT_FALSE(SegmentsAreMinimalEquivalent(rdbg.dbg, err));
    EXPECT_NE(err.find("[1, 8) is not minimal"), std::string::npos) << err;
}

// resolveVertex on a 2-in/2-out junction k-mer of a converted SPG, with every pair resolved (full) or with a one-to-one
// matching. Sample segments must satisfy the definition in the resulting graph.
TEST(SampleSegments, MinimalEquivalentAfterResolveVertex) {
    for (bool full : {true, false}) {
        for (uint32_t seed = 1; seed <= 5; seed++) {
            RandomDBG rdbg(RandomGenome(300, seed), 5);
            ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
            ConvertToSupregraph(rdbg.dbg);
            ag::Vertex *core = nullptr;
            for (ag::Vertex &v : rdbg.dbg.vertices())
                if (core == nullptr && v.size() == 5 && v.inDeg() == 2 && v.outDeg() == 2)
                    core = &v;
            if (core == nullptr)
                continue;
            std::vector<ag::Edge *> in, out;
            for (ag::Edge &e : core->incoming())
                in.push_back(&e);
            for (ag::Edge &e : *core)
                out.push_back(&e);
            ag::VertexResolutionPlan plan(*core);
            for (size_t i = 0; i < 2; i++)
                for (size_t j = 0; j < 2; j++)
                    if (full || i == j)
                        plan.add(*in[i], *out[j]);
            rdbg.dbg.resolveVertex(*core, plan);
            rdbg.dbg.removeMarked();
            std::string err;
            EXPECT_TRUE(SegmentsAreMinimalEquivalent(rdbg.dbg, err))
                                << (full ? "full" : "matching") << " seed=" << seed << "\n" << err;
        }
    }
}

TEST(SampleSegments, MinimalEquivalentAfterMergeAllSPG) {
    size_t merged = 0;
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg::Supregraph spg = spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg));
        if (!ResolveMatching(spg, 5))
            continue;
        size_t before = spg.size();
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        merged += before - spg.size();
        std::string err;
        EXPECT_TRUE(SegmentsAreMinimalEquivalent(spg, err)) << "seed=" << seed << "\n" << err;
        EXPECT_TRUE(spg.checkConsistency()) << "seed=" << seed;
    }
    EXPECT_GT(merged, 0u);
}

TEST(SampleSegments, MinimalEquivalentAfterFlips) {
    size_t total_flipped = 0;
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        ConvertToSupregraph(rdbg.dbg);
        ResolveMatching(rdbg.dbg, 5);
        std::vector<ag::VertexId> cores;
        for (ag::Vertex &v : rdbg.dbg.vertices())
            if (v.inDeg() == 1 && v.incFront().isSuffix() && v.outDeg() != 1)
                cores.push_back(v.getId());
        size_t flipped = 0;
        for (ag::VertexId vid : cores) {
            if (vid->marked() || vid->inDeg() != 1 || !vid->incFront().isSuffix() || vid->outDeg() == 1)
                continue;
            rdbg.dbg.flipVertexLockFree(vid->incFront());
            flipped++;
        }
        rdbg.dbg.removeMarked();
        std::string err;
        EXPECT_TRUE(SegmentsAreMinimalEquivalent(rdbg.dbg, err)) << "seed=" << seed << " flipped=" << flipped << "\n" << err;
        EXPECT_TRUE(rdbg.dbg.checkConsistency()) << "seed=" << seed;
        total_flipped += flipped;
    }
    EXPECT_GT(total_flipped, 0u);
}

TEST(SampleSegments, TrackerCheckRejectsBadSamples) {
    RandomDBG rdbg(RandomGenome(300, 1), 5);
    ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
    ConvertToSupregraph(rdbg.dbg);
    ASSERT_TRUE(rdbg.dbg.checkConsistency());
    ag::Vertex *kmer = nullptr;
    for (ag::Vertex &v : rdbg.dbg.vertices())
        if (kmer == nullptr && v.size() == 5)
            kmer = &v;
    ASSERT_NE(kmer, nullptr);
    kmer->coverage_info.addResolutionSample(0, 5, 0, 0);
    std::stringstream ss;
    EXPECT_FALSE(rdbg.dbg.checkConsistency(ss));
    EXPECT_NE(ss.str().find("a k-mer has samples"), std::string::npos) << ss.str();
    EXPECT_NE(ss.str().find("samples of the rc are not the mirror"), std::string::npos) << ss.str();
}

// flipVertex moves the core's samples to its single in-neighbour u (core is a suffix of u), after u's own; on the rc
// side before them. Checked against the expected coordinates and by the tracker's own check (order, rc mirror).
TEST(SampleSegments, FlipMovesSamples) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &u = graph.addVertex(Sequence("TTAACG"));
    ag::Vertex &core = graph.addVertex(Sequence("AACG"));
    ag::Vertex &w1 = graph.addVertex(Sequence("AACGA"));
    ag::Vertex &w2 = graph.addVertex(Sequence("AACGT"));
    graph.addSPEdge(u, core);
    graph.addSPEdge(core, w1);
    graph.addSPEdge(core, w2);
    for (auto vs : {std::make_tuple(&u, 0, 2), std::make_tuple(&core, 1, 3)}) {
        ag::Vertex &v = *std::get<0>(vs);
        v.coverage_info.addResolutionSample(std::get<1>(vs), std::get<2>(vs), 0, 0);
        v.rc().coverage_info.addResolutionSample(v.size() - std::get<2>(vs), v.size() - std::get<1>(vs), 0, 0);
    }
    ASSERT_TRUE(graph.checkConsistency());

    graph.flipVertexLockFree(core.incFront());
    graph.removeMarked();

    std::vector<std::pair<size_t, size_t>> samples;
    for (ag::CoverageSamples::SampleView s : u.coverage_info.getSamples())
        samples.emplace_back(s.start(), s.finish());
    EXPECT_EQ(samples, (std::vector<std::pair<size_t, size_t>>{{0, 2}, {3, 5}}));
    EXPECT_TRUE(graph.checkConsistency());
}

namespace {
//    The graph (with its listeners' checks), unchanged read sequences, and the carried tracker against a fresh one.
    void CheckReadsAndTracker(ag::AssemblyGraph &graph, ag::AlignedReadStorage &storage, ag::SuffixTracker &tracker,
                              const std::map<std::string, std::string> &before, const std::string &ctx) {
        std::stringstream ss;
        std::string err;
        EXPECT_TRUE(graph.checkConsistency(ss)) << ctx << "\n" << ss.str();
        EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << ctx << "\n" << err;
        ag::SuffixTracker fresh(storage, graph, 0, 1000);
        fresh.fillFromStorage(QuietLogger(), 1);
        EXPECT_TRUE(SuffixContentEqual(tracker, fresh, graph, err)) << ctx << "\n" << err;
    }
}

// Read storage through a real SPG pipeline: conversion, a resolution by the reads' pairs, then MergeAllSPG (mergePath
// and mergePrefixPath). Read sequences never change and the storage stays consistent after every step.
TEST(SpgReadStorage, ConsistentThroughResolutionAndMerging) {
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
//        Declared before the listeners, so they are destroyed first.
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, 40, 60, 1, seed + 100));
        std::map<std::string, std::string> before = ReadSeqSnapshot(storage);
//        A tracker can not be carried through the conversion, so the DBG phase has its own.
        {
            ag::SuffixTracker dbg_tracker(storage, rdbg.dbg, 0, 1000);
            dbg_tracker.fillFromStorage(QuietLogger(), 1);
            ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
            CheckReadsAndTracker(rdbg.dbg, storage, dbg_tracker, before,
                                 "seed=" + std::to_string(seed) + " after DBG merging");
        }
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        ag::SuffixTracker tracker(storage, spg, 0, 1000);
        tracker.fillFromStorage(QuietLogger(), 1);
        auto check = [&](ag::AssemblyGraph &graph, const std::string &step) {
            CheckReadsAndTracker(graph, storage, tracker, before, "seed=" + std::to_string(seed) + " " + step);
        };
        check(spg, "after conversion");
        if (!ResolveByReads(spg, storage, 5))
            continue;
        check(spg, "after resolution");
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        check(spg, "after merging");
    }
}

// Invariant: once the suffix tracker is attached (after conversion), graph operations never shrink read codes at the end
// of either strand; only whole leading edges may disappear. Checked after every resolution, merging and flipping step of
// the supregraph pipeline.
TEST(SpgReadStorage, ReadCodesOnlyExtendedAfterConversion) {
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, 40, 60, 1, seed + 100));
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        ag::SuffixTracker tracker(storage, spg, 0, 1000);
        tracker.fillFromStorage(QuietLogger(), 1);
        auto codes = ReadCodesSnapshot(storage);
        auto check = [&](const std::string &step) {
            std::string err;
            EXPECT_TRUE(ReadCodesOnlyExtended(codes, storage, err)) << "seed=" << seed << " " << step << "\n" << err;
            codes = ReadCodesSnapshot(storage);
        };
        for (size_t i = 0; i < 10 && ResolveByReads(spg, storage, 5); i++)
            check("after resolution " + std::to_string(i));
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        check("after merging paths");
        FlipAll(spg);
        spg.removeMarked();
        check("after flips");
    }
}

// Reads and the suffix tracker through flips: conversion and a resolution by the reads' pairs, then every vertex
// with a single incoming suffix edge and out-degree other than 1 is flipped (FromDBG already flips the initial ones).
TEST(SpgReadStorage, ConsistentThroughFlips) {
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, 40, 60, 1, seed + 100));
        std::map<std::string, std::string> before = ReadSeqSnapshot(storage);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        ag::SuffixTracker tracker(storage, spg, 0, 1000);
        tracker.fillFromStorage(QuietLogger(), 1);
        ResolveByReads(spg, storage, 5);
        std::vector<ag::VertexId> cores;
        for (ag::Vertex &v : spg.vertices())
            if (v.inDeg() == 1 && v.incFront().isSuffix() && v.outDeg() != 1)
                cores.push_back(v.getId());
        size_t flipped = 0;
        for (ag::VertexId vid : cores) {
            if (vid->marked() || vid->inDeg() != 1 || !vid->incFront().isSuffix() || vid->outDeg() == 1)
                continue;
            spg.flipVertexLockFree(vid->incFront());
            flipped++;
        }
        spg.removeMarked();
        CheckReadsAndTracker(spg, storage, tracker, before,
                             "seed=" + std::to_string(seed) + " flipped=" + std::to_string(flipped));
    }
}

// K + A + rc(K) (A not a palindrome) gives two different edges e and e.rc() that both go from K to rc(K), so both
// directions of a read lying on e sit in K's bucket. Converting e must leave the other direction to the call for e.rc()
// (regression: it was rewritten together with the first one and crashed).
TEST(SpgReadStorage, ConversionOfEdgeFromVertexToItsRc) {
    const size_t k = 11;
    std::string kmer = CollapsedGenome(k, 3);
    std::string a = CollapsedGenome(40, 4);
//    rc(A) starts with the complement of A's last nucleotide: make it differ from A's first, so the two strands leave K
//    on different edges.
    for (char c : std::string("ACGT"))
        if (c != Sequence(std::string(1, a[0])).rc().str()[0] && c != a[a.size() - 2])
            a.back() = c;
    std::string s = kmer + a + Sequence(kmer).rc().str();
    SeqDBG g({s}, k);
    std::vector<ag::AlignedRead> reads;
    for (size_t len : std::vector<size_t>{13, 16, 25, s.size()})
        for (size_t start = 0; start + len <= s.size(); start += 2)
            reads.emplace_back("r" + std::to_string(len) + "_" + std::to_string(start), g.align(s.substr(start, len)));
    ag::AlignedReadStorage storage(g.dbg, std::move(reads));
    std::map<std::string, std::string> before = ReadSeqSnapshot(storage);
    ag::MergeAllToEdges(QuietLogger(), 1, g.dbg);
    ag::Vertex *u = nullptr;
    for (ag::Vertex &v : g.dbg.vertices())
        if (v.getSeq().str() == kmer)
            u = &v;
    ASSERT_NE(u, nullptr);
    size_t to_rc = 0;
    for (ag::Edge &e : *u)
        if (e.getFinish() == u->rc() && e != e.rc())
            to_rc++;
    ASSERT_EQ(to_rc, 2u) << "both e and e.rc() must go from K to rc(K)";
    ConvertToSupregraph(g.dbg);
    std::stringstream ss;
    EXPECT_TRUE(storage.checkConsistency(ss)) << ss.str();
    std::string err;
    EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << err;
}
