// DLLAlignmentStorage: contig paths must stay consistent (DLLAlignmentStorage::checkConsistency) through every
// graph operation, on both DBG and supregraph.

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "supregraph/supregraph.hpp"
#include "supregraph/coverage_tracker.hpp"

using namespace ag_test;

namespace {
//    The genome's own path, registered while the graph is still a hashed DBG (as PrepareDLLPathTracker does).
    void AddGenome(ag::DLLAlignmentStorage &dll, RandomDBG &rdbg) {
        dbg::KmerIndex index(rdbg.dbg);
        index.fillAnchors(QuietLogger(), 1, rdbg.dbg, 2);
        Contig genome(Sequence(rdbg.genome), "genome");
        std::vector<ag::AlignmentChain<Contig, ag::Edge>> al = index.carefulAlign(genome);
        ASSERT_FALSE(al.empty());
        dll.addContig(genome, al);
    }
}

TEST(DLLPathStorage, SplitEdge) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(400, seed), 5);
        ag::DLLAlignmentStorage dll(rdbg.dbg);
        AddGenome(dll, rdbg);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed;
        ag::Edge *longest = nullptr;
        for (ag::Edge &e : rdbg.dbg.edges())
            if (longest == nullptr || e.truncSize() > longest->truncSize())
                longest = &e;
        ASSERT_GE(longest->truncSize(), 4u);
        rdbg.dbg.splitEdge(*longest, {ag::EdgePosition(*longest, longest->truncSize() / 2)});
        EXPECT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed;
    }
}

// Two tips overlapping by 10 nucleotides with one mismatch are closed into one edge; the contigs spelling the tips
// keep only the parts the merge did not modify.
TEST(DLLPathStorage, MergeTipsToEdge) {
    const size_t k = 7;
    std::string genome = CollapsedGenome(80, 7);
    std::string left = genome.substr(0, 45), right = genome.substr(35);
    char &m = left[40];
    for (char c : std::string("ACGT"))
        if (c != genome[39] && c != genome[40] && c != genome[41] && c != genome[42]) { m = c; break; }
    if (m == genome[40])
        for (char c : std::string("ACGT"))
            if (c != genome[39] && c != genome[40] && c != genome[41]) { m = c; break; }
    SeqDBG g({left, right}, k);
    ag::DLLAlignmentStorage dll(g.dbg);
    {
        dbg::KmerIndex index(g.dbg);
        index.fillAnchors(QuietLogger(), 1, g.dbg, 1);
        for (const std::pair<std::string, std::string> &c : {std::make_pair(std::string("left"), left), std::make_pair(std::string("right"), right)}) {
            Contig contig(Sequence(c.second), c.first);
            std::vector<ag::AlignmentChain<Contig, ag::Edge>> al = index.carefulAlign(contig);
            dll.addContig(contig, al);
        }
    }
    ag::MergeAllToEdges(QuietLogger(), 1, g.dbg);
    ASSERT_TRUE(dll.checkConsistency(std::cerr));
    ag::Vertex *left_tip = nullptr, *right_tip = nullptr;
    for (ag::Vertex &v : g.dbg.vertices()) {
        if (v.getSeq().str() == left.substr(left.size() - k)) left_tip = &v;
        if (v.getSeq().str() == right.substr(0, k)) right_tip = &v;
    }
    ASSERT_NE(left_tip, nullptr);
    ASSERT_NE(right_tip, nullptr);
    g.dbg.mergeTipsToEdge(left_tip->incFront(), right_tip->front(), AlignmentForm::Equal(10));
    g.dbg.removeMarked();
    EXPECT_TRUE(dll.checkConsistency(std::cerr));
    EXPECT_FALSE(dll.getPath("left").empty());
    EXPECT_FALSE(dll.getPath("right").empty());
}

namespace {
//    Gaps (consecutive fragments that are not linked) and the number of contig positions the path covers.
    std::pair<size_t, size_t> PathShape(ag::DLLAlignmentPath &path) {
        size_t gaps = 0, covered = 0, reach = 0;
        const ag::AlignmentFragment *prev = nullptr;
        for (const ag::AlignmentFragment &f : path) {
            if (prev != nullptr && !prev->isLinked(f))
                gaps++;
            size_t from = std::max<size_t>(reach, f.seg.left);
            if (f.seg.right > from)
                covered += f.seg.right - from;
            reach = std::max<size_t>(reach, f.seg.right);
            prev = &f;
        }
        return {gaps, covered};
    }

//    Resolves every 2-in/2-out k-mer junction the genome passes by exactly the pairs the genome uses there.
    size_t ResolveByGenome(ag::AssemblyGraph &graph, ag::DLLAlignmentPath &path, size_t k) {
        std::vector<std::pair<ag::EdgeId, ag::EdgeId>> used;
        const ag::AlignmentFragment *prev = nullptr;
        for (const ag::AlignmentFragment &f : path) {
            if (prev != nullptr && prev->isLinked(f) && prev->edge->getFinish().size() == k &&
                    prev->edge->getFinish().inDeg() == 2 && prev->edge->getFinish().outDeg() == 2)
                used.emplace_back(prev->edge, f.edge);
            prev = &f;
        }
//        Resolving a core resolves its rc too, so the pairs used at core.rc() go into core's plan as their rc.
        std::unordered_map<ag::VertexId, std::vector<std::pair<ag::EdgeId, ag::EdgeId>>> by_core;
        for (const auto &p : used) {
            ag::Vertex &core = p.first->getFinish();
            if (core.isCanonical())
                by_core[core.getId()].push_back(p);
            else
                by_core[core.rc().getId()].emplace_back(p.second->rc().getId(), p.first->rc().getId());
        }
        size_t resolved = 0;
        for (auto &c : by_core)
            if (!c.first->marked()) {
                ag::VertexResolutionPlan plan(*c.first);
                for (const auto &p : c.second)
                    plan.add(*p.first, *p.second);
                graph.resolveVertex(*c.first, plan);
                resolved++;
            }
        graph.removeMarked();
        return resolved;
    }
}

// Resolving by the genome's own pairs, merging and flipping must not lose any part of the genome's path.
TEST(DLLPathStorage, SupregraphPipelineLossless) {
    size_t resolved = 0;
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::DLLAlignmentStorage dll(rdbg.dbg);
        AddGenome(dll, rdbg);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        std::pair<size_t, size_t> before = PathShape(dll.getPath("genome"));
        resolved += ResolveByGenome(spg, dll.getPath("genome"), 5);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after resolution";
        EXPECT_EQ(PathShape(dll.getPath("genome")), before) << "seed=" << seed << " after resolution";
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after MergeAllSPG";
        EXPECT_EQ(PathShape(dll.getPath("genome")), before) << "seed=" << seed << " after MergeAllSPG";
        FlipAll(spg);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after flips";
        EXPECT_EQ(PathShape(dll.getPath("genome")), before) << "seed=" << seed << " after flips";
    }
    EXPECT_GT(resolved, 0u);
}

// The genome's path through the whole supregraph pipeline: conversion, resolution, MergeAllSPG, flips.
TEST(DLLPathStorage, SupregraphPipeline) {
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
//        The supregraph must outlive the storage attached to it.
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::DLLAlignmentStorage dll(rdbg.dbg);
        AddGenome(dll, rdbg);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after MergeAllToEdges";
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after conversion";
        ResolveMatching(spg, 5);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after resolution";
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after MergeAllSPG";
        FlipAll(spg);
        ASSERT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed << " after flips";
    }
}

// SimpleRemoveUncovered isolates unsupported outer vertices in parallel (locked isolateAndMark) and then merges with
// several threads. Reads cover only the first half of the genome, so much of the graph goes; the genome's path, the
// reads and the coverage all stay consistent.
TEST(DLLPathStorage, ParallelSimpleRemoveUncovered) {
    size_t removed = 0;
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(600, seed), 5);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::DLLAlignmentStorage dll(rdbg.dbg);
        AddGenome(dll, rdbg);
        RandomDBG half(rdbg.genome.substr(0, 300), 5);
        std::vector<ag::AlignedRead> reads;
        {
            dbg::KmerIndex index(rdbg.dbg);
            index.fillAnchors(QuietLogger(), 1, rdbg.dbg, 1);
            for (size_t start = 0; start + 40 <= 300; start += 5) {
                std::string name = "r" + std::to_string(start);
                reads.emplace_back(name, index.align(Sequence(rdbg.genome.substr(start, 40)), name));
            }
        }
        ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        spg::CoverageTracker cov(spg, storage, QuietLogger(), 1);
        size_t before = spg.size();
        removed += spg::SimpleRemoveUncovered(QuietLogger(), 4, spg);
        spg.removeMarked();
        EXPECT_LT(spg.size(), before) << "seed=" << seed;
        std::stringstream ss;
        EXPECT_TRUE(spg.checkConsistency(ss)) << "seed=" << seed << "\n" << ss.str();
        EXPECT_TRUE(storage.checkConsistency(ss)) << "seed=" << seed << "\n" << ss.str();
    }
    EXPECT_GT(removed, 0u);
}

namespace {
//    Positions of all (possibly overlapping) occurrences of pattern in text.
    std::vector<size_t> FindAll(const std::string &text, const std::string &pattern) {
        std::vector<size_t> res;
        for (size_t pos = text.find(pattern); pos != std::string::npos; pos = text.find(pattern, pos + 1))
            res.push_back(pos);
        return res;
    }

//    When the graph spells the genome exactly, every vertex with sample segments is traversed once per occurrence of
//    its label in the genome or its rc (besides partial traversals where the genome starts or ends inside it), and
//    every traversal away from the path ends spells the whole vertex. Vertices occurring at a genome end are skipped:
//    only proper-infix occurrences are guaranteed to stay in the path.
    void CheckTraversals(ag::AssemblyGraph &graph, ag::DLLAlignmentStorage &dll, const std::string &genome,
                         const std::string &stage) {
        const std::string rc_genome = Sequence(genome).rc().str();
        auto traversals = dll.vertexTraversals();
        size_t checked = 0;
        for (ag::Vertex &v : graph.vertices()) {
            if (v.coverage_info.empty())
                continue;
            const std::string label = v.getSeq().str();
            std::vector<size_t> occurrences = FindAll(genome, label);
            for (size_t pos : FindAll(rc_genome, label))
                occurrences.push_back(pos);
            bool at_genome_end = std::any_of(occurrences.begin(), occurrences.end(), [&](size_t pos) {
                return pos == 0 || pos + label.size() == genome.size();
            });
            if (at_genome_end)
                continue;
            checked++;
            auto it = traversals.find(v.getId());
//            A traversal cut by a path end may be a partial copy of the vertex at a genome end, which FindAll does not see.
            size_t found = it == traversals.end() ? 0 : std::count_if(it->second.begin(), it->second.end(),
                    [&](const ag::DLLAlignmentStorage::VertexTraversal &t) {return !t.at_end || t.covered == label.size();});
            std::stringstream details;
            if (it != traversals.end())
                for (const ag::DLLAlignmentStorage::VertexTraversal &t : it->second)
                    details << " (" << t.contig->getInnerId() << ", covered " << t.covered << ", at end " << t.at_end
                            << ")";
            EXPECT_EQ(found, occurrences.size()) << stage << " vertex " << v.getId() << " of length " << label.size()
                                                 << ":" << details.str();
            if (it == traversals.end())
                continue;
            for (const ag::DLLAlignmentStorage::VertexTraversal &t : it->second)
                if (!t.at_end)
                    EXPECT_EQ(t.covered, label.size()) << stage << " vertex " << v.getId();
        }
        EXPECT_GT(checked, 0u) << stage;
    }
}

// Traversals of the tracked genome path give the true multiplicity of every vertex: right after conversion to a
// supregraph and after resolution by the genome's pairs, merges and flips.
TEST(DLLPathStorage, TraversalsAreTrueMultiplicity) {
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), 5);
        const std::string genome = rdbg.genome;
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::DLLAlignmentStorage dll(rdbg.dbg);
        AddGenome(dll, rdbg);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        CheckTraversals(spg, dll, genome, "seed=" + std::to_string(seed) + " after conversion");
        ResolveByGenome(spg, dll.getPath("genome"), 5);
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        FlipAll(spg);
        CheckTraversals(spg, dll, genome, "seed=" + std::to_string(seed) + " after resolution, merges and flips");
    }
}
