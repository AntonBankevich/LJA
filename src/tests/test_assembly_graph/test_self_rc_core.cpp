// A supregraph core that is its own reverse complement: genome X P Y with a palindrome P. Every occurrence of P is also
// an occurrence of its rc, so P's vertex is entered from X and from rc(Y) and left into Y and into rc(X).

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"
#include "supregraph/coverage_tracker.hpp"
#include "supregraph/supregraph.hpp"

using namespace ag_test;

namespace {
    bool SomeVertexSpells(spg::Supregraph &spg, const std::string &genome) {
        std::string rc = (!Sequence(genome)).str();
        for (ag::Vertex &v : spg.vertices()) {
            std::string seq = v.getSeq().str();
            if (seq.find(genome) != std::string::npos || seq.find(rc) != std::string::npos)
                return true;
        }
        return false;
    }
}

// Resolves the self-rc core with the genome's pair (the plan adds its rc pair itself), merges, and checks reads,
// coverage, suffixes, the genome's DLL path and the supregraph invariants.
TEST(SelfRcCore, ResolveAndMergeKeepEverythingConsistent) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        std::string half = RandomGenome(30, seed + 50);
        std::string palindrome = half + (!Sequence(half)).str();
        std::string x = RandomGenome(300, seed);
        std::string y = RandomGenome(300, seed + 100);
        std::string genome = x + palindrome + y;
        RandomDBG rdbg(genome, 21);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::DLLAlignmentStorage dll(rdbg.dbg);
        {
            dbg::KmerIndex index(rdbg.dbg);
            index.fillAnchors(QuietLogger(), 1, rdbg.dbg, 2);
            Contig contig(Sequence(genome), "genome");
            std::vector<ag::AlignmentChain<Contig, ag::Edge>> al = index.carefulAlign(contig);
            ASSERT_FALSE(al.empty()) << "seed=" << seed;
            dll.addContig(contig, al);
        }
        ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, 2 * genome.size() / 10, 100, 1, seed + 200));
        ag::SuffixTracker suffixes(storage, rdbg.dbg, 0, 1000);
        suffixes.fillFromStorage(QuietLogger(), 1);
        std::map<std::string, std::string> before = ReadSeqSnapshot(storage);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        spg::CoverageTracker cov(spg, storage, QuietLogger(), 1);
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        ASSERT_TRUE(spg.checkConsistency(std::cerr)) << "seed=" << seed << " before resolution";

        ag::Vertex *core = nullptr;
        for (ag::Vertex &v : spg.vertices())
            if (v == v.rc() && v.inDeg() == 2 && v.outDeg() == 2)
                core = &v;
        ASSERT_NE(core, nullptr) << "seed=" << seed << ": no self-rc fork";
        EXPECT_NE(core->getSeq().str().find(palindrome), std::string::npos) << "seed=" << seed;
//        The genome passes the core from its X side into its Y side; the other pair is the rc of this one.
        ag::Edge *from_x = nullptr;
        ag::Edge *to_y = nullptr;
        for (ag::Edge &e : core->incoming())
            if (genome.find(e.getStart().getSeq().str()) != std::string::npos)
                from_x = &e;
        for (ag::Edge &e : *core)
            if (genome.find(e.getFinish().getSeq().str()) != std::string::npos)
                to_y = &e;
        ASSERT_TRUE(from_x != nullptr && to_y != nullptr) << "seed=" << seed;
        ag::VertexResolutionPlan plan(*core);
        plan.add(*from_x, *to_y);
        ASSERT_TRUE(plan.allConnected()) << "seed=" << seed << ": the rc pair was not added";
        spg.resolveVertex(*core, plan);
        spg.removeMarked();
        spg::MergeAllSPG(QuietLogger(), 1, spg);

        for (ag::Vertex &v : spg.vertices())
            EXPECT_FALSE(v == v.rc() && v.inDeg() == 2 && v.outDeg() == 2) << "seed=" << seed << ": fork left";
        EXPECT_TRUE(SomeVertexSpells(spg, genome)) << "seed=" << seed;
        std::string err;
        EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << "seed=" << seed << "\n" << err;
        EXPECT_TRUE(storage.checkConsistency(std::cerr)) << "seed=" << seed;
        EXPECT_TRUE(cov.checkConsistency(std::cerr)) << "seed=" << seed;
        EXPECT_TRUE(dll.checkConsistency(std::cerr)) << "seed=" << seed;
        EXPECT_TRUE(spg.checkConsistency(std::cerr)) << "seed=" << seed;
    }
}
