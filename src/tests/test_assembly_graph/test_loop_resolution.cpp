// Loop resolution in the supregraph: a tandem repeat X (R L)^n R Y gives a simple loop at R with body L. Reads are
// shorter than R, so no read crosses the core and the number of copies must come from coverage. k = 21 keeps the
// random parts free of repeats.

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"
#include "supregraph/coverage_tracker.hpp"
#include "supregraph/loop_resolution.hpp"
#include "supregraph/supregraph.hpp"

using namespace ag_test;

namespace {
    std::vector<spg::SimpleLoop> AllSimpleLoops(spg::Supregraph &spg) {
        std::vector<spg::SimpleLoop> res;
        for (ag::Vertex &v : spg.verticesUnique()) {
            std::optional<spg::SimpleLoop> loop = spg::DetectSimpleLoop(v);
            if (loop)
                res.push_back(*loop);
        }
        return res;
    }

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

TEST(LoopResolution, UnrollsTandemRepeatByCoverage) {
    for (size_t copies = 1; copies <= 4; copies++) {
        uint32_t seed = 10 + copies;
        std::string x = RandomGenome(300, seed);
        std::string r = RandomGenome(60, seed + 100);
        std::string l = RandomGenome(80, seed + 200);
        std::string y = RandomGenome(300, seed + 300);
        std::string genome = x;
        for (size_t i = 0; i < copies; i++)
            genome += r + l;
        genome += r + y;
        RandomDBG rdbg(genome, 21);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, genome.size(), 40, 1, seed + 400));
        ag::SuffixTracker suffixes(storage, rdbg.dbg, 0, 1000);
        suffixes.fillFromStorage(QuietLogger(), 1);
        std::map<std::string, std::string> before = ReadSeqSnapshot(storage);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        spg::CoverageTracker cov(spg, storage, QuietLogger(), 1);
        spg::MergeAllSPG(QuietLogger(), 1, spg);
        std::vector<spg::SimpleLoop> loops = AllSimpleLoops(spg);
        ASSERT_EQ(loops.size(), 1u) << "copies=" << copies;
        EXPECT_NEAR(spg::EstimateLoopCopies(loops[0]), double(copies), 0.4) << "copies=" << copies;

        EXPECT_EQ(spg::ResolveLoops(QuietLogger(), spg, suffixes, 10), 1u) << "copies=" << copies;
        spg::MergeAllSPG(QuietLogger(), 1, spg);

        EXPECT_TRUE(AllSimpleLoops(spg).empty()) << "copies=" << copies;
        EXPECT_TRUE(SomeVertexSpells(spg, genome)) << "copies=" << copies;
        std::string err;
        EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << "copies=" << copies << "\n" << err;
        EXPECT_TRUE(storage.checkConsistency(std::cerr)) << "copies=" << copies;
        EXPECT_TRUE(cov.checkConsistency(std::cerr)) << "copies=" << copies;
        EXPECT_TRUE(spg.checkConsistency(std::cerr)) << "copies=" << copies;
    }
}
