#include "gtest/gtest.h"
#include "graph_test_utils.hpp"

using namespace ag_test;

namespace {
    // Builds a hub with `branchCount` independent chains of `branchLength` edges each, and one read
    // per branch ending at a different offset along that branch (so each branch contributes a
    // distinct, independently-mergeable unbranching path plus its own read statistics).
    struct StarFixture {
        GraphBuilder gb;
        std::vector<std::vector<ag::EdgeId>> branches;
        ag::AlignedReadStorage storage;
        ag::SuffixTracker tracker;

        static std::vector<ag::AlignedRead> makeReads(const std::vector<std::vector<ag::EdgeId>> &branches) {
            std::vector<ag::AlignedRead> reads;
            for (size_t b = 0; b < branches.size(); b++) {
                // End each read at a varying offset (between 1 and branchLength edges).
                size_t endAt = 1 + (b % branches[b].size());
                reads.emplace_back("read" + std::to_string(b), MakePath(branches[b], 0, endAt));
            }
            return reads;
        }

        StarFixture(size_t branchCount, size_t branchLength) : gb(8), branches(gb.buildStar(branchCount, branchLength)),
                    storage(gb.graph, makeReads(branches)),
                    tracker(storage, gb.graph, /*min_len=*/0, /*max_len=*/1000) {
            tracker.fillFromStorage(TestLogger(), 1);
        }
    };
}

TEST(GraphModificationTest, SingleThreadedMergeAllToEdgesStaysConsistent) {
    StarFixture f(20, 6);
    ASSERT_TRUE(f.storage.checkConsistency());

    ag::MergeAllToEdges(TestLogger(), 1, f.gb.graph);

    EXPECT_TRUE(f.storage.checkConsistency());
}

// Regression test for the SuffixTracker data race fixed alongside this test: MergePathsToEdges runs
// with real thread parallelism, and before the fix, SuffixTracker::fireMergePathToEdge (and the other
// merge/split/resolve listeners) read edge_data without any lock while fireAddEdge/fireDeleteEdge
// mutated it under a different, unrelated lock (storage's), racing across concurrently-merged
// branches. Rebuilds a fresh graph each iteration since MergeAllToEdges consumes it.
//
// Note: this was checked to still reliably pass (i.e. not reproduce the crash) even when run
// against the old, unlocked SuffixTracker code, at this scale and well beyond (hundreds of
// branches, dozens of iterations). Data races on a small unordered_map critical section are
// timing-sensitive and evidently don't manifest reliably in a fast, small-scale unit test even
// though the race is real - this test is kept as a real exercise of the concurrent code path
// and a sanity net, not as a guaranteed reproducer of that specific bug.
TEST(GraphModificationTest, ParallelMergeAllToEdgesStaysConsistent) {
    for (int iter = 0; iter < 10; iter++) {
        StarFixture f(100, 6);
        ASSERT_TRUE(f.storage.checkConsistency());

        ag::MergeAllToEdges(TestLogger(), 8, f.gb.graph);

        ASSERT_TRUE(f.storage.checkConsistency());
    }
}

// flipVertex: core has a unique incoming suffix edge inc: u->core (core is a suffix of u) and prefix
// out-edges core->w1, core->w2. After the flip u connects directly to w1, w2 by edges spelling
// inc->e_i, core disappears, and every new edge's code is the concatenation of the replaced path's
// codes, in both orientations.
TEST(GraphModificationTest, FlipVertexEdgesAndCodes) {
    ag::AssemblyGraph graph;
    ag::Vertex &u = graph.addVertex(Sequence("TTAACCG"));
    ag::Vertex &core = graph.addVertex(Sequence("AACCG"));
    ag::Vertex &w1 = graph.addVertex(Sequence("AACCGA"));
    ag::Vertex &w2 = graph.addVertex(Sequence("AACCGT"));
    ag::Edge &inc = graph.addSPEdge(u, core);
    ag::Edge &e1 = graph.addSPEdge(core, w1);
    ag::Edge &e2 = graph.addSPEdge(core, w2);
    ASSERT_TRUE(inc.isSuffix());
    ASSERT_TRUE(e1.isPrefix());
    ASSERT_TRUE(e2.isPrefix());
    std::map<int, std::pair<std::string, std::string>> expected; // finish -> (code, rc code)
    for (ag::Edge *e : {&e1, &e2})
        expected[e->getFinish().getInnerId()] = {(inc.getCode() + e->getCode()).str(),
                                            (e->rc().getCode() + inc.rc().getCode()).str()};
    Sequence s1 = u.getSeq() + inc.truncSeq() + e1.truncSeq();
    Sequence s2 = u.getSeq() + inc.truncSeq() + e2.truncSeq();

    ag::Vertex &res = graph.flipVertexLockFree(core.incFront());
    graph.removeMarked();

    ASSERT_EQ(res, u);
    ASSERT_EQ(u.outDeg(), 2u);
    for (ag::Edge &e : u) {
        ASSERT_TRUE(e.getFinish() == w1 || e.getFinish() == w2);
        EXPECT_EQ(e.fullSeq(), e.getFinish() == w1 ? s1 : s2);
        EXPECT_EQ(e.rc().fullSeq(), !e.fullSeq());
        EXPECT_EQ(e.getCode().str(), expected[e.getFinish().getInnerId()].first);
        EXPECT_EQ(e.rc().getCode().str(), expected[e.getFinish().getInnerId()].second);
        EXPECT_EQ(e.getCode()[0], e.truncSeq()[0]);
    }
    EXPECT_EQ(w1.inDeg(), 1u);
    EXPECT_EQ(w2.inDeg(), 1u);
    EXPECT_EQ(graph.size(), 6u); // core and core.rc() removed
    EXPECT_TRUE(graph.checkConsistency());
}

namespace {
    // v0 -p-> v1 -p-> peak -s-> s1 -s-> s2 with sequences AAC, AACG, AACGTA, CGTA, GTA.
    struct PeakChain {
        ag::AssemblyGraph graph;
        std::vector<ag::VertexId> vs;
        std::vector<ag::EdgeId> es;

        PeakChain() {
            for (const char *s : {"AAC", "AACG", "AACGTA", "CGTA", "GTA"})
                vs.push_back(graph.addVertex(Sequence(s)).getId());
            for (size_t i = 0; i + 1 < vs.size(); i++)
                es.push_back(graph.addSPEdge(*vs[i], *vs[i + 1]).getId());
        }
    };
}

// mergePath on a path of prefix edges followed by suffix edges: the peak already spells the path, so no vertex is
// created; every other vertex of the path is removed and the path's outside neighbours connect to the peak. Covers
// prefix+suffix, all-prefix, all-suffix and the one-prefix-one-suffix path.
TEST(GraphModificationTest, MergePathIntoExistingPeak) {
    for (std::pair<size_t, size_t> c : std::vector<std::pair<size_t, size_t>>{{0, 4}, {0, 2}, {2, 4}, {1, 3}}) {
        PeakChain f;
        ASSERT_TRUE(f.graph.checkConsistency());
        ag::GraphPath path = MakePath(f.es, c.first, c.second);

        ag::Vertex &res = f.graph.mergePath(path);
        f.graph.removeMarked();

        const std::string name = "case " + std::to_string(c.first) + "-" + std::to_string(c.second);
        EXPECT_EQ(res, *f.vs[2]) << name;
        EXPECT_EQ(f.graph.size(), 2 * (5 - (c.second - c.first))) << name;
        EXPECT_EQ(res.inDeg(), c.first > 0 ? 1u : 0u) << name;
        if (c.first > 0) {
            EXPECT_EQ(res.incFront().getStart(), *f.vs[c.first - 1]) << name;
            EXPECT_TRUE(res.incFront().isPrefix()) << name;
        }
        EXPECT_EQ(res.outDeg(), c.second < 4 ? 1u : 0u) << name;
        if (c.second < 4) {
            EXPECT_EQ(res.front().getFinish(), *f.vs[c.second + 1]) << name;
            EXPECT_TRUE(res.front().isSuffix()) << name;
        }
        EXPECT_TRUE(f.graph.checkConsistency()) << name;
    }
}
