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
