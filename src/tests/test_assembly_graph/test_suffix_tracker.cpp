#include "gtest/gtest.h"
#include "graph_test_utils.hpp"

using namespace ag_test;

namespace {
    // Builds a chain of 6 vertices (5 edges), two reads:
    //  read0: edges 0,1,2 (ends inside the chain, at edge 2)
    //  read1: edges 0,1,2,3,4 (spans the whole chain)
    // Returns the storage and tracker together with the chain so tests can inspect both.
    struct Fixture {
        GraphBuilder gb;
        std::vector<ag::EdgeId> chain;
        ag::AlignedReadStorage storage;
        ag::SuffixTracker tracker;

        static std::vector<ag::AlignedRead> makeReads(const std::vector<ag::EdgeId> &chain) {
            std::vector<ag::AlignedRead> reads;
            reads.emplace_back("short", MakePath(chain, 0, 3));
            reads.emplace_back("long", MakePath(chain, 0, 5));
            return reads;
        }

        Fixture() : gb(8), chain(gb.buildChain(6)),
                    storage(gb.graph, makeReads(chain)),
                    tracker(storage, gb.graph, /*min_len=*/0, /*max_len=*/1000) {
            tracker.fillFromStorage(TestLogger(), 1);
        }
    };
}

TEST(SuffixTrackerTest, TracksReadEndsAndConsistency) {
    Fixture f;
    EXPECT_TRUE(f.storage.checkConsistency());

    // "short" truly ends at edge 2 (forward direction) -> its rc direction starts at edge2.rc(),
    // so edge2.rc()'s own SuffixRecord (tracking forward continuations after edge2.rc(), i.e. what
    // used to precede edge2) is not what we check here; instead getNumberOfEnds() on edge2 tracks
    // reads whose forward path terminates at edge2, matched against storage->startCnt(edge2.rc()).
    const ag::SuffixRecord &rec2 = f.tracker.getSuffixRecord(*f.chain[2]);
    EXPECT_EQ(rec2.getNumberOfEnds(), f.storage.startCnt(f.chain[2]->rc()));
    EXPECT_EQ(rec2.getNumberOfEnds(), 1u);

    // "long" passes through edge2 without ending there, and continues to edge4 where it truly ends.
    const ag::SuffixRecord &rec4 = f.tracker.getSuffixRecord(*f.chain[4]);
    EXPECT_EQ(rec4.getNumberOfEnds(), f.storage.startCnt(f.chain[4]->rc()));
    EXPECT_EQ(rec4.getNumberOfEnds(), 1u);

    // Edge0 has both reads passing through it (as the first edge of their paths), neither ends there.
    const ag::SuffixRecord &rec0 = f.tracker.getSuffixRecord(*f.chain[0]);
    EXPECT_EQ(rec0.getNumberOfEnds(), 0u);
    // Both reads continue past edge0 into edge1, so countStartsWith(empty) should count both paths.
    EXPECT_EQ(rec0.countStartsWith(ag::GraphPath(f.chain[0]->getFinish())), 2u);
}

TEST(SuffixTrackerTest, InvalidateReadUpdatesRecord) {
    Fixture f;
    ag::AlignedRead &shortRead = f.storage[0];
    ASSERT_EQ(shortRead.getId(), "short");

    f.storage.delayedInvalidateRead(shortRead, "test");
    f.storage.apply(shortRead);

    EXPECT_TRUE(f.storage.checkConsistency());
    const ag::SuffixRecord &rec2 = f.tracker.getSuffixRecord(*f.chain[2]);
    EXPECT_EQ(rec2.getNumberOfEnds(), 0u);
}

TEST(SuffixTrackerTest, MergePathToEdgeAggregatesEnds) {
    Fixture f;
    ag::Edge &merged = f.gb.graph.mergePathToEdge(MakePath(f.chain, 0, 5));

    EXPECT_TRUE(f.storage.checkConsistency());
    const ag::SuffixRecord &rec = f.tracker.getSuffixRecord(merged);
    // Both reads now end within (or at the far end of) the single merged edge: "short" ended
    // partway through what is now `merged`, "long" ends exactly at its far boundary.
    EXPECT_EQ(rec.getNumberOfEnds(), f.storage.startCnt(merged.rc()));
    EXPECT_EQ(rec.getNumberOfEnds(), 2u);
}

// splitEdge has a contract (not enforced by asserts, so it's easy to violate by accident when
// hand-picking split positions, as an earlier version of these tests did): splitting an edge e into
// a path P must never cut a read - i.e. no inner vertex of P may be an inner vertex of any read
// path that touches e. Equivalently, every read touching e must, after the split, either:
//   (a) end within P's first edge,
//   (b) start within P's last edge, or
//   (c) be entirely contained within a single edge of P.
// Each of the three tests below builds a split position that respects the contract for exactly one
// of these cases and checks SuffixTracker/AlignedReadStorage stay consistent across the split.

// Case (a): a read starting on an earlier edge and ending partway into the edge being split, at a
// position inside what becomes the split's first piece.
TEST(SuffixTrackerTest, SplitEdgeReadEndingInFirstPieceStaysConsistent) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(3); // edges: e0 (before), e1 (to be split)

    // Starts at e0, extends 3 bases into e1's truncSeq (out of 8), i.e. ends strictly inside what
    // will become e1's first piece [0,4).
    ag::GraphPath path(*chain[0]);
    path += *chain[1];
    path.setCutRight(5); // trims the last 5 of e1's 8 bases, leaving 3 covered

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r0", path);
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ag::SuffixTracker tracker(storage, gb.graph, 0, 1000);
    tracker.fillFromStorage(TestLogger(), 1);
    ASSERT_TRUE(storage.checkConsistency());
    ASSERT_EQ(tracker.getSuffixRecord(*chain[1]).getNumberOfEnds(), 1u);

    ag::GraphPath split = gb.graph.splitEdge(*chain[1], {ag::EdgePosition(*chain[1], 4)});

    EXPECT_TRUE(storage.checkConsistency());
    ag::Edge &firstPiece = split.frontEdge();
    EXPECT_EQ(tracker.getSuffixRecord(firstPiece).getNumberOfEnds(), storage.startCnt(firstPiece.rc()));
    EXPECT_EQ(tracker.getSuffixRecord(firstPiece).getNumberOfEnds(), 1u);
}

// Case (b): a read starting partway into the edge being split, at a position inside what becomes
// the split's last piece, and continuing on into a later edge.
TEST(SuffixTrackerTest, SplitEdgeReadStartingInLastPieceStaysConsistent) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(3); // edges: e0 (to be split), e1 (after)

    // Skips e0's start vertex (8 bases) plus the first 5 of e0's own 8 bases, i.e. starts strictly
    // inside what will become e0's second piece [4,8), then continues into e1.
    ag::GraphPath path(*chain[0], /*cut_left=*/13, /*cut_right=*/0);
    path += *chain[1];

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r0", path);
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ag::SuffixTracker tracker(storage, gb.graph, 0, 1000);
    tracker.fillFromStorage(TestLogger(), 1);
    ASSERT_TRUE(storage.checkConsistency());
    ASSERT_EQ(tracker.getSuffixRecord(*chain[1]).getNumberOfEnds(), 1u);

    gb.graph.splitEdge(*chain[0], {ag::EdgePosition(*chain[0], 4)});

    EXPECT_TRUE(storage.checkConsistency());
    // The read's end (at e1) is untouched by the split.
    EXPECT_EQ(tracker.getSuffixRecord(*chain[1]).getNumberOfEnds(), storage.startCnt(chain[1]->rc()));
    EXPECT_EQ(tracker.getSuffixRecord(*chain[1]).getNumberOfEnds(), 1u);
}

// Case (c): a read entirely contained within a single edge of the split, touching neither boundary.
TEST(SuffixTrackerTest, SplitEdgeReadContainedInOnePieceStaysConsistent) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(2); // single edge, truncSize == k == 8, fullSize == 16

    std::vector<ag::AlignedRead> reads;
    // Starts at edge0's very start and covers only its first 10 bases (fullSize 16, cut_right 6),
    // so after splitting at position 4 it's entirely contained within the first piece
    // (fullSize 4 + 8 == 12).
    reads.emplace_back("r0", ag::GraphPath(*chain[0], 0, 6));
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ag::SuffixTracker tracker(storage, gb.graph, 0, 1000);
    tracker.fillFromStorage(TestLogger(), 1);

    ASSERT_TRUE(storage.checkConsistency());
    ASSERT_EQ(tracker.getSuffixRecord(*chain[0]).getNumberOfEnds(), 1u);

    ag::GraphPath split = gb.graph.splitEdge(*chain[0], {ag::EdgePosition(*chain[0], 4)});

    EXPECT_TRUE(storage.checkConsistency());
    ag::Edge &firstPiece = split.frontEdge();
    EXPECT_EQ(tracker.getSuffixRecord(firstPiece).getNumberOfEnds(), storage.startCnt(firstPiece.rc()));
    EXPECT_EQ(tracker.getSuffixRecord(firstPiece).getNumberOfEnds(), 1u);
}

// A read entirely on an unrelated branch of the graph (touching the split edge not at all, so
// trivially compliant with splitEdge's contract) is untouched by splitting a distant edge.
TEST(SuffixTrackerTest, SplitEdgeDoesNotDisturbUnrelatedReads) {
    GraphBuilder gb(8);
    std::vector<std::vector<ag::EdgeId>> branches = gb.buildStar(2, 3); // hub + 2 branches of 3 edges each

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("onBranch1", MakePath(branches[1], 0, 3)); // fully on the untouched branch

    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ag::SuffixTracker tracker(storage, gb.graph, 0, 1000);
    tracker.fillFromStorage(TestLogger(), 1);

    ASSERT_TRUE(storage.checkConsistency());
    ASSERT_EQ(tracker.getSuffixRecord(*branches[1][2]).getNumberOfEnds(), 1u);

    gb.graph.splitEdge(*branches[0][1], {ag::EdgePosition(*branches[0][1], 4)});

    EXPECT_TRUE(storage.checkConsistency());
    EXPECT_EQ(tracker.getSuffixRecord(*branches[1][2]).getNumberOfEnds(), 1u);
}
