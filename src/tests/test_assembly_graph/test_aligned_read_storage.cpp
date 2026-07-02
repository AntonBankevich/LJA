#include "gtest/gtest.h"
#include "graph_test_utils.hpp"

using namespace ag_test;

TEST(AlignedReadStorageTest, BasicAddAndConsistency) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(6); // 5 edges

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("read0", MakePath(chain, 0, 3)); // edges 0,1,2
    reads.emplace_back("read1", MakePath(chain, 2, 5)); // edges 2,3,4

    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    EXPECT_TRUE(storage.checkConsistency());

    // read0 starts at edge 0 and ends at edge 2; read1 starts at edge 2 and ends at edge 4.
    EXPECT_EQ(storage.startCnt(*chain[0]), 1u);
    EXPECT_EQ(storage.startCnt(*chain[2]), 1u);
    // Both reads end within edges 2 and 4 respectively - i.e. rc() of those edges is where the
    // reverse-direction reads start.
    EXPECT_EQ(storage.startCnt(chain[2]->rc()), 1u);
    EXPECT_EQ(storage.startCnt(chain[4]->rc()), 1u);
}

TEST(AlignedReadStorageTest, InvalidateReadDropsStartCounts) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(4); // 3 edges

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("read0", MakePath(chain, 0, 3));

    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency());
    ASSERT_EQ(storage.startCnt(*chain[0]), 1u);

    ag::AlignedRead &read = storage[0];
    storage.delayedInvalidateRead(read, "test invalidate");
    storage.apply(read);

    EXPECT_FALSE(read.valid());
    EXPECT_EQ(storage.startCnt(*chain[0]), 0u);
    EXPECT_TRUE(storage.checkConsistency());
}
