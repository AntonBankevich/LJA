// IndependentPathBuckets: every unbranching path in exactly one bucket, no two paths of a bucket sharing a vertex up to
// rc, and the same split whatever the number of threads.

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "assembly_graph/ag_algorithms.hpp"
#include "supregraph/supregraph.hpp"
#include <set>

using namespace ag_test;

TEST(PathBuckets, BucketsAreIndependentAndThreadIndependent) {
    size_t multi_bucket_graphs = 0;
    for (uint32_t seed = 1; seed <= 10; seed++) {
        RandomDBG rdbg(RandomGenome(500, seed), 5);
        rdbg.dbg.disableHashing();
        spg::Supregraph spg = spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg));
        std::vector<ag::GraphPath> paths = ag::AllUnbranchingPaths(QuietLogger(), 1, spg);
        ASSERT_FALSE(paths.empty()) << "seed=" << seed;
        std::vector<std::vector<size_t>> buckets = ag::IndependentPathBuckets(QuietLogger(), 1, paths);
        EXPECT_EQ(ag::IndependentPathBuckets(QuietLogger(), 8, paths), buckets) << "seed=" << seed;
        std::vector<size_t> seen(paths.size(), 0);
        for (const std::vector<size_t> &bucket : buckets) {
            EXPECT_FALSE(bucket.empty()) << "seed=" << seed;
            std::set<ag::ConstVertexId> used;
            for (size_t i : bucket) {
                seen[i]++;
                std::set<ag::ConstVertexId> own;
                for (const ag::Vertex &v : paths[i].vertices())
                    own.insert(v.isCanonical() ? v.getId() : v.rc().getId());
                for (ag::ConstVertexId vid : own)
                    EXPECT_TRUE(used.insert(vid).second) << "seed=" << seed << " path " << i;
            }
        }
        for (size_t i = 0; i < paths.size(); i++)
            EXPECT_EQ(seen[i], 1u) << "seed=" << seed << " path " << i;
        multi_bucket_graphs += buckets.size() > 1;
    }
    EXPECT_GT(multi_bucket_graphs, 0u);
}
