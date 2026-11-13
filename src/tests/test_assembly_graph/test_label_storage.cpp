// multigraph::LabelStorage: every edge's label lists the original edges it was merged from, in order.

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "dbg/multi_graph.hpp"

using namespace ag_test;

namespace {
    std::string IdStr(const ag::EdgeId &id) {
        std::stringstream ss;
        ss << id;
        return ss.str();
    }
}

TEST(LabelStorage, LabelsSpellEdgesThroughMergeAllToEdges) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(1000, seed), 5);
        multigraph::LabelStorage labels(rdbg.dbg);
        std::map<std::string, std::string> original;
        for (ag::Edge &e : rdbg.dbg.edges())
            original[IdStr(e.getId())] = e.truncSeq().str();
        ASSERT_TRUE(rdbg.dbg.checkConsistency()) << "seed=" << seed;
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        EXPECT_TRUE(rdbg.dbg.checkConsistency()) << "seed=" << seed;
        size_t merged = 0;
        for (ag::Edge &e : rdbg.dbg.edges()) {
            const std::vector<ag::EdgeId> &label = labels.getLabel(e);
            std::string spelled;
            for (const ag::EdgeId &id : label) {
                auto it = original.find(IdStr(id));
                ASSERT_NE(it, original.end()) << "label of " << e.getId() << " names a non-original edge " << id;
                spelled += it->second;
            }
            EXPECT_EQ(spelled, e.truncSeq().str()) << "label of edge " << e.getId() << " seed=" << seed;
            if (label.size() > 1)
                merged++;
        }
        EXPECT_GT(merged, 0u);
    }
}
