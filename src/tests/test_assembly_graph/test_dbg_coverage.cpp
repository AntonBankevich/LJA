// dbg::CoverageTracker: edge coverage must always equal what the stored reads give (CoverageTracker::checkConsistency,
// run through the graph's checkConsistency) through every DBG operation it supports.

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"

using namespace ag_test;

namespace {
    std::string RC(const std::string &s) { return Sequence(s).rc().str(); }

//    Every substring of length len at the given step, and its reverse complement.
    std::vector<std::string> Tiles(const std::string &s, size_t len, size_t step) {
        std::vector<std::string> res;
        for (size_t i = 0; i + len <= s.size(); i += step) {
            res.push_back(s.substr(i, len));
            res.push_back(RC(s.substr(i, len)));
        }
        return res;
    }

    std::vector<ag::AlignedRead> AlignAll(SeqDBG &g, const std::vector<std::string> &seqs) {
        std::vector<ag::AlignedRead> reads;
        for (size_t i = 0; i < seqs.size(); i++)
            reads.emplace_back("r" + std::to_string(i), g.align(seqs[i]));
        return reads;
    }

    size_t TotalCoverage(ag::AssemblyGraph &graph) {
        size_t res = 0;
        for (ag::Edge &e : graph.edges())
            res += e.intCov();
        return res;
    }

    ag::Vertex *FindVertex(ag::AssemblyGraph &graph, const std::string &seq) {
        for (ag::Vertex &v : graph.vertices())
            if (v.getSeq().str() == seq)
                return &v;
        return nullptr;
    }
}

TEST(DBGCoverage, RandomThroughMergeAllToEdges) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(1000, seed), 5);
        std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 200, 40, 1, seed);
        dbg::DBGAlignedReadStorage storage(QuietLogger(), 1, rdbg.dbg, std::move(reads), true);
        ASSERT_GT(TotalCoverage(rdbg.dbg), 0u);
        ASSERT_TRUE(rdbg.dbg.checkConsistency()) << "seed=" << seed;
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        EXPECT_TRUE(rdbg.dbg.checkConsistency()) << "seed=" << seed;
    }
}

// S = X + rc(X) contains a self-rc edge (through its centre); every read on it contributes its forward and
// backward (k+1)-mers to the same edge, so it gets double coverage.
TEST(DBGCoverage, SelfRcEdgeGetsDoubleCoverage) {
    std::string x = CollapsedGenome(30, 5);
    std::string s = x + RC(x);
    SeqDBG g({s}, 7);
    dbg::DBGAlignedReadStorage storage(QuietLogger(), 1, g.dbg, AlignAll(g, Tiles(s, 15, 4)), true);
    ASSERT_TRUE(g.dbg.checkConsistency());
    ag::MergeAllToEdges(QuietLogger(), 1, g.dbg);
    size_t self_rc = 0;
    for (ag::Edge &e : g.dbg.edges())
        if (e == e.rc() && e.intCov() > 0)
            self_rc++;
    ASSERT_GT(self_rc, 0u) << "the palindrome must leave a covered self-rc edge";
    EXPECT_TRUE(g.dbg.checkConsistency());
}

// s1 and s2 share J, so the edge from s1's start ends at a vertex with two incoming edges. That edge is split
// in the middle, with reads inside each piece, ending just before the split vertex, and starting at the split
// vertex and continuing through J (the case where the read leaves the split edge).
TEST(DBGCoverage, SplitEdge) {
    const size_t k = 11;
    std::string x = CollapsedGenome(30, 1), j = CollapsedGenome(14, 2), y = CollapsedGenome(20, 3);
    std::string z = CollapsedGenome(20, 4), w = CollapsedGenome(20, 6);
    std::string s1 = x + j + y, s2 = z + j + w;
    SeqDBG g({s1, s2}, k);
//    The split vertex is s1[15, 26) and the first piece spells s1[0, 26), the second s1[15, 41).
    std::vector<std::string> seqs = {s1.substr(0, 20), s1.substr(2, 14), s1.substr(15, 35), s1.substr(16, 14)};
    for (size_t i = 0, n = seqs.size(); i < n; i++)
        seqs.push_back(RC(seqs[i]));
    dbg::DBGAlignedReadStorage storage(QuietLogger(), 1, g.dbg, AlignAll(g, seqs), true);
    ag::MergeAllToEdges(QuietLogger(), 1, g.dbg);
    ASSERT_TRUE(g.dbg.checkConsistency());
    ag::Vertex *start = FindVertex(g.dbg, s1.substr(0, k));
    ASSERT_NE(start, nullptr);
    ASSERT_EQ(start->outDeg(), 1u);
    ag::Edge &e = start->front();
    ASSERT_EQ(e.getFinish().inDeg(), 2u);
    ASSERT_GE(e.truncSize(), 29u) << "the junction is at J (possibly a nucleotide earlier)";
    g.dbg.splitEdge(e, {ag::EdgePosition(e, 15)});
    EXPECT_TRUE(g.dbg.checkConsistency());
}

// Two tips overlapping by 10 nucleotides with one mismatch (so they share no k-mer) are closed into one edge.
TEST(DBGCoverage, MergeTipsToEdge) {
    const size_t k = 7;
    std::string genome = CollapsedGenome(80, 7);
    std::string left = genome.substr(0, 45), right = genome.substr(35);
    char &m = left[40];
    for (char c : std::string("ACGT"))
        if (c != genome[39] && c != genome[40] && c != genome[41] && c != genome[42]) { m = c; break; }
    if (m == genome[40])
        for (char c : std::string("ACGT"))
            if (c != genome[39] && c != genome[40] && c != genome[41]) { m = c; break; }
    ASSERT_NE(m, genome[40]);
    SeqDBG g({left, right}, k);
    std::vector<std::string> seqs = Tiles(left, 15, 3);
    for (const std::string &s : Tiles(right, 15, 3))
        seqs.push_back(s);
    seqs.push_back(left.substr(left.size() - 12));
    seqs.push_back(RC(right.substr(0, 12)));
    dbg::DBGAlignedReadStorage storage(QuietLogger(), 1, g.dbg, AlignAll(g, seqs), true);
    ag::MergeAllToEdges(QuietLogger(), 1, g.dbg);
    ASSERT_TRUE(g.dbg.checkConsistency());
    ag::Vertex *left_tip = FindVertex(g.dbg, left.substr(left.size() - k));
    ag::Vertex *right_tip = FindVertex(g.dbg, right.substr(0, k));
    ASSERT_NE(left_tip, nullptr);
    ASSERT_NE(right_tip, nullptr);
    ASSERT_EQ(left_tip->inDeg(), 1u);
    ASSERT_EQ(right_tip->outDeg(), 1u);
    g.dbg.mergeTipsToEdge(left_tip->incFront(), right_tip->front(), AlignmentForm::Equal(10));
    g.dbg.removeMarked();
    EXPECT_TRUE(g.dbg.checkConsistency());
    EXPECT_GT(TotalCoverage(g.dbg), 0u);
}
