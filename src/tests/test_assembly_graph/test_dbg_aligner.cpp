// dbg::KmerIndex::carefulAlign: every chain matches exactly, and a contig taken from the genome is covered by its
// chains from its very first k-mer, also when it starts inside an edge and is aligned through an anchor.

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"

using namespace ag_test;

namespace {
//    Checks every chain of contig's alignment and returns the contig k-mer start positions they cover.
    std::vector<bool> CheckChains(Contig &contig, const std::vector<ag::AlignmentChain<Contig, ag::Edge>> &als,
                                  size_t k, const std::string &ctx) {
        std::vector<bool> covered(contig.truncSize() + 1, false);
        std::string cseq = contig.getSeq().str();
        for (const ag::AlignmentChain<Contig, ag::Edge> &al : als) {
            const ag::Edge &edge = al.seg_to.contig();
            std::string eseq = edge.getSeq().str();
            EXPECT_EQ(al.seg_from.size(), al.seg_to.size()) << ctx;
            EXPECT_EQ(cseq.substr(al.seg_from.left, al.seg_from.size() + k),
                      eseq.substr(al.seg_to.left, al.seg_to.size() + k))
                    << ctx << " chain " << al.seg_from.left << "-" << al.seg_from.right << " on edge " << edge.getId()
                    << " [" << al.seg_to.left << "," << al.seg_to.right << ")";
            for (size_t p = al.seg_from.left; p < al.seg_from.right; p++)
                covered[p] = true;
        }
        return covered;
    }
}

TEST(DBGAligner, CarefulAlignCoversContigsStartingInsideEdges) {
    const size_t k = 11;
    size_t mid_edge_starts = 0;
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(600, seed), k);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        dbg::KmerIndex index(rdbg.dbg);
        index.fillAnchors(QuietLogger(), 1, rdbg.dbg, 5);
        for (size_t start = 0; start + 150 <= rdbg.genome.size(); start += 37) {
            const std::string ctx = "seed=" + std::to_string(seed) + " start=" + std::to_string(start);
            Contig contig(Sequence(rdbg.genome.substr(start, 150)), "c" + std::to_string(start));
            std::vector<ag::AlignmentChain<Contig, ag::Edge>> als = index.carefulAlign(contig);
            ASSERT_FALSE(als.empty()) << ctx;
            std::vector<bool> covered = CheckChains(contig, als, k, ctx);
//            Every k-mer of the contig except the last one is followed by a nucleotide, so it starts a covered step.
            for (size_t p = 0; p + k < contig.truncSize(); p++)
                EXPECT_TRUE(covered[p]) << ctx << " position " << p << " not covered";
            if (als.front().seg_to.left > 0)
                mid_edge_starts++;
        }
    }
    EXPECT_GT(mid_edge_starts, 0u) << "some contigs must start inside an edge";
}
