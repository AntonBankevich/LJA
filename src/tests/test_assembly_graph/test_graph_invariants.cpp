#include "gtest/gtest.h"
#include "graph_test_utils.hpp"
#include "random_graph_utils.hpp"
#include "supregraph/supregraph.hpp"
#include <sstream>

// Tests for AssemblyGraph::checkConsistency: valid graphs of every shape pass, basic edits keep them
// valid, and deliberately broken graphs are reported with the expected message.

using namespace ag_test;

namespace {
    // Runs checkConsistency, returning whether it passed and what it reported.
    std::pair<bool, std::string> Check(ag::AssemblyGraph &graph) {
        std::stringstream ss;
        bool ok = graph.checkConsistency(ss);
        return {ok, ss.str()};
    }

    ::testing::AssertionResult Consistent(ag::AssemblyGraph &graph) {
        auto res = Check(graph);
        if (res.first)
            return ::testing::AssertionSuccess();
        return ::testing::AssertionFailure() << res.second;
    }

    ::testing::AssertionResult ReportsViolation(ag::AssemblyGraph &graph, const std::string &message) {
        auto res = Check(graph);
        if (res.first)
            return ::testing::AssertionFailure() << "graph unexpectedly consistent, expected: " << message;
        if (res.second.find(message) == std::string::npos)
            return ::testing::AssertionFailure() << "expected \"" << message << "\" in report:\n" << res.second;
        return ::testing::AssertionSuccess();
    }

    // Small supregraph-shaped graph: u -suffix-> core -prefix-> w1, w2.
    struct SPGFixture {
        ag::AssemblyGraph graph;
        ag::Vertex &u = graph.addVertex(Sequence("TTAACCG"));
        ag::Vertex &core = graph.addVertex(Sequence("AACCG"));
        ag::Vertex &w1 = graph.addVertex(Sequence("AACCGA"));
        ag::Vertex &w2 = graph.addVertex(Sequence("AACCGT"));
        ag::Edge &inc = graph.addSPEdge(u, core);
        ag::Edge &e1 = graph.addSPEdge(core, w1);
        ag::Edge &e2 = graph.addSPEdge(core, w2);
    };
}

TEST(GraphInvariants, EmptyGraph) {
    ag::AssemblyGraph graph;
    EXPECT_TRUE(Consistent(graph));
}

TEST(GraphInvariants, ChainAndBasicEdits) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(5);
    EXPECT_TRUE(Consistent(gb.graph));

    gb.graph.removeEdge(*chain[1]);
    EXPECT_TRUE(Consistent(gb.graph));

    ag::Vertex &v = chain[3]->getFinish();
    gb.graph.isolateAndMark(v);
    EXPECT_TRUE(Consistent(gb.graph)); // marked vertices are allowed while isolated
    gb.graph.removeMarked();
    EXPECT_TRUE(Consistent(gb.graph));
    EXPECT_EQ(gb.graph.size(), 8u);
}

TEST(GraphInvariants, PalindromicVertexAndSelfRCEdge) {
    ag::AssemblyGraph graph;
    ag::Vertex &a = graph.addVertex(Sequence("AACGT"));
    ag::Vertex &b = graph.addVertex(Sequence("GACGT"));
    ag::Vertex &p = graph.addVertex(Sequence("ACGT"));
    ASSERT_EQ(p, p.rc());
    ag::Edge &self_rc = graph.addEdge(a, a.rc(), Sequence("AACGTT"));
    ASSERT_EQ(self_rc, self_rc.rc());
    ag::Edge &into_p = graph.addEdge(b, p, Sequence("GACGT")); // p is a suffix of b
    ASSERT_TRUE(into_p.isSuffix());
    EXPECT_TRUE(Consistent(graph));
}

TEST(GraphInvariants, SupregraphShape) {
    SPGFixture f;
    EXPECT_TRUE(Consistent(f.graph));
}

TEST(GraphInvariants, AddVertexPairKeepsIdOrder) {
    ag::AssemblyGraph graph;
    graph.addVertexPair({});
    graph.addVertex(Sequence("AACCG"));
    graph.addVertexPair({});
    graph.addSelfRCVertex({});
    graph.addVertexPair({});
    std::vector<int> ids;
    for (ag::Vertex &v : graph.vertices())
        ids.push_back(v.getInnerId());
    EXPECT_EQ(ids, std::vector<int>({1, -1, 2, -2, 3, -3, 4, 5, -5}));
    // Vertices made by addVertexPair have no sequence yet, but the list order itself must be fine.
    EXPECT_EQ(Check(graph).second.find("not sorted"), std::string::npos);
}

TEST(GraphInvariantsDeathTest, AddVertexPairRejectsOutOfOrderId) {
    ag::AssemblyGraph graph;
    graph.addVertexPair({});
    graph.addVertexPair({});
    EXPECT_DEATH(graph.addVertexPair({}, 1), "");
}

TEST(GraphInvariants, RandomDBGAndSupregraph) {
    RandomDBG rdbg(RandomGenome(1000, 17), 5);
    EXPECT_TRUE(Consistent(rdbg.dbg));
    ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
    EXPECT_TRUE(Consistent(rdbg.dbg));
    ConvertToSupregraph(rdbg.dbg);
    EXPECT_TRUE(Consistent(rdbg.dbg));
}

TEST(GraphInvariants, DetectsSuffixEdgeWithSiblings) {
    SPGFixture f;
    ag::Vertex &x = f.graph.addVertex(Sequence("AACCGTA"));
    f.graph.addEdge(f.u, x, f.u.getSeq() + x.getSeq());
    EXPECT_TRUE(ReportsViolation(f.graph, "has an outgoing suffix edge and other outgoing edges"));
}

TEST(GraphInvariants, DetectsSharedFirstCodeNucleotide) {
    GraphBuilder gb(8);
    ag::Vertex &a = gb.newVertex(); // AAAAAAAA
    ag::Vertex &b = gb.newVertex(); // CAAAAAAA
    gb.newVertex();                 // ACAAAAAA (unused)
    ag::Vertex &d = gb.newVertex(); // CCAAAAAA
    gb.connect(a, b);
    gb.connect(a, d); // a->b and a->d both start with C
    EXPECT_TRUE(ReportsViolation(gb.graph, "has the same first nucleotide"));
}

TEST(GraphInvariants, DetectsMarkingProblems) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(3);
    ag::Vertex &v = chain[0]->getFinish();
    v.mark();
    EXPECT_TRUE(ReportsViolation(gb.graph, "marked flag differs from its rc"));
    v.rc().mark();
    EXPECT_TRUE(ReportsViolation(gb.graph, "is marked but still has outgoing edges"));
    v.unmark();
    v.rc().unmark();
    EXPECT_TRUE(Consistent(gb.graph));
}

TEST(GraphInvariants, DetectsNonCorporealEdge) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(3);
    chain[0]->setCorporeal(false);
    EXPECT_TRUE(ReportsViolation(gb.graph, "is not corporeal"));
    chain[0]->setCorporeal(true);
    EXPECT_TRUE(Consistent(gb.graph));
}

// swapDirection turns a path into exactly its RC() in place, and twice gives the path back.
TEST(GraphInvariants, SwapDirectionIsInPlaceRC) {
    RandomDBG rdbg(RandomGenome(600, 3), 5);
    std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 50, 1, 7);
    for (ag::AlignedRead &read : reads) {
        ag::GraphPath path = read.getPath();
        ag::GraphPath swapped = path;
        swapped.swapDirection();
        EXPECT_TRUE(swapped == path.RC()) << read.getId();
        swapped.swapDirection();
        EXPECT_TRUE(swapped == path) << read.getId();
    }
}

namespace {
//    The sub-path built edge by edge, as subPath did before it sliced the code decks.
    ag::GraphPath ReferenceSubPath(const ag::GraphPath &path, ag::PathPosition from, ag::PathPosition to) {
        ag::GraphPath res(from.getVertex());
        for (ag::PathPosition cur = from; cur != to; ++cur)
            res += cur.nextEdge();
        if (from == path.firstPosition())
            res.setCutLeft(path.leftCut());
        if (to == path.lastPosition())
            res.setCutRight(path.rightCut());
        return res;
    }

    void CheckSubPaths(ag::AlignedReadStorage &storage, const std::string &ctx) {
        size_t checked = 0;
        for (ag::AlignedRead &read : storage) {
            const ag::GraphPath &path = read.getPath();
            if (!read.valid() || path.empty())
                continue;
            std::vector<ag::PathPosition> positions;
            for (ag::PathPosition p = path.firstPosition(); ; ++p) {
                positions.push_back(p);
                if (p == path.lastPosition())
                    break;
            }
            for (size_t i = 0; i < positions.size(); i++)
                for (size_t j = i; j < positions.size(); j++) {
                    EXPECT_TRUE(path.subPath(positions[i], positions[j]) == ReferenceSubPath(path, positions[i], positions[j]))
                            << ctx << " read " << read.getId() << " " << i << "-" << j;
                    checked++;
                }
        }
        EXPECT_GT(checked, 0u) << ctx;
    }
}

// subPath slices the code decks; it must equal the edge-by-edge construction on a DBG and on a supregraph (where
// prefix/suffix edges have their code in only one of the decks).
TEST(GraphInvariants, SubPathMatchesEdgeByEdgeConstruction) {
    RandomDBG rdbg(RandomGenome(400, 5), 5);
    ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, 30, 40, 1, 9));
    ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
    CheckSubPaths(storage, "dbg");
    ConvertToSupregraph(rdbg.dbg);
    CheckSubPaths(storage, "spg");
}

// extend walks the sequence through the graph, comparing the rest of an edge in one go once inside it: extending from
// the genome's first k-mer spells the genome, and a single mismatch invalidates the path.
TEST(GraphInvariants, ExtendSpellsTheSequence) {
    const size_t k = 11;
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(500, seed), k);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        ag::Vertex *start = nullptr;
        for (ag::Vertex &v : rdbg.dbg.vertices())
            if (v.getSeq().str() == rdbg.genome.substr(0, k))
                start = &v;
        ASSERT_NE(start, nullptr);
        ag::GraphPath path(*start);
        path.extend(Sequence(rdbg.genome.substr(k)));
        ASSERT_TRUE(path.valid()) << "seed=" << seed;
        EXPECT_EQ(path.Seq().str(), rdbg.genome) << "seed=" << seed;
        std::string wrong = rdbg.genome.substr(k);
        wrong[wrong.size() / 2] = wrong[wrong.size() / 2] == 'A' ? 'C' : 'A';
        ag::GraphPath bad(*start);
        bad.extend(Sequence(wrong));
        EXPECT_FALSE(bad.valid()) << "seed=" << seed;
    }
}

// A read path may end with a suffix edge (empty trunc) while its right cut lies inside the last vertex; truncSeq then cuts
// into the earlier edges. Mirrored for a prefix edge at the start.
TEST(GraphPathTest, TruncSeqWithCutsBeyondEndEdges) {
    ag::AssemblyGraph graph(2);
    ag::Vertex &x = graph.addVertex(Sequence("TTA"));
    ag::Vertex &u = graph.addVertex(Sequence("TACGG"));
    ag::Vertex &v = graph.addVertex(Sequence("CGG"));
    ag::Edge &xu = graph.addEdge(x, u, Sequence("TTACGG"));
    ag::Edge &uv = graph.addSPEdge(u, v);
    ASSERT_TRUE(uv.isSuffix());
    ag::GraphPath path(xu);
    path += uv;
    path.setCutRight(2);
    EXPECT_EQ(path.truncSeq().str(), "C");
    EXPECT_EQ(path.truncSeq().size(), path.truncLen());

    ag::GraphPath rc = path.RC();
    ASSERT_TRUE(rc.frontEdge().isPrefix());
    EXPECT_EQ(rc.truncSeq().size(), rc.truncLen());
    EXPECT_EQ(rc.Seq().str(), (!path.Seq()).str());
}

// Supregraph: an edge from a vertex to itself is allowed only as the vertex's whole connected component (a cycle
// collapsed into one vertex); next to any other edge it is reported.
TEST(GraphInvariants, SupregraphSelfLoopOnlyAsWholeComponent) {
    ag::AssemblyGraph alone_base(2);
    ag::Vertex &v = alone_base.addVertex(Sequence("ACG"));
    alone_base.addEdge(v, v, Sequence("ACGTTACG"));
    spg::Supregraph alone(std::move(alone_base));
    EXPECT_TRUE(Consistent(alone));

    ag::AssemblyGraph with_exit_base(2);
    ag::Vertex &a = with_exit_base.addVertex(Sequence("ACG"));
    ag::Vertex &b = with_exit_base.addVertex(Sequence("CGA"));
    with_exit_base.addEdge(a, a, Sequence("ACGTTACG"));
    with_exit_base.addEdge(a, b, Sequence("ACGA"));
    spg::Supregraph with_exit(std::move(with_exit_base));
    EXPECT_TRUE(ReportsViolation(with_exit, "to itself, but the vertex has other edges"));
}
