#include "gtest/gtest.h"
#include "graph_test_utils.hpp"
#include "random_graph_utils.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"

using namespace ag_test;

namespace {
//    Adds the resolution sample [start, finish) to v and the mirrored one to v.rc(), as SampleSegmentTracker keeps both
//    orientations.
    void AddSample(ag::Vertex &v, size_t start, size_t finish) {
        v.coverage_info.addResolutionSample(start, finish, 1, 1);
        v.rc().coverage_info.addResolutionSample(v.size() - finish, v.size() - start, 1, 1);
    }
}

TEST(AlignedReadStorageTest, BasicAddAndConsistency) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(6); // 5 edges

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("read0", MakePath(chain, 0, 3)); // edges 0,1,2
    reads.emplace_back("read1", MakePath(chain, 2, 5)); // edges 2,3,4

    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    EXPECT_TRUE(storage.checkConsistency());

    // read0 starts at edge 0 and ends at edge 2; read1 starts at edge 2 and ends at edge 4.
    EXPECT_EQ(storage.getOutgoingReads(chain[0]->getStart()).size(), 1u);
    EXPECT_EQ(storage.getOutgoingReads(chain[2]->getStart()).size(), 1u);
    // Both reads end within edges 2 and 4 respectively - i.e. rc() of those edges is where the
    // reverse-direction reads start.
    EXPECT_EQ(storage.getOutgoingReads(chain[2]->rc().getStart()).size(), 1u);
    EXPECT_EQ(storage.getOutgoingReads(chain[4]->rc().getStart()).size(), 1u);
}

TEST(AlignedReadStorageTest, InvalidateReadDropsStartCounts) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(4); // 3 edges

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("read0", MakePath(chain, 0, 3));

    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency());
    ASSERT_EQ(storage.getOutgoingReads(chain[0]->getStart()).size(), 1u);

    ag::AlignedRead &read = storage[0];
    storage.startReadCorrection();
    storage.delayedInvalidateRead(read, "test invalidate");
    storage.applyCorrections(TestLogger(), 1);

    EXPECT_FALSE(read.valid());
    EXPECT_EQ(storage.getOutgoingReads(chain[0]->getStart()).size(), 0u);
    EXPECT_TRUE(storage.checkConsistency());
}

// A read with an empty path (a substring of one vertex) is filed in reads_inside_vertices under that vertex
// (forward) and its rc (backward), not in starts.
TEST(AlignedReadStorageTest, SubreadIsFiledInsideItsVertex) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(3);
    ag::Vertex &v = chain[0]->getFinish();
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("sub", ag::GraphPath(v, 2, 3));
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));

    EXPECT_TRUE(storage.checkConsistency());
    EXPECT_TRUE(gb.graph.checkConsistency());
    EXPECT_EQ(storage.getSubstringReads(v.getId()).size(), 1u);
    EXPECT_EQ(storage.getSubstringReads(v.rc().getId()).size(), 1u);
    EXPECT_EQ(storage.getOutgoingReads(v).size(), 0u);
}

// Deleting a vertex turns its subreads into legacy reads: the path becomes a LegacyPath that keeps the cuts, and
// the read moves to legacy_reads (checked by checkConsistency).
TEST(AlignedReadStorageTest, DeletedVertexTurnsSubreadIntoLegacy) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(3);
    ag::Vertex &v = chain[0]->getFinish();
    ag::Vertex::id_type vid = v.getInnerId();
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("sub", ag::GraphPath(v, 2, 3));
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency());

    gb.graph.isolateAndMark(v);
    gb.graph.removeMarked();

    const ag::GraphPath &path = storage[0].getPath();
    EXPECT_TRUE(path.isLegacy());
    EXPECT_EQ(path.getStartId().innerId(), vid);
    EXPECT_EQ(path.leftCut(), 2u);
    EXPECT_EQ(path.rightCut(), 3u);
    EXPECT_TRUE(storage.checkConsistency());
    EXPECT_TRUE(gb.graph.checkConsistency());
}

// mergePath absorbs the inner vertices of the path into the new vertex res. A read that was a subread of an inner
// vertex is a substring of res, so it must become a subread of res with the same sequence -- not a legacy read when
// the inner vertex is deleted.
TEST(AlignedReadStorageTest, MergePathKeepsSubreadOfInnerVertex) {
    ag::AssemblyGraph graph(2);
    ag::Vertex &s = graph.addVertex(Sequence("TAAC"));
    ag::Vertex &a = graph.addVertex(Sequence("AAC"));
    ag::Vertex &f = graph.addVertex(Sequence("AACG"));
    ag::Edge &sa = graph.addSPEdge(s, a);
    ag::Edge &af = graph.addSPEdge(a, f);
    AddSample(a, 0, 3);
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("sub", ag::GraphPath(a, 0, 0));
    ag::AlignedReadStorage storage(graph, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency());
    ag::GraphPath merged(sa);
    merged += af;

    ag::Vertex &res = graph.mergePath(merged);
    graph.removeMarked();

    const ag::GraphPath &path = storage[0].getPath();
    ASSERT_TRUE(path.valid() && !path.isLegacy());
    EXPECT_TRUE(path.empty());
    EXPECT_EQ(path.getStart(), res);
    EXPECT_EQ(path.Seq().str(), "AAC");
    EXPECT_TRUE(storage.checkConsistency());
    EXPECT_TRUE(graph.checkConsistency());
}

// Control: a subread of the path's start vertex is not touched by mergePath (start survives).
TEST(AlignedReadStorageTest, ReadOnSurvivingVertexUnaffectedByMergePath) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(6);
    ag::Vertex &start = chain[1]->getStart(); // v1, start of v1 -> v4
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("sub", ag::GraphPath(start, 2, 3));
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency());
    std::string seq = storage[0].getPath().Seq().str();

    gb.graph.mergePath(MakePath(chain, 1, 4));
    gb.graph.removeMarked();

    const ag::GraphPath &path = storage[0].getPath();
    ASSERT_TRUE(path.valid() && !path.isLegacy());
    EXPECT_EQ(path.getStart(), start);
    EXPECT_EQ(path.Seq().str(), seq);
    EXPECT_TRUE(storage.checkConsistency());
    EXPECT_TRUE(gb.graph.checkConsistency());
}

// The storage constructor normalizes read paths; without sample segments it must not change them, cuts included.
TEST(AlignedReadStorageTest, ConstructorKeepsCuts) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(4);
    ag::GraphPath path = MakePath(chain, 0, 3);
    path.setCutLeft(3);
    path.setCutRight(2);
    std::string seq = path.Seq().str();
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r", path);
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));

    EXPECT_EQ(storage[0].getPath().leftCut(), 3u);
    EXPECT_EQ(storage[0].getPath().rightCut(), 2u);
    EXPECT_EQ(storage[0].getPath().Seq().str(), seq);
    EXPECT_TRUE(storage.checkConsistency());
}

// Normalization extends a read path through unbranching steps of any type to every vertex whose first sample segment
// lies inside the read. Hand-built SPGs: k = 3, every read's end vertices have a segment inside the read.
// Suffix edge: the read TAAC is a subread of u = TAACG; u -> w = AACG.
TEST(AlignedReadStorageTest, NormalizeExtendsThroughSuffixEdge) {
    for (bool segment_inside : {true, false}) {
        ag::AssemblyGraph graph(3);
        ag::Vertex &u = graph.addVertex(Sequence("TAACG"));
        ag::Vertex &w = graph.addVertex(Sequence("AACG"));
        graph.addSPEdge(u, w);
        AddSample(u, 0, 2);                   // TA
        AddSample(w, segment_inside ? 0 : 1, segment_inside ? 2 : 4); // AA inside / ACG past the read
        std::vector<ag::AlignedRead> reads;
        reads.emplace_back("r", ag::GraphPath(u, 0, 1)); // TAAC
        ag::AlignedReadStorage storage(graph, std::move(reads));

        const ag::GraphPath &path = storage[0].getPath();
        EXPECT_EQ(path.Seq().str(), "TAAC");
        EXPECT_EQ(path.getFinish(), segment_inside ? w : u);
        EXPECT_TRUE(storage.checkConsistency());
        EXPECT_TRUE(graph.checkConsistency());
    }
}

// Prefix edge, on the left: the read ACGT is a subread of u = AACGT; p = AACG -> u.
TEST(AlignedReadStorageTest, NormalizeExtendsThroughPrefixEdgeAtStart) {
    for (bool segment_inside : {true, false}) {
        ag::AssemblyGraph graph(3);
        ag::Vertex &p = graph.addVertex(Sequence("AACG"));
        ag::Vertex &u = graph.addVertex(Sequence("AACGT"));
        graph.addSPEdge(p, u);
        AddSample(u, 1, 5);                   // ACGT
        AddSample(p, segment_inside ? 1 : 0, segment_inside ? 3 : 2); // AC inside / AA past the read
        std::vector<ag::AlignedRead> reads;
        reads.emplace_back("r", ag::GraphPath(u, 1, 0)); // ACGT
        ag::AlignedReadStorage storage(graph, std::move(reads));

        const ag::GraphPath &path = storage[0].getPath();
        EXPECT_EQ(path.Seq().str(), "ACGT");
        EXPECT_EQ(path.getStart(), segment_inside ? p : u);
        EXPECT_TRUE(storage.checkConsistency());
        EXPECT_TRUE(graph.checkConsistency());
    }
}

// Prefix edge at the end: the read AACG is u = AACG; u -> w = AACGT, w's segment AA lies inside the read.
TEST(AlignedReadStorageTest, NormalizeExtendsThroughPrefixEdgeAtEnd) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &u = graph.addVertex(Sequence("AACG"));
    ag::Vertex &w = graph.addVertex(Sequence("AACGT"));
    graph.addSPEdge(u, w);
    AddSample(u, 0, 4);
    AddSample(w, 0, 2);
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r", ag::GraphPath(u, 0, 0));
    ag::AlignedReadStorage storage(graph, std::move(reads));

    EXPECT_EQ(storage[0].getPath().getFinish(), w);
    EXPECT_EQ(storage[0].getPath().Seq().str(), "AACG");
    EXPECT_TRUE(storage.checkConsistency());
}

// A path that is not in normal form is reported. Both constructors normalize, so the path is put back afterwards: a
// subread of u = TAACG although w's segment lies inside the read.
TEST(AlignedReadStorageTest, PathNotInNormalFormIsReported) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &u = graph.addVertex(Sequence("TAACG"));
    ag::Vertex &w = graph.addVertex(Sequence("AACG"));
    graph.addSPEdge(u, w);
    AddSample(u, 0, 2);
    AddSample(w, 0, 2);
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r", ag::GraphPath(u, 0, 1)); // TAAC
    ag::AlignedReadStorage storage(TestLogger(), 1, graph, std::move(reads));
    storage[0].getPath() = ag::GraphPath(u, 0, 1);

    std::stringstream ss;
    EXPECT_FALSE(storage.checkConsistency(ss));
    EXPECT_NE(ss.str().find("is not in normal form"), std::string::npos) << ss.str();
}

// Normalization passes through vertices without sample segments: u -> x (no segments) -> y, where y's segment lies
// inside the read (a subread of u).
TEST(AlignedReadStorageTest, NormalizeStepsOverSegmentlessVertex) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &u = graph.addVertex(Sequence("TTAACG"));
    ag::Vertex &x = graph.addVertex(Sequence("TAACG"));
    ag::Vertex &y = graph.addVertex(Sequence("AACG"));
    graph.addSPEdge(u, x);
    graph.addSPEdge(x, y);
    AddSample(u, 0, 2);
    AddSample(y, 0, 2);
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r", ag::GraphPath(u, 0, 1)); // TTAAC
    ag::AlignedReadStorage storage(graph, std::move(reads));

    EXPECT_EQ(storage[0].getPath().getFinish(), y);
    EXPECT_EQ(storage[0].getPath().Seq().str(), "TTAAC");
    EXPECT_TRUE(storage.checkConsistency());
}

// An SPG path end must have a sample segment inside the read.
TEST(AlignedReadStorageTest, EndWithoutSegmentInsideReadIsReported) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &u = graph.addVertex(Sequence("TAACG"));
    AddSample(u, 3, 5); // CG, outside the read
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r", ag::GraphPath(u, 0, 2)); // TAA
    ag::AlignedReadStorage storage(TestLogger(), 1, graph, std::move(reads));

    std::stringstream ss;
    EXPECT_FALSE(storage.checkConsistency(ss));
    EXPECT_NE(ss.str().find("has no sample segment inside the read"), std::string::npos) << ss.str();
}

namespace {
//    A small DBG (k = 3) around e = u -> w (a merged edge, AACGGTA): t -> u -> w -> x, plus u -> y. No vertex is its
//    own or another's rc. The reads are minimal; those ending inside e have cuts there, so the reaction must keep cuts.
    struct EdgeToSupreVertexFixture {
        ag::AssemblyGraph graph{3};
        ag::Vertex &t = graph.addVertex(Sequence("TAA"));
        ag::Vertex &u = graph.addVertex(Sequence("AAC"));
        ag::Vertex &w = graph.addVertex(Sequence("GTA"));
        ag::Vertex &x = graph.addVertex(Sequence("TAG"));
        ag::Vertex &y = graph.addVertex(Sequence("ACT"));
        ag::Edge &tu = graph.addEdge(t, u, Sequence("TAAC"));
        ag::Edge &e = graph.addEdge(u, w, Sequence("AACGGTA"));
        ag::Edge &wx = graph.addEdge(w, x, Sequence("GTAG"));
        ag::Edge &uy = graph.addEdge(u, y, Sequence("AACT"));

        ag::GraphPath path(std::initializer_list<ag::Edge *> edges, size_t cut_left = 0, size_t cut_right = 0) {
            ag::GraphPath res(**edges.begin());
            for (auto it = edges.begin() + 1; it != edges.end(); ++it)
                res += **it;
            res.setCutLeft(cut_left);
            res.setCutRight(cut_right);
            return res;
        }

        std::vector<ag::AlignedRead> reads() {
            std::vector<ag::AlignedRead> res;
            res.emplace_back("through", path({&tu, &e, &wx}));
            res.emplace_back("starts", path({&e, &wx}, 2, 0));
            res.emplace_back("ends", path({&tu, &e}, 0, 2));
            res.emplace_back("single", path({&e}, 1, 2));
            res.emplace_back("other", path({&uy}));
            return res;
        }
    };

    std::map<std::string, std::string> SeqsById(ag::AlignedReadStorage &storage) {
        std::map<std::string, std::string> res;
        for (const ag::AlignedRead &read : storage)
            res[read.getId()] = read.getPath().Seq().str();
        return res;
    }

    const ag::AlignedRead &ReadById(ag::AlignedReadStorage &storage, const std::string &id) {
        for (const ag::AlignedRead &read : storage)
            if (read.getId() == id)
                return read;
        VERIFY(false);
        return *storage.begin();
    }
}

// Initial conversion: u and w are k-mers without sample segments, so paths starting/ending with e drop u/w and
// start/end at the new supre-vertex v; a path that is e alone becomes a subread of v. Other reads keep their ends.
TEST(AlignedReadStorageTest, EdgeToSupreVertexDropsSegmentlessEnds) {
    EdgeToSupreVertexFixture f;
    ag::AlignedReadStorage storage(f.graph, f.reads());
    ASSERT_TRUE(storage.checkConsistency());
    std::map<std::string, std::string> before = SeqsById(storage);

    ag::Vertex &v = f.graph.edgeToSupreVertex(f.e);
    f.graph.removeMarked();

    EXPECT_EQ(SeqsById(storage), before);
    EXPECT_EQ(ReadById(storage, "through").getPath().getStart(), f.t);
    EXPECT_EQ(ReadById(storage, "through").getPath().getFinish(), f.x);
    EXPECT_EQ(ReadById(storage, "starts").getPath().getStart(), v);
    EXPECT_EQ(ReadById(storage, "starts").getPath().getFinish(), f.x);
    EXPECT_EQ(ReadById(storage, "ends").getPath().getStart(), f.t);
    EXPECT_EQ(ReadById(storage, "ends").getPath().getFinish(), v);
    EXPECT_TRUE(ReadById(storage, "single").getPath().empty());
    EXPECT_EQ(ReadById(storage, "single").getPath().getStart(), v);
    EXPECT_EQ(ReadById(storage, "other").getPath().getStart(), f.u);
    EXPECT_EQ(ReadById(storage, "other").getPath().getFinish(), f.y);
    EXPECT_TRUE(storage.checkConsistency());
    EXPECT_TRUE(f.graph.checkConsistency());
}

// mergePrefixPath removes the inner vertices of a prefix run; a read starting at one of them (it is a prefix of the
// finish) becomes a read starting at the finish, which receives the inner vertex's samples. Here it becomes a subread.
TEST(AlignedReadStorageTest, MergePrefixPathMovesReadsFromInnerVertices) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &v0 = graph.addVertex(Sequence("AAC"));
    ag::Vertex &v1 = graph.addVertex(Sequence("AACG"));
    ag::Vertex &v2 = graph.addVertex(Sequence("AACGT"));
    ag::Edge &e0 = graph.addSPEdge(v0, v1);
    ag::Edge &e1 = graph.addSPEdge(v1, v2);
    AddSample(v1, 1, 4); // ACG
    AddSample(v2, 2, 5); // CGT
    std::vector<ag::AlignedRead> reads;
    ag::GraphPath read(e1);
    read.setCutLeft(1); // ACGT
    reads.emplace_back("r", read);
    ag::AlignedReadStorage storage(graph, std::move(reads));
    ASSERT_EQ(storage[0].getPath().getStart(), v1);
    ASSERT_TRUE(storage.checkConsistency());
    ag::GraphPath run(e0);
    run += e1;

    graph.mergePrefixPath(run);
    graph.removeMarked();

    ASSERT_TRUE(storage.checkConsistency());
    EXPECT_EQ(storage[0].getPath().getStart(), v2);
    EXPECT_TRUE(storage[0].getPath().empty());
    EXPECT_EQ(storage[0].getPath().Seq().str(), "ACGT");
    EXPECT_TRUE(graph.checkConsistency());
}

// flipVertex removes core, a suffix of its single in-neighbour u: a read starting at core starts at u instead, its left
// cut grown by |u| - |core|, and leaves u by the new edge that replaces core's.
TEST(AlignedReadStorageTest, FlipMovesReadsFromCore) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &u = graph.addVertex(Sequence("TTAACG"));
    ag::Vertex &core = graph.addVertex(Sequence("AACG"));
    ag::Vertex &w1 = graph.addVertex(Sequence("AACGA"));
    ag::Vertex &w2 = graph.addVertex(Sequence("AACGT"));
    graph.addSPEdge(u, core);
    ag::Edge &cw1 = graph.addSPEdge(core, w1);
    graph.addSPEdge(core, w2);
    AddSample(u, 0, 2);  // TT
    AddSample(core, 1, 3); // AC
    AddSample(w1, 2, 5); // CGA
    AddSample(w2, 2, 5); // CGT
    std::vector<ag::AlignedRead> reads;
    ag::GraphPath read(cw1);
    read.setCutLeft(1); // ACGA
    reads.emplace_back("r", read);
    ag::AlignedReadStorage storage(graph, std::move(reads));
    ASSERT_EQ(storage[0].getPath().getStart(), core);
    ASSERT_TRUE(storage.checkConsistency());

    graph.flipVertexLockFree(core.incFront());
    graph.removeMarked();

    ASSERT_TRUE(storage.checkConsistency());
    EXPECT_EQ(storage[0].getPath().getStart(), u);
    EXPECT_EQ(storage[0].getPath().getFinish(), w1);
    EXPECT_EQ(storage[0].getPath().leftCut(), 3u);
    EXPECT_EQ(storage[0].getPath().Seq().str(), "ACGA");
    EXPECT_TRUE(graph.checkConsistency());
}

namespace {
//    u -suffix-> core -prefix-> w1, w2 and w1 -> y1, y2, so that on the rc strand both core.rc and w1.rc have two
//    incoming edges and carry suffix records. The read r = core -> w1 -> y1 (ACGACT) ends, on the rc strand, with the
//    code-less suffix edge w1.rc -> core.rc.
    struct FlipWithTrackerFixture {
        ag::AssemblyGraph graph{3};
        ag::Vertex &u = graph.addVertex(Sequence("TTAACG"));
        ag::Vertex &core = graph.addVertex(Sequence("AACG"));
        ag::Vertex &w1 = graph.addVertex(Sequence("AACGA"));
        ag::Vertex &w2 = graph.addVertex(Sequence("AACGT"));
        ag::Vertex &y1 = graph.addVertex(Sequence("CGACT"));
        ag::Vertex &y2 = graph.addVertex(Sequence("CGAGT"));
        std::unique_ptr<ag::AlignedReadStorage> storage;
        std::unique_ptr<ag::SuffixTracker> tracker;

        FlipWithTrackerFixture() {
            graph.addSPEdge(u, core);
            ag::Edge &cw1 = graph.addSPEdge(core, w1);
            graph.addSPEdge(core, w2);
            ag::Edge &w1y1 = graph.addEdge(w1, y1, Sequence("AACGACT"));
            graph.addEdge(w1, y2, Sequence("AACGAGT"));
            AddSample(u, 0, 2);   // TT
            AddSample(core, 1, 3); // AC
            AddSample(w1, 2, 5);  // CGA
            AddSample(w2, 2, 5);  // CGT
            AddSample(y1, 2, 5);  // ACT
            AddSample(y2, 2, 5);  // AGT
            std::vector<ag::AlignedRead> reads;
            ag::GraphPath read(cw1);
            read += w1y1;
            read.setCutLeft(1); // ACGACT
            reads.emplace_back("r", read);
            storage.reset(new ag::AlignedReadStorage(graph, std::move(reads)));
            tracker.reset(new ag::SuffixTracker(*storage, graph, 0, 1000));
            tracker->fillFromStorage(QuietLogger(), 1);
        }
    };
}

// The read is extended to u. On the rc strand it then ends with the coded edge
// w1.rc -> u.rc instead of the code-less w1.rc -> core.rc, so its suffix after y1.rc -> w1.rc grows and the suffix
// tracker, which only follows graph events, must still match it: a trailing virtually suffix edge does not count for the
// canonical form of a stored suffix.
TEST(AlignedReadStorageTest, FlipExtendingReadKeepsSuffixTracker) {
    FlipWithTrackerFixture f;
    ASSERT_EQ(f.storage->operator[](0).getPath().getStart(), f.core);
    ASSERT_TRUE(f.storage->checkConsistency(std::cerr));
    std::string seq = f.storage->operator[](0).getPath().Seq().str();
    f.graph.flipVertexLockFree(f.core.incFront());
    f.graph.removeMarked();
    EXPECT_EQ(f.storage->operator[](0).getPath().getStart(), f.u);
    EXPECT_EQ(f.storage->operator[](0).getPath().Seq().str(), seq);
    EXPECT_TRUE(f.storage->checkConsistency(std::cerr));
    EXPECT_TRUE(f.graph.checkConsistency());
}

// As above, but on the rc strand the read's suffix after the record edge y1.rc -> a.rc is a coded virtually suffix edge
// a.rc -> w1.rc followed by the code-less suffix edge w1.rc -> core.rc. After the flip the latter becomes the coded
// virtually suffix edge w1.rc -> u.rc. Dropping only the last edge's code from the canonical form would leave a.rc ->
// w1.rc's first letter in it after the flip but not before; all trailing virtually suffix edges must be dropped.
TEST(AlignedReadStorageTest, FlipExtendingReadPastVirtualSuffixEdgeKeepsSuffixTracker) {
    ag::AssemblyGraph graph(3);
    ag::Vertex &u = graph.addVertex(Sequence("TTAACG"));
    ag::Vertex &core = graph.addVertex(Sequence("AACG"));
    ag::Vertex &w1 = graph.addVertex(Sequence("AACGA"));
    ag::Vertex &w2 = graph.addVertex(Sequence("AACGT"));
    ag::Vertex &a = graph.addVertex(Sequence("CGAC"));
    ag::Vertex &y1 = graph.addVertex(Sequence("GACT"));
    ag::Vertex &y2 = graph.addVertex(Sequence("GACG"));
    graph.addSPEdge(u, core);
    ag::Edge &cw1 = graph.addSPEdge(core, w1);
    graph.addSPEdge(core, w2);
//    w1 -> a is w1's only outgoing and a's only incoming edge; a forks into y1, y2.
    ag::Edge &w1a = graph.addEdge(w1, a, Sequence("AACGAC"));
    w1a.setVirtualType(ag::VirtualType::merge);
    ag::Edge &ay1 = graph.addEdge(a, y1, Sequence("CGACT"));
    ay1.setVirtualType(ag::VirtualType::prefix);
    graph.addEdge(a, y2, Sequence("CGACG")).setVirtualType(ag::VirtualType::prefix);
    AddSample(u, 0, 2);   // TT
    AddSample(core, 1, 3); // AC
    AddSample(w1, 2, 5);  // CGA
    AddSample(w2, 2, 5);  // CGT
    AddSample(a, 1, 4);   // GAC
    AddSample(y1, 1, 4);  // ACT
    AddSample(y2, 1, 4);  // ACG
    std::vector<ag::AlignedRead> reads;
    ag::GraphPath read(cw1);
    read += w1a;
    read += ay1;
    read.setCutLeft(1); // ACGACT
    reads.emplace_back("r", read);
    ag::AlignedReadStorage storage(graph, std::move(reads));
    ag::SuffixTracker tracker(storage, graph, 0, 1000);
    tracker.fillFromStorage(QuietLogger(), 1);
    ASSERT_EQ(storage[0].getPath().getStart(), core);
    ASSERT_TRUE(storage.checkConsistency(std::cerr));
    std::string seq = storage[0].getPath().Seq().str();

    graph.flipVertexLockFree(core.incFront());
    graph.removeMarked();

    EXPECT_EQ(storage[0].getPath().getStart(), u);
    EXPECT_EQ(storage[0].getPath().Seq().str(), seq);
    EXPECT_TRUE(storage.checkConsistency(std::cerr));
    EXPECT_TRUE(graph.checkConsistency());
}

// Unaligned reads (invalid paths, as AlignReads produces for reads too short to align) are kept as they are by both
// constructors, which normalize every other path.
TEST(AlignedReadStorageTest, ConstructorsKeepUnalignedReads) {
    for (bool parallel : {false, true}) {
        GraphBuilder gb(8);
        std::vector<ag::EdgeId> chain = gb.buildChain(4);
        std::vector<ag::AlignedRead> reads;
        reads.emplace_back("aligned", MakePath(chain, 0, 3));
        reads.emplace_back("unaligned", ag::GraphPath());
        std::unique_ptr<ag::AlignedReadStorage> storage(parallel ?
                new ag::AlignedReadStorage(TestLogger(), 2, gb.graph, std::move(reads)) :
                new ag::AlignedReadStorage(gb.graph, std::move(reads)));
        EXPECT_TRUE((*storage)[0].valid());
        EXPECT_FALSE((*storage)[1].valid());
        EXPECT_TRUE(storage->checkConsistency()) << "parallel=" << parallel;
    }
}

// mergePath of a loop path J -> a -> b -> J (start == finish; J has another in- and out-edge). Reads starting at the
// inner vertices a and b, including ones that go around the loop again, must keep their sequences: the absorption
// used to stop advancing along the path at its first vertex because that vertex is also the finish. Every vertex
// carries a whole-vertex sample segment, so uncut reads are in normal form.
TEST(AlignedReadStorageTest, MergeLoopPathKeepsReadsStartingInside) {
    ag::AssemblyGraph graph(2);
    ag::Vertex &j = graph.addVertex(Sequence("ACG"));
    ag::Vertex &a = graph.addVertex(Sequence("GTT"));
    ag::Vertex &b = graph.addVertex(Sequence("TCA"));
    ag::Vertex &x = graph.addVertex(Sequence("TTA"));
    ag::Vertex &y = graph.addVertex(Sequence("GAA"));
    ag::Edge &xj = graph.addEdge(x, j, Sequence("TTACG"));
    ag::Edge &ja = graph.addEdge(j, a, Sequence("ACGTT"));
    ag::Edge &ab = graph.addEdge(a, b, Sequence("GTTCA"));
    ag::Edge &bj = graph.addEdge(b, j, Sequence("TCACG"));
    ag::Edge &jy = graph.addEdge(j, y, Sequence("ACGAA"));
    auto make = [](std::initializer_list<ag::Edge *> edges, size_t cut_left, size_t cut_right) {
        ag::GraphPath res(**edges.begin());
        for (auto it = edges.begin() + 1; it != edges.end(); ++it)
            res += **it;
        res.setCutLeft(cut_left);
        res.setCutRight(cut_right);
        return res;
    };
    std::vector<ag::AlignedRead> reads;
    for (ag::Vertex *v : {&j, &a, &b, &x, &y})
        AddSample(*v, 0, v->size());
    reads.emplace_back("wraps", make({&ab, &bj, &ja, &ab}, 0, 0));
    reads.emplace_back("leaves", make({&ab, &bj, &jy}, 0, 0));
    reads.emplace_back("from_b", make({&bj, &ja}, 0, 0));
    reads.emplace_back("inside", make({&ab}, 0, 0));
    reads.emplace_back("enters", make({&xj, &ja, &ab}, 0, 0));
    reads.emplace_back("around", make({&ja, &ab, &bj}, 0, 0));
    ag::AlignedReadStorage storage(graph, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency(std::cerr));
    std::map<std::string, std::string> before = SeqsById(storage);

    ag::GraphPath loop = make({&ja, &ab, &bj}, 0, 0);
    graph.mergePath(loop);
    graph.removeMarked();

    EXPECT_EQ(SeqsById(storage), before);
    for (const ag::AlignedRead &read : storage)
        EXPECT_TRUE(read.getPath().valid() && !read.getPath().isLegacy()) << read.getId();
    EXPECT_TRUE(storage.checkConsistency());
    EXPECT_TRUE(graph.checkConsistency());
}

// resolveVertex on a supregraph core c (suffix edges in, prefix edges out). Reads consisting of the single edge
// c -> w collapse into subreads of w; one of them is stored in the opposite orientation, so the direction starting at
// c is its backward one. Its other direction must be moved from starts[w.rc()], not looked up via read.backward().
TEST(AlignedReadStorageTest, ResolveCollapsesReadsOfBothOrientations) {
    ag::AssemblyGraph graph(2);
    ag::Vertex &c = graph.addVertex(Sequence("ACGG"));
    ag::Vertex &u1 = graph.addVertex(Sequence("TACGG"));
    ag::Vertex &u2 = graph.addVertex(Sequence("GACGG"));
    ag::Vertex &w1 = graph.addVertex(Sequence("ACGGA"));
    ag::Vertex &w2 = graph.addVertex(Sequence("ACGGT"));
    ag::Edge &in1 = graph.addSPEdge(u1, c);
    ag::Edge &in2 = graph.addSPEdge(u2, c);
    ag::Edge &out1 = graph.addSPEdge(c, w1);
    ag::Edge &out2 = graph.addSPEdge(c, w2);
    ASSERT_TRUE(in1.isSuffix() && in2.isSuffix() && out1.isPrefix() && out2.isPrefix());
    for (ag::Vertex *v : {&c, &u1, &u2, &w1, &w2})
        AddSample(*v, 0, v->size());
    auto make = [](std::initializer_list<ag::Edge *> edges) {
        ag::GraphPath res(**edges.begin());
        for (auto it = edges.begin() + 1; it != edges.end(); ++it)
            res += **it;
        return res;
    };
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("through1", make({&in1, &out1}));
    reads.emplace_back("through2", make({&in2, &out2}));
    reads.emplace_back("forward", make({&out1}));
    reads.emplace_back("backward", make({&out2.rc()}));
    ag::AlignedReadStorage storage(graph, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency(std::cerr));
    std::map<std::string, std::string> before = SeqsById(storage);

    ag::VertexResolutionPlan plan(c);
    plan.add(in1, out1);
    plan.add(in2, out2);
    graph.resolveVertex(c, plan);
    graph.removeMarked();

    EXPECT_EQ(SeqsById(storage), before);
    for (const std::string id : {"forward", "backward"}) {
        const ag::GraphPath &path = ReadById(storage, id).getPath();
        EXPECT_TRUE(path.valid() && !path.isLegacy() && path.empty()) << id;
    }
    EXPECT_TRUE(storage.checkConsistency(std::cerr));
    EXPECT_TRUE(graph.checkConsistency());
}
