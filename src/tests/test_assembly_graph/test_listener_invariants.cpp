#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "graph_test_utils.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "dbg/dbg_graph_aligner.hpp"

using namespace ag_test;

namespace {
    constexpr size_t K = 5;
    constexpr size_t GENOME_LEN = 1000;

    // Runs every ag-layer listener's invariant against the current graph/storage/tracker:
    //  read sequences unchanged, AlignedReadStorageMaintenance buckets consistent, SuffixTracker
    //  content equals read-derived content AND equals a from-scratch refill, edge codes navigable.
    void CheckAllInvariants(ag::AssemblyGraph &graph, ag::AlignedReadStorage &storage,
                            ag::SuffixTracker &tracker,
                            const std::map<std::string, std::string> &before,
                            const std::string &ctx) {
        std::string err;
        EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << ctx << "\n" << err;
        EXPECT_TRUE(storage.checkConsistency()) << ctx;
        EXPECT_TRUE(tracker.fireCheckConsistency()) << ctx;
        EXPECT_TRUE(EdgeCodesValid(graph, err)) << ctx << "\n" << err;
        ag::SuffixTracker fresh(storage, graph, 0, 1000);
        fresh.fillFromStorage(QuietLogger(), 1);
        EXPECT_TRUE(SuffixContentEqual(tracker, fresh, graph, err)) << ctx << "\n" << err;
    }
}

// Smoke test: a random genome builds a consistent graph, reads align to non-empty paths, and a
// SuffixTracker fills without crashing. Validates the construction harness before the real invariant
// tests build on it.
TEST(ListenerInvariants, BuildRandomGraphSmoke) {
    RandomDBG rdbg(RandomGenome(GENOME_LEN, 1), K);
    std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 100, 2, 2);

    size_t non_empty = 0;
    for (const ag::AlignedRead &r : reads)
        if (!r.getPath().empty()) non_empty++;
    EXPECT_GT(non_empty, 0u);

    ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));
    EXPECT_TRUE(storage.checkConsistency());

    ag::SuffixTracker tracker(storage, rdbg.dbg, 0, 1000);
    tracker.fillFromStorage(QuietLogger(), 1);
    EXPECT_TRUE(storage.checkConsistency());
}

// Core random test: condense the full de Bruijn graph with MergeAllToEdges and verify that every
// graph-operation listener still stores exactly what it should.
//   - read sequences are unchanged (the central invariant);
//   - AlignedReadStorageMaintenance: storage.checkConsistency() (every read filed under the right
//     start/finish-rc vertex buckets);
//   - SuffixTracker: fireCheckConsistency() (subtracting all read subpaths empties every record) and
//     a direct refill-equality comparison against a tracker rebuilt from scratch on the final graph;
//   - EdgeCodeListener: outgoing edges remain distinguishable by their code's first nucleotide.
TEST(ListenerInvariants, RandomMergeAllPreservesInvariants) {
    for (uint32_t seed = 1; seed <= 5; seed++) {
        RandomDBG rdbg(RandomGenome(GENOME_LEN, seed), K);
        std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 60, 120, 2, seed + 100);

        ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));
        ag::SuffixTracker tracker(storage, rdbg.dbg, 0, 1000);
        tracker.fillFromStorage(QuietLogger(), 1);
        ASSERT_TRUE(storage.checkConsistency()) << "seed=" << seed;

        std::map<std::string, std::string> before = ReadSeqSnapshot(storage);

        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);

        std::string err;
        EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << "seed=" << seed << "\n" << err;
        EXPECT_TRUE(storage.checkConsistency()) << "seed=" << seed;
        EXPECT_TRUE(tracker.fireCheckConsistency()) << "seed=" << seed;
        EXPECT_TRUE(EdgeCodesValid(rdbg.dbg, err)) << "seed=" << seed << "\n" << err;

        ag::SuffixTracker fresh(storage, rdbg.dbg, 0, 1000);
        fresh.fillFromStorage(QuietLogger(), 1);
        EXPECT_TRUE(SuffixContentEqual(tracker, fresh, rdbg.dbg, err)) << "seed=" << seed << "\n" << err;
    }
}

// ---- Tiny hand-built examples (extra-small, exact) ----------------------------------------------

// A read that ends strictly inside the run being merged. mergePathToEdge collapses the interior
// edges 1..3 (v1->v4) of a 5-edge chain while read "short" (edges 0,1,2) ends at v3, an inner vertex
// that the merge deletes - exercising the read-ends-inside-merged-run path through every listener.
TEST(ListenerInvariants, TinyMergePathToEdge) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(6); // v0..v5, edges 0..4
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("short", MakePath(chain, 0, 3)); // v0->v3, ends inside the merged run
    reads.emplace_back("long", MakePath(chain, 0, 5));  // v0->v5, spans everything
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ag::SuffixTracker tracker(storage, gb.graph, 0, 1000);
    tracker.fillFromStorage(TestLogger(), 1);
    ASSERT_TRUE(storage.checkConsistency());

    std::map<std::string, std::string> before = ReadSeqSnapshot(storage);

    gb.graph.mergePathToEdge(MakePath(chain, 1, 4)); // merge interior run v1->v4
    gb.graph.removeMarked();

    CheckAllInvariants(gb.graph, storage, tracker, before, "after mergePathToEdge");
}

// splitEdge respecting its no-cut-reads contract: the read ends inside what becomes the first piece,
// so the split point never falls within a read span. Checks the full listener invariant set (the
// existing split tests only check storage consistency).
TEST(ListenerInvariants, TinySplitEdge) {
    GraphBuilder gb(8);
    std::vector<ag::EdgeId> chain = gb.buildChain(3); // e0 (before), e1 (to be split)
    ag::GraphPath path(*chain[0]);
    path += *chain[1];
    path.setCutRight(5); // ends 3 bases into e1, inside what becomes e1's first piece [0,4)
    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r0", path);
    ag::AlignedReadStorage storage(gb.graph, std::move(reads));
    ag::SuffixTracker tracker(storage, gb.graph, 0, 1000);
    tracker.fillFromStorage(TestLogger(), 1);
    ASSERT_TRUE(storage.checkConsistency());

    std::map<std::string, std::string> before = ReadSeqSnapshot(storage);

    gb.graph.splitEdge(*chain[1], {ag::EdgePosition(*chain[1], 4)});

    CheckAllInvariants(gb.graph, storage, tracker, before, "after splitEdge");
}

// ---- Vertex resolution (supregraph) -------------------------------------------------------------

// resolveVertex on a PROPER 2-in / 2-out de Bruijn junction (built from two designed sequences that
// share exactly the middle 5-mer ACGTA; every edge is a real (k+1)-mer step, so the splitting
// property holds in both directions). The core's four incident edges are corporeal (non suffix/prefix)
// single-nucleotide DBG edges - the case the AB_dev refactoring must support ("supregraph may have
// non-suffix, non-prefix edges; edgeToSupreVertex must run freely"). ExtendAdjacent turns each into a
// supre-vertex + suffix/prefix edge and resolveVertex splits the core along the one-to-one plan each
// read traverses.
//
// Checks the full listener invariant set: read sequences unchanged, AlignedReadStorageMaintenance
// buckets consistent, and the graph structurally consistent (Vertex::checkConsistency).
TEST(ListenerInvariants, ResolveVertexSupregraphConsistent) {
    SeqDBG g({"CACGTAC", "GACGTAG"}, 5); // shared 5-mer ACGTA is the only junction
    ag::GraphPath p1 = g.align("CACGTAC"); // CACGT -> ACGTA -> CGTAC
    ag::GraphPath p2 = g.align("GACGTAG"); // GACGT -> ACGTA -> CGTAG
    ASSERT_EQ(EdgeCount(p1), 2u);
    ASSERT_EQ(EdgeCount(p2), 2u);
    ag::Vertex &core = p1.frontEdge().getFinish();
    ASSERT_EQ(&core, &p2.frontEdge().getFinish());
    ASSERT_TRUE(core.isCore());
    ASSERT_EQ(core.inDeg(), 2u);
    ASSERT_EQ(core.outDeg(), 2u);

    std::vector<ag::AlignedRead> reads;
    reads.emplace_back("r1", p1);
    reads.emplace_back("r2", p2);
    ag::AlignedReadStorage storage(g.dbg, std::move(reads));
    ASSERT_TRUE(storage.checkConsistency());
    std::map<std::string, std::string> before = ReadSeqSnapshot(storage);

    // Pair each read's incoming edge with its outgoing edge (a clean one-to-one matching).
    // Precondition asserts: the stored read paths must still run through `core` exactly as p1/p2 did
    // (if a build sees these fire instead of the plan.add below, the read paths were altered - e.g.
    // normalized/invalidated by storage construction - rather than plan.add itself being at fault).
    for (size_t i = 0; i < 2; i++) {
        ASSERT_TRUE(storage[i].valid() && !storage[i].getPath().empty()) << "read " << i;
        ASSERT_EQ(EdgeCount(storage[i].getPath()), 2u) << "read " << i;
        ASSERT_EQ(&storage[i].getPath().frontEdge().getFinish(), &core) << "read " << i;
        ASSERT_EQ(&storage[i].getPath().backEdge().getStart(), &core) << "read " << i;
    }
    ag::VertexResolutionPlan plan(core);
    plan.add(storage[0].getPath().frontEdge(), storage[0].getPath().backEdge());
    plan.add(storage[1].getPath().frontEdge(), storage[1].getPath().backEdge());
    g.dbg.resolveVertex(core, plan);
    g.dbg.removeMarked();

    for (ag::Vertex &v : g.dbg.vertices())
        v.checkConsistency();
    std::string err;
    EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << "after resolveVertex\n" << err;
    EXPECT_TRUE(storage.checkConsistency()) << "after resolveVertex";
}

// Build a random DBG, condense it, then convert it into a real supregraph via edgeToSupreVertex on
// every corporeal edge (exactly Supregraph::FromDBG). Verifies the ag-layer listeners that operate on
// the resulting supregraph WITHOUT re-walking read paths:
//   - read sequences unchanged (GraphPath::Seq, which does walk - and handles the new suffix edges);
//   - AlignedReadStorageMaintenance read-path buckets consistent;
//   - EdgeCodeListener codes navigable.
// SuffixTracker is intentionally NOT attached here - see DISABLED_SuffixTrackerOnSupregraph below for
// why (its processPath cannot yet traverse the suffix edges the conversion introduces).
TEST(ListenerInvariants, DBGToSupregraphConversionPreservesInvariants) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(400, seed), K);
        std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 90, 2, seed + 200);

        ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));
        ASSERT_TRUE(storage.checkConsistency()) << "seed=" << seed;

        std::map<std::string, std::string> before = ReadSeqSnapshot(storage);

        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg); // condense, like a real DBG
        ConvertToSupregraph(rdbg.dbg);                   // DBG -> supregraph

        std::string ctx = "after DBG->SPG conversion (seed=" + std::to_string(seed) + ")";
        std::string err;
        EXPECT_TRUE(ReadSeqsUnchanged(before, storage, err)) << ctx << "\n" << err;
        EXPECT_TRUE(storage.checkConsistency()) << ctx;
        EXPECT_TRUE(EdgeCodesValid(rdbg.dbg, err)) << ctx << "\n" << err;
    }
}

// SuffixTracker built on an already-converted supregraph: fillFromStorage walks read paths through the
// conversion's suffix edges via processPath, and fireCheckConsistency removes them again. With the
// trailing-suffix-edge handling in processPath, add and remove are symmetric and consistency holds -
// i.e. SuffixTracker::processPath now works on a supregraph (regression test for that fix).
TEST(ListenerInvariants, SuffixTrackerFillsOnSupregraph) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(400, seed), K);
        std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 90, 2, seed + 200);
        ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));

        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        ConvertToSupregraph(rdbg.dbg);

        // Tracker built on the already-converted supregraph: add (fill) and remove (check) are both
        // processPath-on-supregraph, so this exercises processPath's suffix-edge handling end to end.
        ag::SuffixTracker tracker(storage, rdbg.dbg, 0, 1000);
        tracker.fillFromStorage(QuietLogger(), 1);
        EXPECT_TRUE(tracker.fireCheckConsistency()) << "seed=" << seed;
        EXPECT_TRUE(storage.checkConsistency()) << "seed=" << seed;
    }
}

// Carrying a live SuffixTracker THROUGH the DBG->SPG conversion must land on the same records as a
// from-scratch fill on the converted supregraph. The tracker is filled on the condensed DBG, then
// edgeToSupreVertex runs on every corporeal edge with the tracker attached - each conversion routes
// through SuffixTracker::fireEdgeToSupreVertex, which copies the edge's suffixes onto both boundary
// edges of the new supre-vertex. A second tracker filled from scratch on the final supregraph is the
// reference; the two must store exactly the same suffixes on every edge, and the carried-through
// tracker must be self-consistent against the reads.
TEST(ListenerInvariants, SuffixTrackerThroughConversionMatchesFresh) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(400, seed), K);
        std::vector<ag::AlignedRead> reads = ExtractAndAlignReads(rdbg, 40, 90, 2, seed + 200);
        ag::AlignedReadStorage storage(rdbg.dbg, std::move(reads));

        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg); // condense to a realistic DBG
        ag::SuffixTracker carried(storage, rdbg.dbg, 0, 1000);
        carried.fillFromStorage(QuietLogger(), 1);       // fill BEFORE conversion

        ConvertToSupregraph(rdbg.dbg);                   // carries `carried` through via fireEdgeToSupreVertex

        // Reference tracker filled from scratch on the already-converted supregraph.
        ag::SuffixTracker fresh(storage, rdbg.dbg, 0, 1000);
        fresh.fillFromStorage(QuietLogger(), 1);

        std::string err;
        EXPECT_TRUE(SuffixContentEqual(carried, fresh, rdbg.dbg, err)) << "seed=" << seed << "\n" << err;
        EXPECT_TRUE(carried.fireCheckConsistency()) << "carried self-consistent, seed=" << seed;
    }
}

// DLLAlignmentStorage ("honest path tracking") is designed to keep a contig's alignment valid even as
// parts of the path are deleted. Track the random genome's own path through the graph, then subject the
// graph to merges AND edge deletions, and check after each that the tracked alignment stays valid (no
// fragment ever references a deleted edge/vertex; the path just loses the deleted stretches).
TEST(ListenerInvariants, DLLAlignmentToleratesDeletions) {
    RandomDBG rdbg(RandomGenome(400, 1), K);
    ag::DLLAlignmentStorage dll(rdbg.dbg);

    // Align the genome to its own graph and register that path with the tracker (as PrepareDLLPathTracker
    // does in the real pipeline: carefulAlign -> addContig), while the graph is still a hashed DBG.
    dbg::KmerIndex index(rdbg.dbg);
    index.fillAnchors(QuietLogger(), 1, rdbg.dbg, 2);
    Contig genome(Sequence(rdbg.genome), "genome");
    auto al = index.carefulAlign(genome);
    ASSERT_FALSE(al.empty());
    dll.addContig(genome, al);
    ASSERT_TRUE(dll.checkConsistency()) << "initial";

    // Merge operation: unbranching runs collapse; the tracked path must follow onto the merged edges.
    ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
    EXPECT_TRUE(dll.checkConsistency()) << "after MergeAllToEdges";

    // Delete a selection of edges the path may run through. Each deletion fires fireDeleteEdge, which
    // must drop exactly the fragments on that edge and leave the rest of the path valid.
    std::vector<ag::EdgeId> eids;
    for (ag::Edge &e : rdbg.dbg.edgesUnique())
        eids.push_back(e.getId());
    size_t deleted = 0;
    for (size_t i = 0; i < eids.size(); i += 3) { // every third edge
        ag::EdgeId eid = eids[i];
        if (!eid->getStart().isCore() || !eid->getFinish().isCore())
            continue; // keep it simple: only delete between non-degenerate vertices
        rdbg.dbg.removeEdge(*eid);
        deleted++;
        ASSERT_TRUE(dll.checkConsistency()) << "after deleting edge " << i;
    }
    EXPECT_GT(deleted, 0u);
}

// mergePath collapses an interior sub-path into a single vertex. Four pre-existing bugs, all fixed:
//   (1) Locker deadlock (mergePathLockFree uses addEdgeLockFree, assembly_graph.cpp:94,102);
//   (2) EdgeCodeListener/SuffixTracker::fireMergePath look up the new (still non-corporeal) boundary
//       edges with getNonCorporealOutgoing instead of getOutgoing;
//   (3) mergePathLockFree isolates the path's vertices from a snapshot taken before the edges were
//       added, instead of re-walking path.vertices() over the mutated graph;
//   (4) mergePathLockFree's outgoing block iterated path.getStart() instead of path.getFinish(),
//       so res was re-attached to an interior path vertex (res -> B) instead of the finish's
//       successor; isolateAndMark(B) then deleted that non-corporeal edge (violating the "only
//       corporeal edges are deleted" pipeline invariant), which drifted _outDeg and SEGVd, and the
//       real outgoing edge was never created so the graph was severed (merged vertex outDeg 0).
//       Fixed by iterating path.getFinish().
TEST(ListenerInvariants, MergePathOnInteriorSubpath) {
    SeqDBG g({RandomGenome(28, 1)}, K);
    std::vector<ag::GraphPath> paths = ag::AllUnbranchingPaths(QuietLogger(), 1, g.dbg);
    ag::GraphPath chosen;
    for (ag::GraphPath &p : paths)
        if (EdgeCount(p) >= 4) { chosen = p; break; }
    ASSERT_TRUE(chosen.valid());
    std::vector<ag::Edge *> edges;
    for (ag::Edge &e : chosen.edges())
        edges.push_back(&e);
    ag::GraphPath sub(*edges[1]); // interior 2-edge sub-path
    sub += *edges[2];
    std::string subSeq = sub.Seq().str();
    // The merged vertex absorbs the start's predecessors and the finish's successors, so its
    // in/out degrees must equal those of the first/last vertices of the merged path.
    size_t expIn = sub.getStart().inDeg();
    size_t expOut = sub.getFinish().outDeg();

    ag::Vertex &U = g.dbg.mergePath(sub);
    g.dbg.removeMarked();

    EXPECT_EQ(U.getSeq().str(), subSeq) << "merged vertex spells the merged sub-path";
    EXPECT_EQ(U.inDeg(), expIn) << "merged vertex inDeg equals first path vertex inDeg";
    EXPECT_EQ(U.outDeg(), expOut) << "merged vertex outDeg equals last path vertex outDeg";
    EXPECT_GT(U.outDeg(), 0u) << "merge must not sever the graph";
    for (ag::Vertex &v : g.dbg.vertices())
        v.checkConsistency();
}
