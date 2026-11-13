// Read correction (rerouting and invalidating reads): the graph/read mode switch, every read listener after
// corrections, and the DBG cleaning pipeline that invalidates most reads and then removes uncovered parts of edges.

#include "gtest/gtest.h"
#include "random_graph_utils.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
#include "dbg/graph_modification.hpp"
#include "error_correction/read_cleaning.hpp"
#include "supregraph/supregraph.hpp"
#include "assembly_graph/data_structures/read_logger.hpp"
#include <fstream>

using namespace ag_test;

namespace {
    constexpr size_t K = 5;

//    A random DBG with aligned reads in a storage that tracks coverage and suffixes.
    struct CorrectionFixture {
        RandomDBG rdbg;
        dbg::DBGAlignedReadStorage storage;

        CorrectionFixture(std::string genome, size_t k, size_t read_cnt, size_t read_len, uint32_t seed) :
                rdbg(std::move(genome), k),
                storage(QuietLogger(), 1, rdbg.dbg, ExtractAndAlignReads(rdbg, read_cnt, read_len, 1, seed + 50), true) {
            storage.trackSuffixes(QuietLogger(), 1, rdbg.dbg, 0, 1000);
        }
        explicit CorrectionFixture(uint32_t seed) : CorrectionFixture(RandomGenome(600, seed), K, 150, 40, seed) {}
    };

//    The storage with every read listener, the graph with every graph listener, and the suffix tracker against one
//    rebuilt from the current reads.
    void CheckAll(ag::AssemblyGraph &graph, ag::AlignedReadStorage &reads, ag::SuffixTracker &suffixes,
                  const std::string &ctx) {
        std::stringstream rs, gs;
        EXPECT_TRUE(reads.checkConsistency(rs)) << ctx << "\n" << rs.str();
        EXPECT_TRUE(graph.checkConsistency(gs)) << ctx << "\n" << gs.str();
        ag::SuffixTracker fresh(reads, graph, 0, 1000);
        fresh.fillFromStorage(QuietLogger(), 1);
        std::string err;
        EXPECT_TRUE(SuffixContentEqual(suffixes, fresh, graph, err)) << ctx << "\n" << err;
    }

    size_t ValidReads(ag::AlignedReadStorage &reads) {
        size_t res = 0;
        for (ag::AlignedRead &read : reads)
            if (read.valid())
                res++;
        return res;
    }
}

// ---- The mode switch ------------------------------------------------------------------------------------------------

TEST(ReadCorrectionDeathTest, ReadOperationsNeedReadMode) {
    CorrectionFixture f(1);
    ag::AlignedReadStorage &reads = f.storage.getReads();
    ag::AlignedRead &read = reads[0];
    EXPECT_DEATH(reads.delayedInvalidateRead(read, "test"), "");
    EXPECT_DEATH(reads.rerouteRead(read, read.getPath(), "test"), "");
    EXPECT_DEATH(reads.applyCorrections(QuietLogger(), 1), "");
    reads.startReadCorrection();
    EXPECT_DEATH(reads.startReadCorrection(), "");
}

TEST(ReadCorrectionDeathTest, GraphOperationsNeedGraphMode) {
    CorrectionFixture f(1);
    ag::AlignedReadStorage &reads = f.storage.getReads();
    reads.startReadCorrection();
    EXPECT_TRUE(reads.inReadCorrection());
    EXPECT_DEATH(ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg), "inReadCorrection");
    reads.applyCorrections(QuietLogger(), 1);
    EXPECT_FALSE(reads.inReadCorrection());
//    Back in graph mode: graph operations work again.
    ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
    CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), "after returning to graph mode");
}

// ---- Rerouting on a DBG ---------------------------------------------------------------------------------------------

// Reads are trimmed to inner subpaths, moved onto their reverse complement (another start vertex) or invalidated;
// every listener matches the reads after the corrections, and graph operations afterwards keep it so.
TEST(ReadCorrection, RerouteAndInvalidateOnDBG) {
    for (uint32_t seed = 1; seed <= 3; seed++) {
        CorrectionFixture f(seed);
        ag::AlignedReadStorage &reads = f.storage.getReads();
        ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), "seed=" + std::to_string(seed) + " before");
        size_t trimmed = 0, moved = 0, invalidated = 0;
        reads.startReadCorrection();
        for (size_t i = 0; i < reads.size(); i++) {
            ag::AlignedRead &read = reads[i];
            if (!read.valid() || read.getPath().empty())
                continue;
            const ag::GraphPath &path = read.getPath();
            if (i % 4 == 0) {
                reads.delayedInvalidateRead(read, "test invalidate");
                invalidated++;
            } else if (i % 4 == 1 && EdgeCount(path) >= 3) {
                reads.rerouteRead(read, path.subPath(path.firstPosition() + 1, path.lastPosition() - 1), "test trim");
                trimmed++;
            } else if (i % 4 == 2) {
                reads.rerouteRead(read, path.RC(), "test rc");
                moved++;
            }
        }
        reads.applyCorrections(QuietLogger(), 1);
        EXPECT_GT(trimmed, 0u);
        EXPECT_GT(moved, 0u);
        EXPECT_GT(invalidated, 0u);
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), "seed=" + std::to_string(seed) + " after corrections");
        std::map<std::string, std::string> before = ReadSeqSnapshot(reads);
        ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), "seed=" + std::to_string(seed) + " after merging again");
        std::string err;
        EXPECT_TRUE(ReadSeqsUnchanged(before, reads, err)) << "seed=" << seed << " " << err;
    }
}

// Several correction rounds in a row, each applied before the next starts.
TEST(ReadCorrection, RepeatedRounds) {
    CorrectionFixture f(4);
    ag::AlignedReadStorage &reads = f.storage.getReads();
    ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
    for (size_t round = 0; round < 3; round++) {
        reads.startReadCorrection();
        for (size_t i = round; i < reads.size(); i += 3) {
            ag::AlignedRead &read = reads[i];
            if (read.valid() && !read.getPath().empty())
                reads.rerouteRead(read, read.getPath().RC(), "round");
        }
        reads.applyCorrections(QuietLogger(), 1);
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), "round " + std::to_string(round));
    }
}

// ---- Subreads on a supregraph ---------------------------------------------------------------------------------------

// Reads become subreads of their start vertex (path -> substring of one vertex) and then go back to their original
// paths (substring -> path): both directions of the start-record move.
TEST(ReadCorrection, PathToSubreadAndBackOnSupregraph) {
    size_t to_subread = 0;
    for (uint32_t seed = 1; seed <= 3; seed++) {
        RandomDBG rdbg(RandomGenome(300, seed), K);
        std::unique_ptr<spg::Supregraph> spg_holder;
        ag::AlignedReadStorage reads(rdbg.dbg, ExtractAndAlignReads(rdbg, 60, 40, 1, seed + 60));
        ag::SuffixTracker suffixes(reads, rdbg.dbg, 0, 1000);
        suffixes.fillFromStorage(QuietLogger(), 1);
        ag::MergeAllToEdges(QuietLogger(), 1, rdbg.dbg);
        rdbg.dbg.disableHashing();
        spg_holder.reset(new spg::Supregraph(spg::Supregraph::FromDBG(QuietLogger(), 1, std::move(rdbg.dbg))));
        spg::Supregraph &spg = *spg_holder;
        CheckAll(spg, reads, suffixes, "seed=" + std::to_string(seed) + " after conversion");
        std::vector<ag::GraphPath> original;
        for (ag::AlignedRead &read : reads)
            original.push_back(read.getPath());
        reads.startReadCorrection();
        for (ag::AlignedRead &read : reads)
            if (read.valid() && !read.getPath().empty()) {
                ag::Vertex &v = read.getPath().getStart();
                reads.rerouteRead(read, ag::GraphPath(v, read.getPath().leftCut(), 0), "to subread");
            }
        reads.applyCorrections(QuietLogger(), 1);
        for (ag::AlignedRead &read : reads)
            if (read.valid() && read.getPath().empty())
                to_subread++;
        CheckAll(spg, reads, suffixes, "seed=" + std::to_string(seed) + " after making subreads");
        reads.startReadCorrection();
        for (size_t i = 0; i < reads.size(); i++)
            if (reads[i].valid() && original[i].valid())
                reads.rerouteRead(reads[i], original[i], "back");
        reads.applyCorrections(QuietLogger(), 1);
        CheckAll(spg, reads, suffixes, "seed=" + std::to_string(seed) + " after restoring paths");
    }
    EXPECT_GT(to_subread, 0u);
}

// ---- The DBG cleaning pipeline --------------------------------------------------------------------------------------

namespace {
//    A random genome with a few repeat units, each written over several random places: at a large k its DBG has some
//    junctions (it is not a single path) but long edges between them.
    std::string RepeatGenome(size_t len, uint32_t seed, size_t units = 16, size_t unit_len = 150, size_t copies = 3) {
        std::string genome = RandomGenome(len, seed);
        std::mt19937 rng(seed + 1000);
        std::uniform_int_distribution<size_t> pos(0, len - unit_len);
        for (size_t u = 0; u < units; u++) {
            std::string unit = RandomGenome(unit_len, seed * 100 + u);
            for (size_t c = 0; c < copies; c++)
                genome.replace(pos(rng), unit_len, unit);
        }
        return genome;
    }

//    Counts the split points of every splitEdge (each split fires for both strands).
    struct SplitCounter : public ag::ResolutionListener {
        size_t points = 0;
        explicit SplitCounter(ag::AssemblyGraph &graph) : ag::ResolutionListener(graph, "SplitCounter") {}
        void fireSplitEdge(ag::Edge &, const ag::RAGraphPath &split) override {points += split.size() - 1;}
    };

    size_t TotalEdgeLength(ag::AssemblyGraph &graph) {
        size_t res = 0;
        for (ag::Edge &e : graph.edgesUnique())
            res += e.truncSize();
        return res;
    }
}

// A 15 kb genome with repeats at k = 15 (about 50 junctions) and 200-long reads at about 10x. Reads are invalidated in random order until the
// remaining ones total about 1.5 genome lengths, so coverage has many gaps; RemoveUncovered then splits edges at the
// borders of what the remaining reads cover, removes the uncovered parts and condenses the graph. Every listener stays
// consistent and the remaining reads still spell the same sequences.
TEST(ReadCorrection, InvalidateMostThenRemoveUncovered) {
    const size_t genome_len = 15000;
    for (uint32_t seed = 1; seed <= 3; seed++) {
        const std::string ctx = "seed=" + std::to_string(seed);
        CorrectionFixture f(RepeatGenome(genome_len, seed), 15, 750, 200, seed);
        ag::AlignedReadStorage &reads = f.storage.getReads();
        ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
        size_t junctions = 0;
        for (ag::Vertex &v : f.rdbg.dbg.verticesUnique())
            if (v.isJunction())
                junctions++;
        EXPECT_GT(junctions, 4u) << ctx << ": the graph must not be a single path";
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), ctx + " before");

        size_t remaining = 0;
        for (ag::AlignedRead &read : reads)
            if (read.valid())
                remaining += read.getPath().Seq().size();
        std::vector<size_t> order(reads.size());
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), std::mt19937(seed));
        reads.startReadCorrection();
        for (size_t i : order)
            if (remaining * 2 > genome_len * 3 && reads[i].valid()) {
                remaining -= reads[i].getPath().Seq().size();
                reads.delayedInvalidateRead(reads[i], "test");
            }
        reads.applyCorrections(QuietLogger(), 1);
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), ctx + " after invalidation");

        std::map<std::string, std::string> before = ReadSeqSnapshot(reads);
        size_t length_before = TotalEdgeLength(f.rdbg.dbg);
        SplitCounter splits(f.rdbg.dbg);
        dbg::RemoveUncovered(QuietLogger(), 1, f.rdbg.dbg, {&reads});
        EXPECT_GT(splits.points / 2, 20u) << ctx << ": RemoveUncovered must split edges at many points";
        EXPECT_LT(TotalEdgeLength(f.rdbg.dbg), length_before) << ctx << ": uncovered parts must be removed";
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), ctx + " after RemoveUncovered");
        std::string err;
        EXPECT_TRUE(ReadSeqsUnchanged(before, reads, err)) << ctx << " " << err;
//        Every remaining edge between two junctions is covered by a remaining read.
        for (ag::Edge &e : f.rdbg.dbg.edges())
            if (e.getStart().isJunction() && e.getFinish().isJunction())
                EXPECT_GT(e.intCov(), 0u) << ctx << " edge " << e.getId();
    }
}

// ---- ReadLogger -----------------------------------------------------------------------------------------------------

namespace {
    std::string ReadFile(const std::experimental::filesystem::path &path) {
        std::ifstream is(path);
        std::stringstream ss;
        ss << is.rdbuf();
        return ss.str();
    }
}

// Every rerouted read gets an initial/corrected record and every invalidated one an "invalidated" record (an
// invalidation used to be logged as a rerouting to the invalid path).
TEST(ReadLogger, LogsReroutesAndInvalidations) {
    CorrectionFixture f(2);
    ag::AlignedReadStorage &reads = f.storage.getReads();
    ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
    std::experimental::filesystem::path log_file = std::experimental::filesystem::temp_directory_path() / "lja_read_logger_test.txt";
    std::vector<std::string> rerouted, invalidated;
    {
        ag::ReadLogger logger(reads, 1, log_file);
        reads.startReadCorrection();
        for (size_t i = 0; i < reads.size(); i++) {
            ag::AlignedRead &read = reads[i];
            if (!read.valid() || read.getPath().empty())
                continue;
            if (i % 3 == 0) {
                reads.delayedInvalidateRead(read, "test_invalidate");
                invalidated.push_back(read.getId());
            } else if (i % 3 == 1) {
                reads.rerouteRead(read, read.getPath().RC(), "test_reroute");
                rerouted.push_back(read.getId());
            }
        }
        reads.applyCorrections(QuietLogger(), 1);
    }
    std::string log = ReadFile(log_file);
    std::experimental::filesystem::remove(log_file);
    ASSERT_FALSE(rerouted.empty());
    ASSERT_FALSE(invalidated.empty());
    for (const std::string &id : rerouted) {
        EXPECT_NE(log.find(id + " test_reroute "), std::string::npos) << id;
        EXPECT_NE(log.find(id + "  initial  "), std::string::npos) << id;
        EXPECT_NE(log.find(id + " corrected "), std::string::npos) << id;
    }
    for (const std::string &id : invalidated) {
        EXPECT_NE(log.find(id + " invalidated test_invalidate\n"), std::string::npos) << id;
        EXPECT_NE(log.find(id + "    final    "), std::string::npos) << id;
    }
}

namespace {
//    Invalidates the first two reads from two threads, one each.
    void InvalidateFromTwoThreads(ag::AlignedReadStorage &reads) {
        reads.startReadCorrection();
#pragma omp parallel for num_threads(2) schedule(static, 1)
        for (size_t i = 0; i < 2; i++)
            reads.delayedInvalidateRead(reads[i], "test");
    }
}

TEST(ReadLoggerDeathTest, MoreThreadsThanBuilt) {
    CorrectionFixture f(2);
    ag::AlignedReadStorage &reads = f.storage.getReads();
    std::experimental::filesystem::path log_file = std::experimental::filesystem::temp_directory_path() / "lja_read_logger_death.txt";
    EXPECT_DEATH({
        ag::ReadLogger logger(reads, 1, log_file);
        InvalidateFromTwoThreads(reads);
    }, "logs.size");
    std::experimental::filesystem::remove(log_file);
}

// Rerouting reads away and back leaves zero-count suffix entries behind; records clean them up once they make up a
// third of the record (they used to pile up forever, as the cleanup was unreachable).
TEST(ReadCorrection, SuffixRecordsDropZeroEntries) {
    CorrectionFixture f(5);
    ag::AlignedReadStorage &reads = f.storage.getReads();
    ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
    for (size_t round = 0; round < 6; round++) {
        reads.startReadCorrection();
        for (ag::AlignedRead &read : reads)
            if (read.valid() && !read.getPath().empty())
                reads.rerouteRead(read, read.getPath().RC(), "round");
        reads.applyCorrections(QuietLogger(), 1);
        for (ag::Edge &e : f.rdbg.dbg.edges())
            if (f.storage.getSuffixes().hasSuffixRecord(e)) {
                const ag::SuffixRecord &rec = f.storage.getSuffixes().getSuffixRecord(e);
                size_t zeros = 0;
                for (const std::pair<Sequence, int> &p : rec)
                    if (p.second == 0)
                        zeros++;
                EXPECT_LE(zeros, rec.getSuffixes().size() / 3 + 1) << "round " << round << " edge " << e.getId();
            }
        CheckAll(f.rdbg.dbg, reads, f.storage.getSuffixes(), "round " + std::to_string(round));
    }
}

// A moved DBGAlignedReadStorage takes its trackers along: they listen to and read the new storage, and only the new
// storage deletes them (the defaulted move used to leave both owning the same trackers).
TEST(ReadCorrection, MovedStorageKeepsItsTrackers) {
    CorrectionFixture f(6);
    dbg::DBGAlignedReadStorage moved(std::move(f.storage));
    ag::AlignedReadStorage &reads = moved.getReads();
    ag::MergeAllToEdges(QuietLogger(), 1, f.rdbg.dbg);
    CheckAll(f.rdbg.dbg, reads, moved.getSuffixes(), "after merging");
    reads.startReadCorrection();
    for (size_t i = 0; i < reads.size(); i += 4)
        if (reads[i].valid() && !reads[i].getPath().empty())
            reads.delayedInvalidateRead(reads[i], "test");
    reads.applyCorrections(QuietLogger(), 1);
    CheckAll(f.rdbg.dbg, reads, moved.getSuffixes(), "after corrections");
    EXPECT_FALSE(f.storage.tracksSuffixes());
}
