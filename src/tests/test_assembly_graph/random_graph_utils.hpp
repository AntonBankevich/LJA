#pragma once

// Helpers for building random assembly graphs from a random genome and extracting aligned reads,
// for listener-invariant tests. A random genome of length ~1000 with small k (~5) produces a graph
// with many repeat-induced branches - a deliberately challenging exercise for the graph-modification
// listeners (SuffixTracker, AlignedReadStorageMaintenance, EdgeCodeListener, ...).
//
// Construction builds the *full* de Bruijn graph: every distinct genomic k-mer becomes a vertex and
// every consecutive k-mer pair a single-nucleotide edge. MergeAllToEdges then condenses unbranching
// runs into real multi-nucleotide edges, so the operations under test (merge/split) act on a
// realistic graph. Odd k (e.g. 5) is used so no k-mer can equal its own reverse complement (no
// palindromic / self-rc vertices).

#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/ag_algorithms.hpp"
#include "assembly_graph/data_structures/read_alignment_storage.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_algorithms.hpp"
#include "dbg/dbg_graph_aligner.hpp"
#include "common/logging.hpp"
#include <random>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>

namespace ag_test {

    inline logging::Logger &QuietLogger() {
        static logging::Logger logger(false);
        return logger;
    }

    // Deterministic random genome over {A,C,G,T}.
    inline std::string RandomGenome(size_t len, uint32_t seed) {
        static const char bases[] = {'A', 'C', 'G', 'T'};
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> dist(0, 3);
        std::string g(len, 'A');
        for (size_t i = 0; i < len; i++)
            g[i] = bases[dist(rng)];
        return g;
    }

    // Deterministic random genome without homopolymers (mergeTipsToEdge requires a properly collapsed junction).
    std::string CollapsedGenome(size_t len, uint32_t seed);

    // Owns a hasher and the SparseDBG built from a genome (every k-mer a vertex). SparseDBG is neither
    // copyable nor movable, so this holder is built in place.
    struct RandomDBG {
        hashing::RollingHash hasher;
        dbg::SparseDBG dbg;
        std::string genome;

        RandomDBG(std::string genome_, size_t k)
                : hasher(k), dbg(hasher), genome(std::move(genome_)) {
            Sequence gseq(genome);
            std::vector<hashing::htype> hashes;
            for (const hashing::KWH &kwh : hasher.kmers(gseq))
                hashes.push_back(kwh.hash());
            std::sort(hashes.begin(), hashes.end());
            hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
            for (hashing::htype h : hashes)
                dbg.addKmerVertex(h);
            std::vector<StringContig> contigs;
            contigs.emplace_back(genome, "genome");
            dbg::FillSparseDBGEdges(dbg, contigs.begin(), contigs.end(), QuietLogger(), 1, k);
        }
    };

    // Builds a SparseDBG from a handful of DESIGNED sequences (every k-mer becomes a vertex). Use this
    // for small, exact examples that must still be *proper* de Bruijn graphs: because every edge is a
    // genuine (k+1)-mer step, the splitting property holds automatically (a vertex's outgoing edges have
    // distinct first code nucleotides, in both the forward and the rc direction) - unlike GraphBuilder,
    // which joins arbitrary vertices and can produce code collisions at junctions.
    struct SeqDBG {
        hashing::RollingHash hasher;
        dbg::SparseDBG dbg;

        SeqDBG(const std::vector<std::string> &seqs, size_t k) : hasher(k), dbg(hasher) {
            std::vector<hashing::htype> hashes;
            for (const std::string &s : seqs) {
                Sequence gs(s);
                for (const hashing::KWH &kwh : hasher.kmers(gs))
                    hashes.push_back(kwh.hash());
            }
            std::sort(hashes.begin(), hashes.end());
            hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
            for (hashing::htype h : hashes)
                dbg.addKmerVertex(h);
            std::vector<StringContig> contigs;
            for (size_t i = 0; i < seqs.size(); i++)
                contigs.emplace_back(seqs[i], "s" + std::to_string(i));
            dbg::FillSparseDBGEdges(dbg, contigs.begin(), contigs.end(), QuietLogger(), 1, k);
        }

        ag::GraphPath align(const std::string &s, size_t w = 1) {
            dbg::KmerIndex index(dbg);
            index.fillAnchors(QuietLogger(), 1, dbg, w);
            return index.align(Sequence(s), "q");
        }
    };

    // Converts a (condensed) SparseDBG into a supregraph exactly the way Supregraph::FromDBG does: turn every
    // corporeal (non prefix/suffix) edge into a supre-vertex via edgeToSupreVertex. The graph's own
    // SampleSegmentTracker gives every new supre-vertex between two k-mers its kpomer chunk, as in the real
    // pipeline. Hashing is disabled first, as in the pipeline: supre-vertices are longer than k.
    inline void ConvertToSupregraph(dbg::SparseDBG &graph) {
        graph.disableHashing();
        std::vector<ag::EdgeId> eids;
        for (ag::Edge &e : graph.edgesUnique())
            eids.push_back(e.getId());
        for (ag::EdgeId eid : eids)
            if (!eid->isPrefix() && !eid->isSuffix())
                graph.edgeToSupreVertex(*eid);
    }

    // Resolves one 2-in/2-out k-mer junction of a supregraph by an arbitrary one-to-one matching (so a path may use
    // a pair the plan does not have). Returns false if there is no such junction.
    bool ResolveMatching(ag::AssemblyGraph &graph, size_t k);

    // Resolves a junction k-mer (at least 2 in, 2 out) with exactly the (in, out) pairs the stored reads use (the
    // multiplexing contract: no stored read is broken). Each read path through the pair, in either orientation, adds 1
    // to its support. Returns false if there is no such junction.
    bool ResolveByReads(ag::AssemblyGraph &graph, ag::AlignedReadStorage &storage, size_t k);

    // Flips every vertex with a single incoming suffix edge and out-degree other than 1; returns how many.
    size_t FlipAll(ag::AssemblyGraph &graph);

    inline size_t EdgeCount(const ag::GraphPath &path) {
        size_t n = 0;
        for (ag::Edge &e : path.edges()) { (void) e; n++; }
        return n;
    }

    // Extracts `count` reads of length `read_len` as substrings of the genome at random positions and
    // aligns each to the (current) graph, returning them as ag::AlignedReads. Must be called before the
    // graph is modified, since alignments reference the graph as it is now.
    inline std::vector<ag::AlignedRead> ExtractAndAlignReads(const RandomDBG &rdbg, size_t count,
                                                             size_t read_len, size_t w, uint32_t seed) {
        dbg::KmerIndex index(const_cast<dbg::SparseDBG &>(rdbg.dbg));
        index.fillAnchors(QuietLogger(), 1, const_cast<dbg::SparseDBG &>(rdbg.dbg), w);
        std::mt19937 rng(seed);
        VERIFY(rdbg.genome.size() > read_len);
        std::uniform_int_distribution<size_t> dist(0, rdbg.genome.size() - read_len);
        std::vector<ag::AlignedRead> reads;
        for (size_t i = 0; i < count; i++) {
            size_t pos = dist(rng);
            Sequence rseq(rdbg.genome.substr(pos, read_len));
            std::string name = "read" + std::to_string(i);
            reads.emplace_back(name, index.align(rseq, name));
        }
        return reads;
    }

    // ---- Invariant helpers --------------------------------------------------------------------

    // read id -> nucleotide sequence spelled by the read's path (respecting end cuts). The central
    // invariant under test: no graph operation may ever change this mapping for any read.
    const std::string LEGACY_READ = "<legacy>";

    inline std::map<std::string, std::string> ReadSeqSnapshot(ag::AlignedReadStorage &storage) {
        std::map<std::string, std::string> res;
        for (ag::AlignedRead &r : storage) {
            // Include subreads (empty graph path, i.e. a substring of a single vertex) as long as they
            // still spell a sequence -- a merge can legitimately turn a read into a subread of res.
            // A subread of a deleted vertex becomes legacy: its vertex is gone, so it spells nothing any more.
            if (r.valid() && r.getPath().isLegacy())
                res[r.getId()] = LEGACY_READ;
            else if (r.valid() && !r.getPath().Seq().empty())
                res[r.getId()] = r.getPath().Seq().str();
        }
        return res;
    }

    // Every read present before must still be present with an identical spelled sequence, or have become legacy.
    inline bool ReadSeqsUnchanged(const std::map<std::string, std::string> &before,
                                  ag::AlignedReadStorage &storage, std::string &err) {
        std::map<std::string, std::string> after = ReadSeqSnapshot(storage);
        for (const auto &kv : before) {
            auto it = after.find(kv.first);
            if (it == after.end()) {
                err = "read " + kv.first + " disappeared / became invalid";
                return false;
            }
            if (it->second != kv.second && it->second != LEGACY_READ) {
                err = "read " + kv.first + " sequence changed:\n  before=" + kv.second + "\n  after =" + it->second;
                return false;
            }
        }
        return true;
    }

    // Normalized content of one SuffixRecord: suffix sequence -> total positive count.
//    A stored suffix may end anywhere within its last edge's code (it was stored when that edge was shorter), so each is
//    brought to its canonical form first: the codes up to the first character of the last edge it reaches.
    inline std::map<std::string, int> NormalizeSuffixRecord(const ag::SuffixRecord &rec) {
        std::map<std::string, int> m;
        for (const std::pair<Sequence, int> &p : rec.getSuffixes()) {
            if (p.second <= 0)
                continue;
            ag::GraphPath path(rec.getEdge().getFinish(), p.first);
            Sequence codes(path.getFSplits());
            m[codes.Subseq(0, codes.size() - path.backEdge().getCode().size() + 1).str()] += p.second;
        }
        return m;
    }

    // "Refilled from the read collection remains the same": the live SuffixTracker (carried through
    // the graph operations) must store exactly what a tracker built fresh from the same storage/graph
    // stores, edge by edge.
    inline bool SuffixContentEqual(ag::SuffixTracker &live, ag::SuffixTracker &fresh,
                                   ag::AssemblyGraph &graph, std::string &err) {
        for (ag::Edge &e : graph.edges()) {
//            The live tracker may keep records the fresh one would not create; every fresh record must be there too.
            if (!fresh.hasSuffixRecord(e))
                continue;
            if (!live.hasSuffixRecord(e)) {
                std::stringstream ss;
                ss << "edge " << e.getId() << " has a fresh record but no live one";
                err = ss.str();
                return false;
            }
            std::map<std::string, int> a = NormalizeSuffixRecord(live.getSuffixRecord(e));
            std::map<std::string, int> b = NormalizeSuffixRecord(fresh.getSuffixRecord(e));
            if (a != b) {
                std::stringstream ss;
                ss << "SuffixRecord mismatch on edge " << e.getId() << "\n live:\n" << live.getSuffixRecord(e).str()
                   << " fresh:\n" << fresh.getSuffixRecord(e).str();
                err = ss.str();
                return false;
            }
        }
        return true;
    }

    // EdgeCodeListener invariant: getOutgoing(c) matches an outgoing edge by getCode()[0], so within a
    // vertex all (corporeal, non-suffix) outgoing edges must carry a non-empty code with a distinct
    // first nucleotide.
    inline bool EdgeCodesValid(ag::AssemblyGraph &graph, std::string &err) {
        for (ag::Vertex &v : graph.vertices()) {
            std::set<unsigned char> seen;
            for (ag::Edge &e : v) {
                if (e.isSuffix())
                    continue;
                if (e.getCode().empty()) {
                    std::stringstream ss;
                    ss << "empty code on non-suffix edge " << e.getId();
                    err = ss.str();
                    return false;
                }
                unsigned char c = e.getCode()[0];
                if (seen.count(c)) {
                    err = "two outgoing edges share first code nucleotide at a vertex";
                    return false;
                }
                seen.insert(c);
            }
        }
        return true;
    }
}
