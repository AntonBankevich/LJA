#include "sequences/read_loader.hpp"
#include "dbg/aln_reads_reader.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
#include "assembly_graph/visualization.hpp"
#include "random_graph_utils.hpp"
#include <gtest/gtest.h>
#include <zlib.h>
#include <mutex>
#include <random>

using namespace ag_test;

namespace {
    namespace fs = std::experimental::filesystem;
    typedef std::vector<std::pair<std::string, std::string>> Reads; // (id, sequence)

    const CompressionParameters compressing{true, 6, 10, 1};

    fs::path WriteFile(const std::string &name, const std::string &text) {
        fs::path path = fs::temp_directory_path() / ("lja_read_loader_test_" + name);
        if (endsWith(name, ".gz")) {
            gzFile file = gzopen(path.c_str(), "wb");
            gzwrite(file, text.data(), unsigned(text.size()));
            gzclose(file);
        } else {
            std::ofstream os(path);
            os << text;
        }
        return path;
    }

//    The reads a loader must give for the written ones, computed independently of the loader: split as described
//    for io::ReadLoader, then compressed.
    Reads Expected(const Reads &written, size_t min_read_size, size_t overlap,
                   const CompressionParameters &compression = CompressionParameters::None()) {
        Reads res;
        for (const auto &read : written) {
            const std::string &seq = read.second;
            for (size_t start = 0;; start += min_read_size - overlap) {
                bool last = seq.size() <= 2 * min_read_size - overlap + start;
                std::string piece = seq.substr(start, last ? std::string::npos : min_read_size);
                bool whole = start == 0 && last;
                res.emplace_back(whole ? read.first : read.first + "_" + std::to_string(start),
                                 Sequence(piece).compressed(compression).str());
                if (last)
                    break;
            }
        }
        return res;
    }

//    The reads by index, checking that every index comes once.
    Reads Load(io::ReadLoader &loader, const io::Library &lib) {
        std::mutex mutex;
        std::map<size_t, std::pair<std::string, std::string>> loaded;
        size_t cnt = loader.process(lib, [&](size_t pos, Contig &contig) {
            std::lock_guard<std::mutex> guard(mutex);
            bool inserted = loaded.emplace(pos, std::make_pair(contig.getInnerId(), contig.getSeq().str())).second;
            EXPECT_TRUE(inserted) << "index " << pos << " given twice";
        });
        Reads res;
        for (auto &p : loaded) {
            EXPECT_EQ(p.first, res.size()) << "indices are not consecutive";
            res.push_back(std::move(p.second));
        }
        EXPECT_EQ(cnt, res.size());
        return res;
    }

    Reads Load(const io::Library &lib, size_t threads, size_t min_read_size, size_t overlap, size_t buffer_size,
               const CompressionParameters &compression = CompressionParameters::None()) {
        io::ReadLoader loader(QuietLogger(), threads, compression, io::Splitter(min_read_size, overlap), buffer_size);
        return Load(loader, lib);
    }

    std::string RandomRead(std::mt19937 &rng, size_t max_len) {
        static const char bases[] = {'A', 'C', 'G', 'T'};
        std::uniform_int_distribution<size_t> len_dist(1, max_len);
        std::uniform_int_distribution<int> base(0, 3);
        std::uniform_int_distribution<int> run(1, 4);
        std::uniform_int_distribution<int> kind(0, 9);
        size_t len = len_dist(rng);
        std::string res;
//        Homopolymer runs and some dimer repeats, so that compression has something to do.
        while (res.size() < len) {
            char a = bases[base(rng)];
            if (kind(rng) == 0) {
                char b = bases[base(rng)];
                for (int i = 0, n = 3 + run(rng) * 2; i < n; i++)
                    res += i % 2 == 0 ? a : b;
            } else {
                res += std::string(run(rng), a);
            }
        }
        res.resize(len);
        return res;
    }

    std::string LowerSome(std::mt19937 &rng, std::string s) {
        for (char &c : s)
            if (rng() % 10 == 0)
                c = char(std::tolower(c));
        return s;
    }

    std::string Lines(std::mt19937 &rng, const std::string &s, const std::string &eol) {
        std::uniform_int_distribution<size_t> width(1, 80);
        size_t w = width(rng) < 20 ? width(rng) : s.size();
        std::string res;
        for (size_t i = 0; i < s.size(); i += w)
            res += s.substr(i, w) + eol;
        return res;
    }

//    Read i with a comment or not; its id is read<i>.
    std::string Header(std::mt19937 &rng, size_t i) {
        std::string res = "read" + std::to_string(i);
        if (rng() % 2 == 0)
            res += " some comment " + std::to_string(rng() % 100);
        return res;
    }

    std::string RandomFasta(uint32_t seed, size_t n, size_t max_len, bool crlf, Reads &written) {
        std::mt19937 rng(seed);
        const std::string eol = crlf ? "\r\n" : "\n";
        std::string res;
        for (size_t i = 0; i < n; i++) {
            std::string seq = RandomRead(rng, max_len);
            written.emplace_back("read" + std::to_string(i), seq);
            res += ">" + Header(rng, i) + eol + Lines(rng, LowerSome(rng, seq), eol);
            if (rng() % 5 == 0)
                res += eol;
        }
        if (rng() % 2 == 0)
            res.pop_back();
        return res;
    }

    std::string RandomFastq(uint32_t seed, size_t n, size_t max_len, bool crlf, Reads &written) {
        std::mt19937 rng(seed);
        const std::string eol = crlf ? "\r\n" : "\n";
        std::uniform_int_distribution<int> qual('!', '~');
        std::string res;
        for (size_t i = 0; i < n; i++) {
            std::string seq = RandomRead(rng, max_len);
            written.emplace_back("read" + std::to_string(i), seq);
            std::string q;
            for (size_t j = 0; j < seq.size(); j++)
                q += char(qual(rng));
//            Quality lines starting with '@' or '+' must not be taken for record starts.
            if (rng() % 3 == 0)
                q[0] = rng() % 2 == 0 ? '@' : '+';
            res += "@" + Header(rng, i) + eol + Lines(rng, LowerSome(rng, seq), eol) + "+" + eol + Lines(rng, q, eol);
        }
        if (rng() % 2 == 0)
            res.pop_back();
        return res;
    }

    std::string RandomGfa(uint32_t seed, size_t n, size_t max_len, Reads &written) {
        std::mt19937 rng(seed);
        std::string res = "H\tVN:Z:1.0\n";
        for (size_t i = 0; i < n; i++) {
            std::string seq = RandomRead(rng, max_len);
            written.emplace_back("seg" + std::to_string(i), seq);
            res += "S\tseg" + std::to_string(i) + "\t" + LowerSome(rng, seq);
            res += rng() % 2 == 0 ? "\tLN:i:" + std::to_string(seq.size()) + "\n" : "\n";
            if (rng() % 3 == 0)
                res += "L\tseg" + std::to_string(i) + "\t+\tseg0\t-\t0M\n";
        }
        return res;
    }

    io::Library RandomLibrary(uint32_t seed, Reads &written) {
        return {WriteFile("a.fasta", RandomFasta(seed, 200, 300, false, written)),
                WriteFile("b.fa.gz", RandomFasta(seed + 1, 200, 300, true, written)),
                WriteFile("c.fastq", RandomFastq(seed + 2, 200, 300, false, written)),
                WriteFile("d.fq.gz", RandomFastq(seed + 3, 200, 300, true, written)),
                WriteFile("e.gfa", RandomGfa(seed + 4, 200, 300, written))};
    }

//    The loader gives the expected reads of the written ones for any threads and buffer sizes.
    void CheckReads(const io::Library &lib, const Reads &written, size_t min_read_size, size_t overlap,
                    const CompressionParameters &compression = CompressionParameters::None()) {
        const Reads expected = Expected(written, min_read_size, overlap, compression);
        for (size_t threads : {1, 4}) {
            for (size_t buffer_size : {size_t(7), size_t(64), size_t(1000), size_t(1) << 20}) {
                ASSERT_EQ(Load(lib, threads, min_read_size, overlap, buffer_size, compression), expected)
                                            << threads << " threads, buffer " << buffer_size;
            }
        }
    }
}

TEST(ReadLoader, GivesWrittenReads) {
    Reads written;
    io::Library lib = RandomLibrary(1, written);
    CheckReads(lib, written, size_t(-1) / 2, size_t(-1) / 8);
}

TEST(ReadLoader, SplitsLongReads) {
    Reads written;
    io::Library lib = RandomLibrary(11, written);
    CheckReads(lib, written, 50, 10);
}

TEST(ReadLoader, CompressesReads) {
    Reads written;
    io::Library lib = RandomLibrary(21, written);
    CheckReads(lib, written, size_t(-1) / 2, size_t(-1) / 8, compressing);
    CheckReads(lib, written, 50, 10, compressing);
}

//    The compression parameters are the global ones unless set.
TEST(ReadLoader, TakesGlobalCompression) {
    io::Library lib = {WriteFile("global.fasta", ">a\nAAACCCGT\n")};
    CompressionParameters saved = CompressionParameters::global;
    CompressionParameters::global = compressing;
    io::ReadLoader global_loader(QuietLogger(), 1);
    io::ReadLoader raw_loader(QuietLogger(), 1, CompressionParameters::None());
    Reads global = Load(global_loader, lib);
    Reads raw = Load(raw_loader, lib);
    CompressionParameters::global = saved;
    ASSERT_EQ(global, Reads({{"a", "ACGT"}}));
    ASSERT_EQ(raw, Reads({{"a", "AAACCCGT"}}));
}

//    One long read cut into many pieces that threads pack side by side in the same text, on one line and wrapped.
TEST(ReadLoader, SplitsVeryLongRead) {
    std::mt19937 rng(41);
    std::string seq;
    while (seq.size() < 200000)
        seq += RandomRead(rng, 300);
    std::string wrapped;
    for (size_t i = 0; i < seq.size(); i += 61)
        wrapped += seq.substr(i, 61) + "\n";
    io::Library lib = {WriteFile("long1.fasta", ">long one line\n" + seq + "\n"),
                       WriteFile("long2.fa.gz", ">long wrapped\n" + wrapped)};
    Reads written = {{"long", seq}, {"long", seq}};
    for (const CompressionParameters &compression : {CompressionParameters::None(), compressing}) {
        const Reads expected = Expected(written, 1000, 360, compression);
        ASSERT_GT(expected.size(), 500);
        for (size_t threads : {1, 8}) {
            for (size_t buffer_size : {size_t(1000), size_t(1) << 20})
                ASSERT_EQ(Load(lib, threads, 1000, 360, buffer_size, compression), expected);
        }
    }
}

//    With one thread, tasks run in the order of the reads.
TEST(ReadLoader, OneThreadKeepsOrder) {
    Reads written;
    io::Library lib = RandomLibrary(31, written);
    std::vector<size_t> order;
    io::ReadLoader(QuietLogger(), 1, CompressionParameters::None(), io::Splitter(50, 10), 1000).process(lib,
            [&order](size_t pos, Contig &) {
        order.push_back(pos);
    });
    ASSERT_EQ(order.size(), Expected(written, 50, 10).size());
    for (size_t i = 0; i < order.size(); i++)
        ASSERT_EQ(order[i], i);
}

TEST(ReadLoader, ReadsAll) {
    Reads written;
    io::Library lib = RandomLibrary(51, written);
    std::vector<Contig> all = io::ReadLoader(QuietLogger(), 4).readAll(lib);
    Reads expected = Expected(written, size_t(-1) / 2, size_t(-1) / 8);
    ASSERT_EQ(all.size(), expected.size());
    for (size_t i = 0; i < all.size(); i++) {
        ASSERT_EQ(all[i].getInnerId(), expected[i].first);
        ASSERT_EQ(all[i].getSeq().str(), expected[i].second);
    }
}

//    Reads aligned to a graph, given by a .paths file, are the sequences of their paths, not compressed.
TEST(ReadLoader, ReadsPaths) {
    const size_t k = 15;
    RandomDBG rdbg(RandomGenome(5000, 61), k);
    ag::AlignedReadStorage storage(rdbg.dbg, ExtractAndAlignReads(rdbg, 100, 200, 5, 62));
    Reads written;
    for (const ag::AlignedRead &read : storage)
        written.emplace_back(read.getId(), read.getPath().Seq().str());
    fs::path gfa = fs::temp_directory_path() / "lja_read_loader_test_graph.gfa";
    ag::Printer printer;
    printer.setEdgeInfo(ag::EdgeInfo({&ag::GetEdgeNameForSaving}, {}, {}));
    printer.printGFA(gfa, ag::Component(rdbg.dbg), true);
    fs::path paths = fs::temp_directory_path() / "lja_read_loader_test_reads.paths";
    storage.printReadPaths(QuietLogger(), fs::temp_directory_path() / "lja_read_loader_test_reads.aln", gfa, paths, k);
    io::Library lib = {paths, WriteFile("paths_mixed.fasta", ">x\nACGTT\n")};
    written.emplace_back("x", "ACGTT");
    for (size_t threads : {1, 4}) {
        dbg::ReadLoader loader(QuietLogger(), threads, compressing, io::Splitter(50, 10), 1000);
        Reads expected = Expected(written, 50, 10);
        expected.back().second = "ACGT";
        ASSERT_EQ(Load(loader, lib), expected);
    }
}

//    The text only a piece reads must hold the piece packed.
TEST(ReadLoader, RejectsLargeOverlap) {
    io::Library lib = {WriteFile("overlap.fasta", ">a\nACGT\n")};
    EXPECT_DEATH(Load(lib, 1, 1000, 370, 1000), "");
}

//    Sequences handed to tasks lie in the read buffer, so a task keeping one is an error.
TEST(ReadLoader, RejectsKeptSequence) {
    io::Library lib = {WriteFile("kept.fasta", ">a\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>b\nACGT\n")};
    std::vector<Sequence> kept;
    io::ReadLoader loader(QuietLogger(), 1);
    EXPECT_DEATH(loader.process(lib, [&kept](size_t, Contig &contig) {kept.push_back(contig.getSeq());}), "");
    size_t cnt = loader.process(lib, [&kept](size_t, Contig &contig) {kept.push_back(contig.getSeq().copy());});
    ASSERT_EQ(cnt, 2);
    ASSERT_EQ(kept.back().str(), "ACGT");
}

TEST(ReadLoader, PreallocatedNuclBuffer) {
    std::vector<NuclBuffer::ST> words(2, 0);
    llvm::IntrusiveRefCntPtr<ManagedNuclBuffer> buffer(new ManagedNuclBuffer(40));
    ASSERT_TRUE(buffer->ownsWords());
    buffer->usePreallocated(words.data(), 40);
    ASSERT_FALSE(buffer->ownsWords());
    buffer->setNucl(33, 2);
    ASSERT_EQ(words[1], NuclBuffer::ST(2) << 2u);
    {
        Sequence seq(buffer);
        ASSERT_EQ(seq.size(), 40);
        ASSERT_EQ(buffer->UseCount(), 2);
        ASSERT_EQ(seq[33], 2);
    }
    ASSERT_EQ(buffer->UseCount(), 1);
    buffer = nullptr;
    ASSERT_EQ(words[1], NuclBuffer::ST(2) << 2u);
}

//    Compressing a sequence directly, including an rc view, as characters would be compressed.
TEST(ReadLoader, CompressesSequence) {
    std::mt19937 rng(71);
    for (size_t i = 0; i < 100; i++) {
        std::string s = RandomRead(rng, 300);
        Sequence seq(s);
        std::string rc = (!seq).str();
        struct Chars {
            std::string &s;
            char get(size_t i) const {return s[i];}
            void set(size_t i, char c) {s[i] = c;}
        } chars{s}, rc_chars{rc};
        s.resize(CompressInPlace(chars, s.size(), compressing));
        rc.resize(CompressInPlace(rc_chars, rc.size(), compressing));
        ASSERT_EQ(seq.compressed(compressing).str(), s);
        ASSERT_EQ((!seq).compressed(compressing).str(), rc);
    }
}

//    AlignReads gives each read, in order, its alignment by the k-mer index.
TEST(ReadLoader, AlignReads) {
    const size_t k = 15;
    const size_t w = 5;
    RandomDBG rdbg(RandomGenome(20000, 51), k);
    std::mt19937 rng(52);
    std::uniform_int_distribution<size_t> len_dist(10, 400);
    std::string fastq;
    Reads written;
    for (size_t i = 0; i < 300; i++) {
        size_t len = len_dist(rng);
        std::string seq = rdbg.genome.substr(rng() % (rdbg.genome.size() - len), len);
        if (rng() % 2 == 0)
            seq = Sequence(seq).rc().str();
        written.emplace_back("read" + std::to_string(i), seq);
        fastq += "@read" + std::to_string(i) + "\n" + seq + "\n+\n" + std::string(len, 'I') + "\n";
    }
    std::vector<ag::AlignedRead> loaded = dbg::AlignReads(QuietLogger(), 4, {WriteFile("align.fastq", fastq)},
                                                          rdbg.dbg, w);
    dbg::KmerIndex index(rdbg.dbg);
    index.fillAnchors(QuietLogger(), 1, rdbg.dbg, w);
    ASSERT_EQ(loaded.size(), written.size());
    size_t aligned = 0;
    for (size_t i = 0; i < written.size(); i++) {
        ASSERT_EQ(loaded[i].getId(), written[i].first);
        if (written[i].second.size() < index.minReadLen()) {
            ASSERT_FALSE(loaded[i].valid());
            continue;
        }
        ASSERT_TRUE(loaded[i].getPath() == index.align(Sequence(written[i].second), written[i].first));
        aligned++;
    }
    ASSERT_GT(aligned, 250);
}

TEST(ReadLoader, SkipsEmptyReads) {
    io::Library lib = {WriteFile("empty.fasta", ">a\nACGT\n>b\n>c x\nggt\n"),
                       WriteFile("empty.fastq", "@d\nACGT\n+\nIIII\n@e\n\n+\n\n@f\nGG\n+\nII\n"),
                       WriteFile("empty.gfa", "S\tg\t\nS\th\tAC\n")};
    Reads expected = {{"a", "ACGT"}, {"c", "GGT"}, {"d", "ACGT"}, {"f", "GG"}, {"h", "AC"}};
    for (size_t buffer_size : {size_t(3), size_t(1000)})
        ASSERT_EQ(Load(lib, 2, size_t(-1) / 2, size_t(-1) / 8, buffer_size), expected);
}
