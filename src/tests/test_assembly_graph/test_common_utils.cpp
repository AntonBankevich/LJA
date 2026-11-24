// Small low-level utilities: number printing and parsing, sequence reading.

#include "gtest/gtest.h"
#include "common/hash_utils.hpp"
#include "sequences/read_loader.hpp"
#include "sequences/sequence.hpp"
#include "alignment/ksw_aligner.hpp"
#include <fstream>
#include <sstream>

// A hash of 0 (e.g. the poly-A k-mer) used to be printed as nothing, which shifted every later entry of a saved hash
// list when it was read back.
TEST(HashUtils, HashesRoundTripThroughStreams) {
    hashing::htype big = (hashing::htype(0x0123456789abcdefULL) << 64) | hashing::htype(0xfedcba9876543210ULL);
    std::vector<hashing::htype> values = {0, 1, 9, 10, 1234567890123456789ULL, big, 0};
    std::stringstream ss;
    for (hashing::htype v : values)
        ss << v << " " << true << "\n";
    for (hashing::htype expected : values) {
        hashing::htype v;
        bool flag;
        ss >> v >> flag;
        EXPECT_TRUE(v == expected);
        EXPECT_TRUE(flag);
    }
    std::stringstream zero;
    zero << hashing::htype(0);
    EXPECT_EQ(zero.str(), "0");
}

namespace {
    std::vector<Contig> ReadAll(const std::string &content, const std::string &extension) {
        std::experimental::filesystem::path file = std::experimental::filesystem::temp_directory_path() /
                ("lja_seqio_test" + extension);
        {
            std::ofstream os(file);
            os << content;
        }
        logging::Logger logger(false);
        std::vector<Contig> res = io::ReadLoader(logger, 1).readAll({file});
        std::experimental::filesystem::remove(file);
        return res;
    }
}

// Windows line endings and blank lines at the end of the file (an empty line used to be read with back()).
TEST(SeqIO, ReadsFilesWithCarriageReturnsAndTrailingBlankLines) {
    for (const std::string &tail : {std::string(""), std::string("\n"), std::string("\n\n"), std::string("\r\n")}) {
        std::vector<Contig> fasta = ReadAll(">a\r\nACGT\r\nTT\r\n>b\r\nGGG\r\n" + tail, ".fasta");
        ASSERT_EQ(fasta.size(), 2u) << "tail=" << tail.size();
        EXPECT_EQ(fasta[0].getInnerId(), "a");
        EXPECT_EQ(fasta[0].getSeq().str(), "ACGTTT");
        EXPECT_EQ(fasta[1].getSeq().str(), "GGG");
        std::vector<Contig> fastq = ReadAll("@a\r\nACGT\r\n+\r\nIIII\r\n@b\r\nGG\r\n+\r\nII\r\n" + tail, ".fastq");
        ASSERT_EQ(fastq.size(), 2u) << "tail=" << tail.size();
        EXPECT_EQ(fastq[0].getSeq().str(), "ACGT");
        EXPECT_EQ(fastq[1].getSeq().str(), "GG");
    }
}

// Sequence views in reverse-complement orientation (rtl): copy(), asNumber() and endsWith() used to be wrong for them.
TEST(SequenceTest, RcViewsAndBounds) {
    Sequence s("ACGTTGCAAC");
    Sequence rc = !s;
    EXPECT_EQ(rc.copy().str(), rc.str()) << "copy of an rc view keeps its orientation";
    EXPECT_EQ(rc.asNumber(), Sequence(rc.str()).asNumber());
    EXPECT_EQ(s.Subseq(2, 7).asNumber(), Sequence(s.str().substr(2, 5)).asNumber());
    EXPECT_FALSE(Sequence("ACG").endsWith(Sequence("TACG"))) << "a longer sequence is not a suffix";
    EXPECT_TRUE(Sequence("TACG").endsWith(Sequence("ACG")));
    EXPECT_TRUE(rc.endsWith(Sequence(rc.str().substr(4))));
}

// KSWAligner::score with the default parameters (match +1, mismatch -5, gap open 10, extension 2): match runs used to
// be charged two gap penalties each.
TEST(KSWAligner, ScoreChargesGapsOnlyForGaps) {
    KSWAligner aligner;
    EXPECT_EQ(aligner.score("ACGTACGT", "ACGTACGT", AlignmentForm("8M")), 8);
    EXPECT_EQ(aligner.score("ACGTACGT", "ACGTTCGT", AlignmentForm("8M")), 7 - 5);
//    The query has one extra nucleotide (I), the target two (D).
    EXPECT_EQ(aligner.score("ACGTACGT", "ACGTTACGT", AlignmentForm("4M1I4M")), 8 - 10);
    EXPECT_EQ(aligner.score("ACGTGGACGT", "ACGTACGT", AlignmentForm("4M2D4M")), 8 - (10 + 2));
}
