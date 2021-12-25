//
// Created by Andrey Bzikadze on 12/10/21.
//

#include "repeat_resolution/mdbg_seq.hpp"
#include "repeat_resolution/mdbg_topology.hpp"
#include "gtest/gtest.h"

using namespace repeat_resolution;
using namespace dbg;

/*
TEST(RC, Basic) {
  Edge
  ASSERT_EQ(MDBGSeq(Sequence("AATTCCGG")).RC(), MDBGSeq(Sequence("CCGGAATT")));
  ASSERT_EQ(MDBGSeq().RC(), MDBGSeq());

  // ASSERT_TRUE(MDBGSeq("AATTCCGG").IsCanonical());
  // ASSERT_FALSE(MDBGSeq("CCGGAATT").IsCanonical());
  // ASSERT_TRUE(MDBGSeq("ACGT").IsCanonical());
}
 */

TEST(EdgeSegment, Basic) {
    Sequence seq("AATTCCGG");
    uint64_t k = 2;
    Vertex st(0);
    st.seq = seq.Prefix(k);
    dbg::Edge edge(&st, nullptr, seq.Suffix(seq.size() - k));
    EdgeSegment segm(&edge, 0, seq.size());
    ASSERT_TRUE(segm.LeftFull());
    ASSERT_TRUE(segm.RightFull());
    ASSERT_EQ(segm.ToSequence(), seq);

    int trim = 1;
    segm.TrimLeft(trim);
    ASSERT_FALSE(segm.LeftFull());
    ASSERT_TRUE(segm.start==trim);
    ASSERT_TRUE(segm.ToSequence()==seq.Subseq(trim));

    EdgeSegment segm2(&edge, 0, trim);
    ASSERT_TRUE(segm2.LeftFull());
    ASSERT_TRUE(segm2.ToSequence()==seq.Prefix(trim));

    segm2.ExtendRight(segm);
    ASSERT_TRUE(segm2.LeftFull());
    ASSERT_TRUE(segm2.RightFull());

    segm2.TrimRight(trim);
    ASSERT_FALSE(segm2.RightFull());
    ASSERT_TRUE(segm2.end==seq.size() - trim);
    ASSERT_TRUE(segm2.ToSequence()==seq.Prefix(seq.size() - trim));

    EdgeSegment segm3(&edge, seq.size() - trim, seq.size());
    ASSERT_TRUE(segm3.ToSequence()==seq.Suffix(trim));
    segm2.ExtendRight(segm3);
    ASSERT_TRUE(segm2.RightFull());
    ASSERT_TRUE(segm2.ToSequence()==seq);
}

TEST(MDBGSeq, SingleSegment) {
    Sequence s1("ACGT");
    uint64_t k = 2;
    Vertex st(0);
    st.seq = s1.Prefix(k);
    dbg::Edge edge1(&st, nullptr, s1.Subseq(k));
    MDBGSeq seq1(&edge1, 0, s1.size());
    ASSERT_EQ(seq1.ToSequence(), s1);
    ASSERT_EQ(seq1.Size(), s1.size());
    ASSERT_FALSE(seq1.Empty());
    ASSERT_EQ(seq1.ContainerSize(), 1);

    int trim{1};
    seq1.TrimLeft(trim);
    ASSERT_EQ(seq1.ToSequence(), s1.Subseq(trim));
    seq1.TrimRight(trim);
    ASSERT_EQ(seq1.ToSequence(), s1.Subseq(trim, s1.size() - 1));

    MDBGSeq seq1L(&edge1, 0, trim);
    seq1.Prepend(std::move(seq1L));
    MDBGSeq seq1R(&edge1, s1.size() - 1, s1.size());
    seq1.Append(std::move(seq1R));
    ASSERT_EQ(seq1.ToSequence(), s1);
    ASSERT_EQ(seq1.Size(), s1.size());
    ASSERT_FALSE(seq1.Empty());
    ASSERT_EQ(seq1.ContainerSize(), 1);
}

TEST(MDBGSeq, SeveralSegments) {
    uint64_t k = 2;

    Sequence s1("ACGT");
    Vertex st1(0);
    st1.seq = s1.Prefix(k);
    dbg::Edge edge1(&st1, nullptr, s1.Subseq(k));
    MDBGSeq seq1(&edge1, 0, s1.size());

    Sequence s2("ACC");
    Vertex st2(0);
    st2.seq = s2.Prefix(k);
    dbg::Edge edge2(&st2, nullptr, s2.Subseq(k));
    MDBGSeq seq2(&edge2, 0, s2.size());
    seq1.Append(std::move(seq2));

    Sequence s3("TGA");
    Vertex st3(0);
    st3.seq = s3.Prefix(k);
    dbg::Edge edge3(&st3, nullptr, s3.Subseq(k));
    MDBGSeq seq3(&edge3, 0, s3.size());
    seq1.Prepend(std::move(seq3));

    Sequence concat = Sequence::Concat({s3, s1, s2});
    ASSERT_EQ(seq1.ToSequence(), concat);
    ASSERT_EQ(seq1.Size(), concat.size());
    ASSERT_EQ(seq1.ContainerSize(), 3);

    for (size_t i = 0; i < concat.size(); ++i) {
        for (size_t j = i + 1; j < concat.size(); ++j) {
            // std::cout << i << " " << j << "\n";
            ASSERT_EQ(seq1.Substr(i, j - i).ToSequence(), concat.Subseq(i, j));
        }
    }

    int trim1 = 1;
    seq1.TrimRight(trim1);
    ASSERT_EQ(seq1.ToSequence(), concat.Subseq(0, concat.size() - trim1));
    ASSERT_EQ(seq1.Size(), concat.size() - trim1);
    ASSERT_EQ(seq1.ContainerSize(), 3);

    int trim2 = 2;
    seq1.TrimRight(trim2);
    ASSERT_EQ(seq1.ToSequence(),
              concat.Subseq(0, concat.size() - trim1 - trim2));
    ASSERT_EQ(seq1.Size(), concat.size() - trim1 - trim2);
    ASSERT_EQ(seq1.ContainerSize(), 2);

    seq1.TrimLeft(trim1);
    ASSERT_EQ(seq1.ToSequence(),
              concat.Subseq(trim1, concat.size() - trim1 - trim2));
    ASSERT_EQ(seq1.Size(), concat.size() - trim1 - trim2 - trim1);
    ASSERT_EQ(seq1.ContainerSize(), 2);

    seq1.TrimLeft(trim2);
    ASSERT_EQ(seq1.ToSequence(),
              concat.Subseq(trim1 + trim2, concat.size() - trim1 - trim2));
    ASSERT_EQ(seq1.Size(), concat.size() - trim1 - trim2 - trim1 - trim2);
    ASSERT_EQ(seq1.ContainerSize(), 1);
}