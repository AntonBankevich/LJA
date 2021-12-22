//
// Created by Andrey Bzikadze on 12/10/21.
//

#include "mdbg_seq.hpp"

using namespace repeat_resolution;

// ---------- EdgeSegment ----------

EdgeSegment::EdgeSegment(const dbg::Edge *edge, const uint64_t start,
                         const uint64_t end)
    : edge{edge}, start{start}, end{end} {
    VERIFY(start <= end);
    VERIFY(edge!=nullptr);
    VERIFY(end <= edge->size() + GetStK());
}

bool EdgeSegment::RightFull() const {
    VERIFY(edge!=nullptr);
    size_t k = edge->start()->seq.size();
    return end==GetStK() + edge->size();
}

void EdgeSegment::TrimLeft(const uint64_t size) {
    VERIFY(Size() >= size);
    start += size;
}

void EdgeSegment::TrimRight(const uint64_t size) {
    VERIFY(Size() >= size);
    end -= size;
}

uint64_t EdgeSegment::Size() const {
    VERIFY(start <= end);
    return end - start;
}

size_t EdgeSegment::EdgeSize() const {
    VERIFY(edge!=nullptr);
    return edge->size() + GetStK();
}

void EdgeSegment::ExtendRight(const EdgeSegment &segment) {
    VERIFY(edge!=nullptr and segment.edge!=nullptr)
    VERIFY(*edge==*segment.edge);
    VERIFY(end==segment.start);
    end = segment.end;
}

EdgeSegment EdgeSegment::RC() const {
    VERIFY(edge!=nullptr);
    return EdgeSegment(&edge->rc(), edge->size() + GetStK() - end,
                       edge->size() + GetStK() - start);
}

Sequence EdgeSegment::ToSequence() const {
    VERIFY(edge!=nullptr);
    return edge->suffix(0).Subseq(start, end);
}

double EdgeSegment::Cov() const {
    VERIFY(edge!=nullptr);
    return edge->getCoverage();
}

bool EdgeSegment::operator==(const EdgeSegment &rhs) const {
    return edge==rhs.edge and start==rhs.start and end==rhs.end;
}

// ---------- MDBGSeq ----------

double MDBGSeq::CovListEdgeSegm(const std::list<EdgeSegment> &segms) {
    if (segms.empty()) {
        return 0;
    }

    double cov{0};
    for (const EdgeSegment &segm : segms) {
        cov += segm.Cov()*segm.Size();
    }
    return cov;
}

void MDBGSeq::Swap(MDBGSeq &lhs, MDBGSeq &rhs) {
    std::swap(lhs.segms, rhs.segms);
    std::swap(lhs.cov, rhs.cov);
    std::swap(lhs.size, rhs.size);
}

MDBGSeq::MDBGSeq(const dbg::Edge *edge,
                 const uint64_t start,
                 const uint64_t end) :
    cov{edge->getCoverage()*(end - start)}, size{end - start} {
    segms.emplace_back(edge, start, end);
}

MDBGSeq::MDBGSeq(std::list<EdgeSegment> segms) : segms{std::move(segms)},
                                                 cov{CovListEdgeSegm(this->segms)} {
    size = 0;
    for (const EdgeSegment &segm : this->segms) {
        size += segm.Size();
    }
}

[[nodiscard]] Sequence MDBGSeq::ToSequence() const {
    std::vector<Sequence> sec_vec;
    for (const EdgeSegment segm : segms) {
        sec_vec.emplace_back(segm.ToSequence());
    }
    return Sequence::Concat(sec_vec);
}

[[nodiscard]] size_t MDBGSeq::Size() const {
    // {
    //     size_t size_{0};
    //     for (const EdgeSegment segm : segms) {
    //         size_ += segm.Size();
    //     }
    //     VERIFY(size==size_);
    // }
    return size;
}

[[nodiscard]] size_t MDBGSeq::ContainerSize() const { return segms.size(); }

[[nodiscard]] MDBGSeq MDBGSeq::RC() const {
    std::list<EdgeSegment> segms_rc;
    for (const EdgeSegment segm : segms) {
        segms_rc.emplace_front(segm.RC());
    }
    return MDBGSeq(std::move(segms_rc));
}

[[nodiscard]] bool MDBGSeq::IsCanonical() const {
    const Sequence seq = ToSequence();
    return seq <= !seq;
}

[[nodiscard]] bool MDBGSeq::Empty() const { return segms.empty(); }

void MDBGSeq::Append(MDBGSeq mdbg_seq) {
    if (mdbg_seq.Empty()) {
        return;
    }
    if (Empty()) {
        Swap(*this, mdbg_seq);
        return;
    }
    EdgeSegment &back = segms.back();
    EdgeSegment &front = mdbg_seq.segms.front();

    cov += mdbg_seq.cov;

    if (back.edge==front.edge and back.end==front.start) {
        back.ExtendRight(front);
        size += front.Size();

        mdbg_seq.TrimLeft(front.Size());
    }

    size += mdbg_seq.Size();
    segms.splice(segms.end(), std::move(mdbg_seq.segms));

    // std::cout << CovListEdgeSegm(segms) << " " << cov << "\n";
    // VERIFY(std::fabs(CovListEdgeSegm(segms) - cov) < 0.0001);
}

void MDBGSeq::Prepend(MDBGSeq mdbg_seq) {
    MDBGSeq temp;
    Swap(temp, *this);
    mdbg_seq.Append(std::move(temp));
    Swap(mdbg_seq, *this);
}

void MDBGSeq::TrimLeft(uint64_t size_) {
    VERIFY(size_ <= Size());
    size -= size_;
    while (size_ > 0) {
        EdgeSegment &front = segms.front();
        if (front.Size() <= size_) {
            size_ -= front.Size();
            cov -= front.Cov()*front.Size();
            segms.pop_front();
        } else {
            front.TrimLeft(size_);
            cov -= front.Cov()*size_;
            size_ = 0;
        }
    }

    // std::cout << cov << " " << CovListEdgeSegm(segms) << "\n";
    // VERIFY(std::fabs(CovListEdgeSegm(segms) - cov) < 0.0001);
}

void MDBGSeq::TrimRight(uint64_t size_) {
    VERIFY(size_ <= Size());
    size -= size_;
    while (size_ > 0) {
        EdgeSegment &back = segms.back();
        if (back.Size() <= size_) {
            size_ -= back.Size();
            cov -= back.Cov()*back.Size();
            segms.pop_back();
        } else {
            back.TrimRight(size_);
            cov -= back.Cov()*size_;
            size_ = 0;
        }
    }
    // std::cout << cov << " " << CovListEdgeSegm(segms) << "\n";
    // VERIFY(std::fabs(CovListEdgeSegm(segms) - cov) < 0.00001);
}

[[nodiscard]] MDBGSeq MDBGSeq::Substr(uint64_t pos, const uint64_t len) const {
    VERIFY(pos + len <= Size());
    if (len==0 or Empty()) {
        return MDBGSeq();
    }

    auto left = segms.begin();
    while (left->Size() <= pos) {
        pos -= left->Size();
        ++left;
        VERIFY(left!=segms.end());
    }
    auto right = left;
    uint64_t end = pos + len;
    while (right->Size() < end) {
        end -= right->Size();
        ++right;
        VERIFY(right!=segms.end());
    }
    if (left==right) {
        return MDBGSeq(
            {EdgeSegment(left->edge, left->start + pos, left->start + end)});
    }

    std::list<EdgeSegment> res;
    res.emplace_back(left->edge, left->start + pos, left->end);
    for (++left; left!=right; ++left) {
        res.emplace_back(*left);
    }
    res.emplace_back(right->edge, right->start, right->start + end);
    MDBGSeq subseq(std::move(res));
    VERIFY(subseq.Size()==len);
    return subseq;
}

[[nodiscard]] double MDBGSeq::Cov() const {
    // std::cout << CovListEdgeSegm(segms) << " " << cov << "\n";
    // VERIFY(std::fabs(CovListEdgeSegm(segms) - cov) < 0.0001);
    // return cov/size;
    double min_cov{std::numeric_limits<double>::max()};
    for (const EdgeSegment &segment : segms) {
        min_cov = std::min(min_cov, segment.Cov());
    }
    return min_cov;
}

[[nodiscard]] bool MDBGSeq::operator==(const MDBGSeq &rhs) const {
    return segms==rhs.segms;
}
