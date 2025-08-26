#pragma once

#include "common/oneline_utils.hpp"
#include "common/output_utils.hpp"
#include "common/string_utils.hpp"
#include "nucl_buffer.hpp"
#include "nucl.hpp"
#include "IntrusiveRefCntPtr.h"
#include "common/verify.hpp"
#include <functional>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
//TODO: rewrite everything with bitset to improve binary operations
class ManagedNuclBuffer final : public NuclBuffer, public llvm::ThreadSafeRefCountedBase<ManagedNuclBuffer> {
public:
    explicit ManagedNuclBuffer(size_t nucls) : NuclBuffer(nucls) {}
    ManagedNuclBuffer(size_t nucls, ST *buf) : NuclBuffer(nucls, buf) {}
};

class Sequence {
    // Type to store Seq in Sequences

    size_t from_;
    size_t size_;
    bool rtl_; // Right to left + complimentary (?)
    llvm::IntrusiveRefCntPtr<ManagedNuclBuffer> data_;

    Sequence(size_t size, int)
            : from_(0), size_(size), rtl_(false), data_(new ManagedNuclBuffer(size_)) {}

    //Low level constructor. Handle with care.
    Sequence(const Sequence &seq, size_t from, size_t size, bool rtl)
            : from_(from), size_(size), rtl_(rtl), data_(seq.data_) {
        VERIFY_MSG(from + size <= data_->size(), from << " " << size << " " << data_->size());
        VERIFY_MSG(from <= data_->size(), from << " " << size << " " << data_->size());
        VERIFY_MSG(size <= data_->size(), from << " " << size << " " << data_->size());
    }

    Sequence maxFreeExtension() const {
        if(rtl_)
            return {*this, 0, from_ + size(), rtl_};
        else
            return {*this, from_, data_->size() - from_, rtl_};

    }
public:
    /**
     * Sequence initialization (arbitrary size string)
     *
     * @param s ACGT or 0123-string
     */
    explicit Sequence(const char *s, bool rc = false)
            : Sequence(strlen(s), 0) {
        data_->InitFromNucls(s, rc);
    }

    explicit Sequence(const std::string &s, bool rc = false)
            : Sequence(s.size(), 0) {
        data_->InitFromNucls(s, rc);
    }

    explicit Sequence(const std::vector<char> &s, bool rc = false)
            : Sequence(s.size(), 0) {
        data_->InitFromNucls(s, rc);
    }

    explicit Sequence(char c) : Sequence(std::vector<char>({c})){
    }

    explicit Sequence(const std::vector<unsigned char> &s, bool rc = false)
            : Sequence(s.size(), 0) {
        data_->InitFromNucls(s, rc);
    }

    explicit Sequence(char *s, bool rc = false)
            : Sequence(strlen(s), 0) {
        data_->InitFromNucls(s, rc);
    }

    template<class S>
    explicit Sequence(const S &s, bool rc = false) : Sequence(s.size(), 0) {
        data_->InitFromNucls(s, rc);
    }

    template<class I>
    explicit Sequence(I begin, I end, bool rc = false) : Sequence(size_t(end - begin), int(0)) {
        data_->InitFromNucls(begin, end, end - begin, rc);
    }

    Sequence()
            : Sequence(size_t(0), 0) {
        data_->InitZero();
    }

    Sequence(const Sequence &s)
            : Sequence(s, s.from_, s.size_, s.rtl_) {}

    Sequence(Sequence && s) noexcept
            : from_(s.from_), size_(s.size_), rtl_(s.rtl_), data_(std::move(s.data_)) {}

    static Sequence Concat(const std::vector<Sequence> &v) {
        std::stringstream ss;
        for(const auto &seq : v) {
            ss << seq.str();
        }
        return Sequence(ss.str());
    }

    Sequence &operator=(const Sequence &rhs) {
        if (&rhs == this)
            return *this;

        from_ = rhs.from_;
        size_ = rhs.size_;
        rtl_ = rhs.rtl_;
        data_ = rhs.data_;

        return *this;
    }

    Sequence &operator=(Sequence &&other) = default;

    Sequence copy() const {
        Sequence res = Sequence(size_, 0);
        res.data_->InitFromNucls(*this, this->rtl_);
        return std::move(res);
;    }

    unsigned char operator[](const size_t index) const {
        VERIFY_MSG(index < size_, itos(index) + " " + itos(size_));
        if (rtl_) {
            size_t i = from_ + size_ - 1 - index;
            return complement(data_->getNucl(i));
        } else {
            size_t i = from_ + index;
            return data_->getNucl(i);
        }
    }

    size_t asNumber() const {
        size_t res = 0;
        if (rtl_) {
            for(size_t i = from_ + size_ - 1; i + 1 >= from_ + 1; i++) {
                res = (res << 2u) + (complement(data_->getNucl(i)));
            }
        } else {
            for(size_t i = from_; i < from_ + size_; i++) {
                res = (res << 2u) + data_->getNucl(i);
            }
        }
        return res;
    }


    bool operator==(const Sequence &that) const {
        if (size_ != that.size_)
            return false;

        if (data_ == that.data_ && from_ == that.from_ && rtl_ == that.rtl_)
            return true;

        for (size_t i = 0; i < size_; ++i) {
            if (this->operator[](i) != that[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator<(const Sequence &other) const {
        for (size_t i = 0; i < size_; ++i) {
            if (i == other.size())
                return true;
            else if (this->operator[](i) != other[i]) {
                return this->operator[](i) < other[i];
            }
        }
        return false;
    }

    bool operator>(const Sequence &other) const {
        return other < *this;
    }

    bool operator<=(const Sequence &other) const {
        return !(other < *this);
    }

    bool operator>=(const Sequence &other) const {
        return !(other > *this);
    }

    bool operator!=(const Sequence &that) const {
        return !(operator==(that));
    }

    bool isCanonical() const {
        return *this <= this->rc();
    }

    /**
     * @param from inclusive
     * @param to exclusive;
     */
    inline Sequence Subseq(size_t from, size_t to) const;

    inline Sequence Subseq(size_t from) const; // up to size_ by default

    Sequence operator+(const Sequence &s) const;

    inline Sequence operator*(size_t mult) const;

    inline Sequence Prefix(size_t count) const;

    inline Sequence Suffix(size_t count) const;

    inline unsigned char lastNucl() const {return operator[](size() - 1);}
    inline unsigned char firstNucl() const {return operator[](0);}

    Sequence dicompress() const {
        if(size() <= 5)
            return *this;
        std::vector<unsigned char> res = {operator[](0), operator[](1),
                                          operator[](2), operator[](3), operator[](4)};
        for(size_t i = 5; i < size(); i++) {
            unsigned char next = operator[](i);
            if(res.size() >= 5 && next == res[res.size() - 2]
                        && res[res.size() - 1] == res[res.size() - 3]
                        && res[res.size() - 2] == res[res.size() - 4]
                        && res[res.size() - 3] == res[res.size() - 5]){
                res.pop_back();
            } else {
                res.emplace_back(next);
            }
            VERIFY(res.back() == next);
        }
        return Sequence(res);
    }

    inline std::string str() const;

    size_t size() const {
        return size_;
    }

    bool empty() const {
        return size() == 0;
    }

    template<class S>
    bool subseqMatch(const S &other, size_t this_start, size_t other_start, size_t len) const {
        if(size() < this_start + len || other.size() < other_start + len)
            return false;
        for(size_t i = 0; i < len; i++)
            if(this->operator[](this_start + i) != other[other_start + i])
                return false;
        return true;
    }

    template<class S>
    bool startsWith(const S & other) const {
        return subseqMatch(other, 0, 0, other.size());
    }

    bool containsAtPosition(const Sequence & other, size_t pos) const {
        return subseqMatch(other, pos, 0, other.size_);
    }

    bool endsWith(const Sequence & other) const {
        return subseqMatch(other, size() - other.size_, 0, other.size_);
    }

    bool nonContradicts(const Sequence & other) const {
        return subseqMatch(other, 0, 0, std::min(size(), other.size_));
    }

    template<class Seq>
    bool contains(const Seq &s, size_t offset = 0) const {
        VERIFY_DEV(offset + s.size() <= size());

        for (size_t i = 0, e = s.size(); i != e; ++i)
            if (operator[](offset + i) != s[i])
                return false;

        return true;
    }

    Sequence rc() const {
        return Sequence(*this, from_, size_, !rtl_);
    }

    Sequence operator!() const {
        return rc();
    }

    size_t commonPrefix(const Sequence & other) const {
        size_t res = 0;
        while(res < size() && res < other.size() && this->operator[](res) == other[res])
            res += 1;
        return res;
    }

    Sequence makeSequence() {
        return *this;
    }
};

inline std::ostream &operator<<(std::ostream &os, const Sequence &s);

/**
 * getStart of Sequence is Seq with preferred size
 */

// O(1)
//including from, excluding to
//safe if not #DEFINE NDEBUG
Sequence Sequence::Subseq(size_t from, size_t to) const {
    VERIFY(from <= to);
    if (rtl_) {
        return Sequence(*this, from_ + size_ - to, to - from, true);
    } else {
        return Sequence(*this, from_ + from, to - from, false);
    }
}

//including from, excluding to
Sequence Sequence::Subseq(size_t from) const {
    return Subseq(from, size_);
}

Sequence Sequence::Prefix(size_t count) const {
    return Subseq(0, count);
}

Sequence Sequence::Suffix(size_t count) const {
    return Subseq(size_ - count);
}



std::string Sequence::str() const {
    VERIFY(size_ < 1000000000000ull);
    std::string res(size_, '-');
    for (size_t i = 0; i < size_; ++i) {
        res[i] = nucl(this->operator[](i));
    }
    return res;
}

std::ostream &operator<<(std::ostream &os, const Sequence &s) {
    os << s.str();
    return os;
}

class SequenceBuilder {
    std::vector<char> buf_;
public:
    template<typename S>
    SequenceBuilder &append(const S &s) {
        for (size_t i = 0; i < s.size(); ++i) {
            buf_.push_back(s[i]);
        }
        return *this;
    }

    template<typename S>
    SequenceBuilder &appendAll(S begin, S end) {
        while(begin != end) {
            append(*begin);
            ++begin;
        }
        return *this;
    }

    SequenceBuilder &append(char c) {
        buf_.push_back(c);
        return *this;
    }

    Sequence BuildSequence() {
        return Sequence(buf_);
    }

    void reserve(size_t n) {
        buf_.reserve(n);
    }

    size_t size() const {
        return buf_.size();
    }

    void clear() {
        return buf_.clear();
    }

    unsigned char operator[](const size_t index) const {
        return buf_[index];
    }

    std::string str() const {
        std::string s(buf_.size(), '-');
        for (size_t i = 0; i < s.size(); ++i) {
            s[i] = nucl(buf_[i]);
        }
        return s;
    }
};

Sequence Sequence::operator*(size_t mult) const {
    SequenceBuilder sb;
    for(size_t i = 0; i < mult; i++) {
        sb.append(*this);
    }
    return sb.BuildSequence();
}


class CompositeSequence {
private:
    std::vector<Sequence> sequences_;
    size_t left_;
    size_t right_;
    size_t size_;
public:
    CompositeSequence(std::vector<Sequence> sequences, size_t from, size_t to) :
            sequences_(std::move(sequences)), left_(from), right_(0) {
        size_ = 0;
        for(Sequence & sequence : sequences_) {
            size_ += sequence.size();
        }
        VERIFY(left_ + right_ <= size_);
        size_ -= left_ + right_;
        right_ = size_ - to;
    }

    explicit CompositeSequence(std::vector<Sequence> sequences) :
            sequences_(std::move(sequences)), left_(0), right_(0) {
        size_ = 0;
        for(Sequence & sequence : sequences_) {
            size_ += sequence.size();
        }
        VERIFY(left_ + right_ <= size_);
        size_ -= left_ + right_;
    }

    size_t size() const {
        return size_;
    }

    unsigned char operator[](size_t index) const {
        VERIFY_MSG(index < size_, itos(index) + " " + itos(size_));
        index += left_;
        size_t cur = 0;
        while(cur < sequences_.size() && index >= sequences_[cur].size()) {
            index -= sequences_[cur].size();
            cur += 1;
        }
        return sequences_[cur][index];
    }

    CompositeSequence operator!() const {
        std::function<Sequence(const Sequence &)> rc = [this](const Sequence &seq) {
            return !seq;
        };
        return {oneline::map(sequences_.rbegin(), sequences_.rend(), rc), right_, left_};
    }

//    TODO reduce sequence vector size
    CompositeSequence Subseq(size_t from, size_t to) const {
        size_t left = from + left_;
        size_t right = size_ - to + right_;
        return {sequences_, left, right};
    }
};
