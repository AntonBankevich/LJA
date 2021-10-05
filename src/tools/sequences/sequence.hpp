#pragma once

#include "common/oneline_utils.hpp"
#include "common/output_utils.hpp"
#include "nucl.hpp"
#include "IntrusiveRefCntPtr.h"
#include "common/verify.hpp"
#include <functional>
#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <sstream>

template<size_t N, size_t base = 2>
struct log_ {
    const static size_t value = 1 + log_<N / base, base>::value;
};

template<size_t base>
struct log_<1, base> {
    const static size_t value = 0;
};

template<size_t base>
struct log_<0, base> {
    const static size_t value = 0;
};

class Sequence {
    // Type to store Seq in Sequences
    typedef u_int64_t ST;
    // Number of bits in ST
    const static size_t STBits = sizeof(ST) << 3u;
    // Number of nucleotides in ST
    const static size_t STN = (STBits >> 1u);
    // Number of bits in STN (for faster div and mod)
    const static size_t STNBits = log_<STN, 2>::value;

    class ManagedNuclBuffer final : public llvm::ThreadSafeRefCountedBase<ManagedNuclBuffer> {
    public:
        explicit ManagedNuclBuffer(size_t nucls) : _data(new ST[Sequence::DataSize(nucls)]) {
        }

        ManagedNuclBuffer(size_t nucls, ST *buf) : _data(new ST[Sequence::DataSize(nucls)]) {
            std::uninitialized_copy(buf, buf + Sequence::DataSize(nucls), data());
        }

    private:
        ST *_data;
    public:
        const ST *data() const { return _data; }

        ST *data() { return _data; }

        ~ManagedNuclBuffer() {
            delete[] _data;
        }

    };

    size_t from_;
    size_t size_;
    bool rtl_; // Right to left + complimentary (?)
    llvm::IntrusiveRefCntPtr<ManagedNuclBuffer> data_;

    static size_t DataSize(size_t size) {
        return (size + STN - 1) >> STNBits;
    }

    template<typename S>
    void InitFromNucls(const S &s, bool rc = false) {
        size_t bytes_size = DataSize(size_);
        ST *bytes = data_->data();
        if(size_ > 0 && (!(is_dignucl(s[0]) || is_nucl(s[0])))) {
            std::cerr << "Bad nucleotide sequence " << size_ << " " << s << std::endl;
        }
        VERIFY(size_ == 0 || is_dignucl(s[0]) || is_nucl(s[0]));

        // Which symbols does our string contain : 0123 or ACGT?
        bool digit_str = size_ == 0 || is_dignucl(s[0]);

        // data -- one temporary variable corresponding to the i-th array element
        // and some counters
        ST data = 0;
        size_t cnt = 0;
        size_t cur = 0;

        if (rc) {
            for (int i = (int) size_ - 1; i >= 0; --i) {
                //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
                char c = complement(digit_str ? s[(unsigned) i] : dignucl(s[(unsigned) i]));

                data = data | (ST(c) << cnt);
                cnt += 2;

                if (cnt == STBits) {
                    bytes[cur++] = data;
                    cnt = 0;
                    data = 0;
                }
            }
        } else {
            for (size_t i = 0; i < size_; ++i) {
                //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
                char c = digit_str ? s[i] : dignucl(s[i]);

                data = data | (ST(c) << cnt);
                cnt += 2;

                if (cnt == STBits) {
                    bytes[cur++] = data;
                    cnt = 0;
                    data = 0;
                }
            }
        }

        if (cnt != 0)
            bytes[cur++] = data;

        for (; cur < bytes_size; ++cur)
            bytes[cur] = 0;
    }

    Sequence(size_t size, int)
            : from_(0), size_(size), rtl_(false), data_(new ManagedNuclBuffer(size_)) {}

    //Low level constructor. Handle with care.
    Sequence(const Sequence &seq, size_t from, size_t size, bool rtl)
            : from_(from), size_(size), rtl_(rtl), data_(seq.data_) {}

public:
    /**
     * Sequence initialization (arbitrary size string)
     *
     * @param s ACGT or 0123-string
     */
    explicit Sequence(const char *s, bool rc = false)
            : Sequence(strlen(s), 0) {
        InitFromNucls(s, rc);
    }

    explicit Sequence(const std::string &s, bool rc = false)
            : Sequence(s.size(), 0) {
        InitFromNucls(s, rc);
    }

    explicit Sequence(const std::vector<char> &s, bool rc = false)
            : Sequence(s.size(), 0) {
        InitFromNucls(s, rc);
    }

    explicit Sequence(const std::vector<unsigned char> &s, bool rc = false)
            : Sequence(s.size(), 0) {
        InitFromNucls(s, rc);
    }

    explicit Sequence(char *s, bool rc = false)
            : Sequence(strlen(s), 0) {
        InitFromNucls(s, rc);
    }

    Sequence()
            : Sequence(size_t(0), 0) {
        memset(data_->data(), 0, DataSize(size_));
    }

    Sequence(const Sequence &s)
            : Sequence(s, s.from_, s.size_, s.rtl_) {}

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
        return Sequence(str());
    }

    unsigned char operator[](const size_t index) const {
        VERIFY(index < size_);
        const ST *bytes = data_->data();
        if (rtl_) {
            size_t i = from_ + size_ - 1 - index;
            return complement((bytes[i >> STNBits] >> ((i & (STN - 1u)) << 1u)) & 3u);
        } else {
            size_t i = from_ + index;
            return (bytes[i >> STNBits] >> ((i & (STN - 1u)) << 1u)) & 3u;
        }
    }

    size_t asNumber() const {
        size_t res = 0;
        const ST *bytes = data_->data();
        if (rtl_) {
            for(size_t i = from_ + size_ - 1; i + 1 >= from_ + 1; i++) {
                res = (res << 2u) + (complement((bytes[i >> STNBits] >> ((i & (STN - 1u)) << 1u)) & 3u));
            }
        } else {
            for(size_t i = from_; i < from_ + size_; i++) {
                res = (res<< 2u) + ((bytes[i >> STNBits] >> ((i & (STN - 1u)) << 1u)) & 3u);
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

    bool operator<=(const Sequence &other) const {
        return !(other < *this);
    }

    bool operator!=(const Sequence &that) const {
        return !(operator==(that));
    }

    /**
     * @param from inclusive
     * @param to exclusive;
     */
    inline Sequence Subseq(size_t from, size_t to) const;

    inline Sequence Subseq(size_t from) const; // up to size_ by default

    inline Sequence operator+(const Sequence &s) const;

    inline Sequence operator*(size_t mult) const;

    inline Sequence Prefix(size_t count) const;

    inline Sequence Suffix(size_t count) const;

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

    inline std::string err() const;

    size_t size() const {
        return size_;
    }

    bool empty() const {
        return size() == 0;
    }

    bool startsWith(const Sequence & other) const {
        return (other.size() <= size()) && (Subseq(0, other.size()) == other);
    }

    bool endsWith(const Sequence & other) const {
        return (other.size() <= size()) && (Subseq(size() - other.size(), size()) == other);
    }

    bool nonContradicts(const Sequence & other) const {
        size_t ms = std::min(size(), other.size());
        return Subseq(0, ms) == other.Subseq(0, ms);
    }

    template<class Seq>
    bool contains(const Seq &s, size_t offset = 0) const {
        VERIFY_DEV(offset + s.size() <= size());

        for (size_t i = 0, e = s.size(); i != e; ++i)
            if (operator[](offset + i) != s[i])
                return false;

        return true;
    }

    Sequence operator!() const {
        return Sequence(*this, from_, size_, !rtl_);
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
 * start of Sequence is Seq with preferred size
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


/**
 * @todo optimize sequence copy
 */
Sequence Sequence::operator+(const Sequence &s) const {
    if (data_ == s.data_ && rtl_ == s.rtl_ &&
            (
                (!rtl_ && this->from_ + size_ == s.from_ ) ||
                (rtl_ && this->from_ == s.from_ + s.size_)
            )
        )
    {
        return Sequence(*this, std::min(from_, s.from_), size_ + s.size_, rtl_);
    } else {
        return Sequence(str() + s.str());
    }
}

std::string Sequence::str() const {
    VERIFY(size_ < 1000000000000ull);
    std::string res(size_, '-');
    for (size_t i = 0; i < size_; ++i) {
        res[i] = nucl(this->operator[](i));
    }
    return res;
}

std::string Sequence::err() const {
    std::ostringstream oss;
    oss << "{ *data=" << data_->data() <<
        ", from_=" << from_ <<
        ", size_=" << size_ <<
        ", rtl_=" << int(rtl_) << " }";
    return oss.str();
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
        VERIFY(index < size_);
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
