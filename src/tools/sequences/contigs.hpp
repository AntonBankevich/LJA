#pragma once

#include "sequence.hpp"
#include "nucl.hpp"
#include "IntrusiveRefCntPtr.h"
#include "common/string_utils.hpp"
#include "common/verify.hpp"
#include <algorithm>
#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <sstream>
#include <unordered_map>

using std::string;
using std::string;
namespace basic {
    inline string Reverse(const string &s){
        if (s[0] == '-')
            return s.substr(1);
        else
            return "-" + s;
    }
}

class RawSequence {
public:
    size_t id;
    const std::string seq;
    RawSequence(size_t _id, std::string _seq): id(_id), seq(std::move(_seq)) {
    }
};

//TODO: store cut_left and cut_right instead of left and right and merge with IdSegment
template<class T>
class Segment{
    T *contig_ptr;
public:
    size_t left;
    size_t right;
    Segment(T &contig_, size_t left_, size_t right_) : left(left_), right(right_), contig_ptr(&contig_){
        VERIFY(0 <= left and left <= right and right <= contig_ptr->truncSize())
    }

    Segment(T &contig) : contig_ptr(&contig), left(0), right(contig.truncSize()) {
        VERIFY(0 <= left and left <= right and right <= contig_ptr->truncSize())
    }

    Segment() : contig_ptr(nullptr), left(left), right(right) {}

    bool valid() const {return contig_ptr == nullptr;}
    T &contig() const {return *contig_ptr;}
    size_t size() const {return right - left;}
    size_t cutLeft() const {return left;}
    size_t cutRight() const {return contig_ptr->truncSize() - right;}

    Sequence truncSeq() const {return contig_ptr->truncSeq().Subseq(left, right);}
    Sequence fullSeq() const {
        size_t k = contig_ptr->getStartSize();
        if(left >= k) {
            return contig_ptr->truncSeq().Subseq(left - k, right);
        } else {
            return contig_ptr->getSeq().Subseq(left, right + k);
        }
    }

    size_t dist(const Segment<T> &other) const {
        VERIFY(contig_ptr == other.contig_ptr);
        if (right <= other.left)
            return other.left - right;
        else if (other.right <= left)
            return left - other.right;
        else
            return 0;
    }

    Segment<T> RC() const {
        return Segment(contig_ptr->rc(), contig_ptr->truncSize() - right, contig_ptr->truncSize() - left);
    }

    bool inter(const Segment &other) const {
        return contig_ptr == other.contig_ptr and not (right <= other.left or left >= other.right);
    }

    int interSize(const Segment &other) const {
        if (not inter(other))
            return -1;
        else
            return int(std::min(right, other.right) - std::max(left, other.left));
    }

    bool operator<(const Segment<T> &other) const {
        return contig() < other.contig() ||
                    (contig() == other.contig() &&
                            (left < other.left || left == other.left && right < other.right));
    }

    bool operator>(const Segment<T> &other) const {
        return contig() > other.contig() ||
               (contig() == other.contig() &&
                (left > other.left || left == other.left && right > other.right));
    }

    Segment<T> operator+(const Segment<T> &other) const {
        VERIFY(contig_ptr == other.contig_ptr);
        VERIFY(right == other.left);
        return {*contig_ptr, left, other.right};
    }

    bool operator==(const Segment<T> &other) const {
        return contig() == other.contig() && left == other.left && right == other.right;
    }

    bool operator!=(const Segment<T> &other) const {
        return contig() != other.contig() || left != other.left || right != other.right;
    }

    Segment<T> shrinkRightBy(size_t len) const {
        VERIFY(len <= size());
        return {*contig_ptr, left, right - len};
    }

    Segment<T> shrinkLeftBy(size_t len) const {
        VERIFY(len <= size());
        return {*contig_ptr, left + len, right};
    }

    Segment<T> shrinkRightToLen(size_t len) const {
        return {*contig_ptr, left, left + len};
    }

    Segment<T> shrinkLeftToLen(size_t len) const {
        return {*contig_ptr, right - len, right};
    }

    Segment<T> extendBy(size_t len) const {
        return {*contig_ptr, left - std::min(left, len), std::min(contig_ptr->truncSize(), right + len)};
    }

    Segment<T> extendRight(size_t len) const {
        return {*contig_ptr, left, std::min(contig_ptr->truncSize(), right + len)};
    }

    Segment<T> extendLeft(size_t len) const {
        return {*contig_ptr, left - std::min(left, len), right};
    }

    Segment<T> unite(const Segment<T> &other) const {
        return {*contig_ptr, std::min(left, other.left), std::max(right, other.right)};
    }

    Segment<T> nest(const Segment<T> &other) const {return {other.contig(), other.left + left, other.left + right};}

    std::string coordinaresStr() const {
        std::stringstream ss;
        ss << "[" << left << ":";
        if (right > contig().truncSize() * 3 / 4)
            ss << contig().truncSize() << "-" << (contig().truncSize() - right);
        else
            ss << right;
        ss << "]";
        return ss.str();
    }
};

template<class T>
class Position{
    T *contig_ptr;
    size_t pos;
public:
    Position() : contig_ptr(nullptr), pos(0){}

    Position(T &contig, size_t pos) : contig_ptr(&contig), pos(pos){
        VERIFY(0 <= pos && pos <= contig_ptr->fullSize())
    }

    T &contig() const {return *contig_ptr;}
    size_t getPos() const {return pos;}
    Sequence Suffix() const {return contig_ptr->getSeq().Subseq(pos);}
    Sequence Prefix() const {return contig_ptr->getSeq().Subseq(0, pos);}

    size_t dist(Position<T> &other) const {
        VERIFY(contig_ptr == other.contig_ptr);
        return pos < other.pos ? other.pos - pos : pos - other.pos;
    }

    size_t dist(Segment<T> other) const {
        VERIFY(contig_ptr == other.contig_ptr);
        if(pos < other.left) return other.left - pos;
        if(pos > other.right) return pos - other.right;
        return 0;
    }

    Position<T> RC() const {return {contig_ptr->rc(), contig_ptr->fullSize() - pos};}

    bool operator<(const Position<T> &other) const {
        return contig() < other.contig() ||
               (contig() == other.contig() && (pos < other.pos));
    }

    bool operator>(const Position<T> &other) const {
        return contig() > other.contig() ||
               (contig() == other.contig() && (pos > other.pos));
    }

    bool operator==(const Position<T> &other) const {
        return contig() == other.contig() && pos == other.pos;
    }

    bool operator!=(const Position<T> &other) const {
        return contig() != other.contig() || pos != other.pos;
    }
};

namespace std {
    template<class T>
    struct hash<Position<T>>{
    size_t operator()(const Position<T> &value) const noexcept {
        typedef typename T::id_type id_type;
        return std::hash<id_type>()(value.contig().getInnerId()) + 7 * value.getPos();
    }
};
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Segment<T>& seg) {
    return os << seg.contig().getInnerId() << seg.coordinaresStr();
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Position<T>& seg) {
    return os << seg.contig().getInnerId() << "[" << seg.getPos() << "]";
}

template<class T, typename I>
class NamedSequence {
private:
    Sequence seq;
    I id;
public:
    typedef I id_type;

    NamedSequence(Sequence _seq, I _id) : seq(std::move(_seq)), id(std::move(_id)) {
    }

    const I &getInnerId() const {
        return id;
    }

    const Sequence &getSeq() const {
        return seq;
    }

    const Sequence &truncSeq() const {
        return seq;
    }

    size_t getStartSize() const {return 0;}

    Segment<T> asSegment() const {
        return Segment<T>(*this, 0u, truncSize());
    }

    Segment<T> segment(size_t left, size_t right) const {
        return Segment<T>(*(static_cast<const T*>(this)), left, right);
    }

    Segment<T> suffix(int pos) const {
        if (pos < 0)
            pos = truncSize() + pos;
        if (pos < 0)
            pos = 0;
        if (pos > truncSize())
            pos = truncSize();
        return Segment<T>(*this, pos, truncSize());
    }

    Segment<T> prefix(size_t len) const {
        len = min(len, truncSize());
        return Segment<T>(*this, 0, len);
    }

    size_t truncSize() const {
        return getSeq().size();
    }

    size_t fullSize() const {
        return getSeq().size();
    }

//    NamedSequence &rc() const {
//        return *_rc;
//    }

    bool operator==(const NamedSequence &other) const {
        return getInnerId() == other.getInnerId();
    }

    char operator[](size_t ind) const {
        return getSeq()[ind];
    }

    string str() const {
        return getSeq().str();
    }

    bool isNull() const {
        return getSeq().empty();
    }
};

class Contig: public NamedSequence<Contig, std::string> {
public:
    Contig(): NamedSequence(Sequence(), ""){
    }

    Contig(Sequence _seq, string _id): NamedSequence(std::move(_seq), std::move(_id)) {
    }

//    Contig(const Sequence &_seq, const string &_id, Contig *_rc): NamedSequence(_seq, _id, _rc) {
//    }

    Contig(const string &_seq, string _id): NamedSequence(Sequence(_seq), std::move(_id)) {
    }

    Contig RC() const {
        if(getInnerId()[0] == '-')
            return {!getSeq(), getInnerId().substr(1, getInnerId().size() - 1)};
        else
            return {!getSeq(), "-" + getInnerId()};
    }

    bool operator<(const Contig &other) {
        return getInnerId() < other.getInnerId();
    }
//    Contig(const string &_seq, const string &_id, Contig *_rc): NamedSequence(Sequence(_seq), _id, _rc) {
//    }
};

inline std::ostream& operator<<(std::ostream& os, const Contig& contig) {
    os << contig.getInnerId();
    return os;
}

class StringContig {
private:
    static std::string extractId(const std::string &s) {
        if(s.find(' ') == size_t(-1))
            return s;
        else
            return s.substr(0, s.find(' '));
    }

    static std::string extractComment(const std::string &s) {
        if(s.find(' ') == size_t(-1))
            return "";
        else
            return s.substr(s.find(' ') + 1);
    }

    static std::string makeUpperCase(std::string &&s) {
        for(size_t i = 0; i < s.size(); i++) {
            if('a' <= s[i] && s[i] <= 'z')
                s[i] += 'A' - 'a';
        }
        return std::move(s);
    }
public:
    std::string id;
    std::string comment;
    std::string seq;
    static bool homopolymer_compressing;
    static size_t min_dimer_to_compress;
    static size_t max_dimer_size;
    static size_t dimer_step;

    static void SetDimerParameters(const std::string &s) {
        std::vector<std::string> vals = split(s, ",");
        StringContig::min_dimer_to_compress = std::stoull(vals[0]);
        StringContig::max_dimer_size = std::stoull(vals[1]);
        StringContig::dimer_step = std::stoull(vals[2]);
    }

    StringContig() : id(""), comment(""), seq("") {
    }

    StringContig(std::string _seq, std::string _id) : id(extractId(_id)), comment(extractComment(_id)), seq(makeUpperCase(std::move(_seq))) {
    }

    StringContig(StringContig && other) = default;

    StringContig(const StringContig & other) = default;

    StringContig& operator=(StringContig && other) = default;

    void compress();

//    void atCompress() {
//        size_t cur = 0;
//        size_t at_len = 0;
//        bool start = true;
//        for(size_t i = 0; i < seq.size(); i++) {
//            if(i >= 2) {
//                if (seq[i] == seq[cur - 2])
//                    at_len += 1;
//                else
//                    start = false;
//            }
//        }
//        seq.erase(std::unique(seq.begin(), seq.end()), seq.end());
//    }

    Contig makeContig() {
        compress();
        return Contig(Sequence(seq), id);
    }

//    Contig makeCompressedContig() {
//        compress();
//        return makeContig();
//    }

    Sequence makeSequence() {
        compress();
        return Sequence(seq);
    }

//    Sequence makeCompressedSequence() {
//        compress();
//        return makeSequence();
//    }

    bool isNull() const {
        return id.empty() && seq.empty();
    }

    size_t size() const {
        return seq.size();
    }
};