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

template<class T>
class Segment{
    T *contig_ptr;
public:
    size_t left;
    size_t right;
    Segment(T &contig_, size_t left_, size_t right_) : left(left_), right(right_), contig_ptr(&contig_){
        VERIFY(0 <= left and left <= right and right <= contig_ptr->size())
    }

    T &contig() const {
        return *contig_ptr;
    }

    size_t size() const {
        return right - left;
    }

    Sequence seq() const {
        return contig_ptr->seq.Subseq(left, right);
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
        return Segment(contig_ptr->rc(), contig_ptr->size() - right, contig_ptr->size() - left);
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

    Segment<T> shrinkRight(size_t len) const {
        return {*contig_ptr, left, right - len};
    }

    Segment<T> shrinkLeft(size_t len) const {
        return {*contig_ptr, left + len, right};
    }

    Segment<T> unite(const Segment<T> &other) const {
        return {*contig_ptr, std::min(left, other.left), std::max(right, other.right)};
    }

    std::string coordinaresStr() const {
        std::stringstream ss;
        ss << "[" << left << ":";
        if (right > contig().size() * 3 / 4)
            ss << contig().size() << "-" << (contig().size() - right);
        else
            ss << right;
        ss << "]";
        return ss.str();
    }
};

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Segment<T>& seg) {
    return os << seg.contig().getId() << seg.coordinaresStr();
}

template<class T>
class NamedSequence {
public:
    string id;
    Sequence seq;
//protected:
//    T * _rc;
public:
//    NamedSequence(const Sequence &_seq, string _id, T *_rc) : seq(_seq), id(std::move(_id)), _rc(_rc){
//    }

    NamedSequence(const Sequence &_seq, string _id) : seq(_seq), id(std::move(_id)){
//        _rc = new T(!seq, basic::Reverse(id), static_cast<T*>(this));
    }

    std::string getId() const {
        return id;
    }

    Segment<T> asSegment() const {
        return Segment<T>(*this, 0u, size());
    }

    Segment<T> segment(size_t left, size_t right) const {
        return Segment<T>(*(static_cast<const T*>(this)), left, right);
    }

    Segment<T> suffix(size_t pos) const {
        if (pos < 0)
            pos = size() + pos;
        if (pos < 0)
            pos = 0;
        if (pos > size())
            pos = size();
        return Segment<T>(*this, pos, size());
    }

    Segment<T> prefix(size_t len) const {
        len = min(len, size());
        return Segment<T>(*this, 0, len);
    }

    size_t size() const {
        return seq.size();
    }

//    NamedSequence &rc() const {
//        return *_rc;
//    }

    bool operator==(const NamedSequence &other) const {
        return id == other.id;
    }

    char operator[](size_t ind) const {
        return nucl(seq[ind]);
    }

    string str() const {
        return seq.str();
    }

    bool isNull() const {
        return seq.empty();
    }
};

class Contig: public NamedSequence<Contig> {
public:
    Contig(): NamedSequence(Sequence(), ""){
    }

    Contig(const Sequence &_seq, const string &_id): NamedSequence(_seq, _id) {
    }

//    Contig(const Sequence &_seq, const string &_id, Contig *_rc): NamedSequence(_seq, _id, _rc) {
//    }

    Contig(const string &_seq, const string &_id): NamedSequence(Sequence(_seq), _id) {
    }

    Contig RC() const {
        if(id[0] == '-')
            return Contig(!seq, id.substr(1, id.size() - 1));
        else
            return Contig(!seq, "-" + id);
    }

    bool operator<(const Contig &other) {
        return getId() < other.getId();
    }
//    Contig(const string &_seq, const string &_id, Contig *_rc): NamedSequence(Sequence(_seq), _id, _rc) {
//    }
};

inline std::ostream& operator<<(std::ostream& os, const Contig& contig) {
    os << contig.getId();
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

    StringContig(std::string && _seq, std::string &&_id) : id(extractId(_id)), comment(extractComment(_id)), seq(_seq) {
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



template <class T>
class SequenceCollection {
private:
    std::unordered_map<std::string, T *> items;
public:
    SequenceCollection() {
    }

    SequenceCollection(const std::vector<T*> & sequences) {
        for(T *item: sequences){
            items[item->id] = item;
        }
    }
};