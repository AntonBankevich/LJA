//
// Created by anton on 9/22/20.
//

#include "sequence.hpp"

ManagedNuclBuffer *ManagedNuclBuffer::Preallocated(ST *words, size_t nucls) {
    auto *res = new ManagedNuclBuffer();
    res->usePreallocated(words, nucls);
    return res;
}

void ManagedNuclBuffer::usePreallocated(ST *words, size_t nucls) {
    const int users = UseCount();
    VERIFY(users <= 1);
    if (owns)
        freeWords();
    setWords(words, nucls);
    owns = false;
}

ManagedNuclBuffer::~ManagedNuclBuffer() {
    if (!owns)
        setWords(nullptr, 0);
}

Sequence Sequence::compressed(const CompressionParameters &parameters) const {
    Sequence res = copy();
    struct Nucls {
        ManagedNuclBuffer &buffer;
        unsigned char get(size_t i) const {return buffer.getNucl(i);}
        void set(size_t i, unsigned char c) {buffer.setNucl(i, c);}
    } nucls{*res.data_};
    return res.Subseq(0, CompressInPlace(nucls, res.size(), parameters));
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
        return {*this, std::min(from_, s.from_), size_ + s.size_, rtl_};
    }
    if(size() + s.size_ > 100) {
        if (this->maxFreeExtension().Subseq(size()).startsWith(s)) {
            return this->maxFreeExtension().Subseq(0, size() + s.size_);
        }
        if ((!s).maxFreeExtension().Subseq(s.size_).startsWith(!*this)) {
            return !((!s).maxFreeExtension().Subseq(0, size_ + s.size_));
        }
    }
    return Sequence(str() + s.str());
}

