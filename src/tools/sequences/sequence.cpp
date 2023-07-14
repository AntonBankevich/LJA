//
// Created by anton on 9/22/20.
//

#include "sequence.hpp"
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

