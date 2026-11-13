#include <common/verify.hpp>
#include <vector>
#include <common/oneline_utils.hpp>
#include "nucl_buffer.hpp"
size_t NuclDeck::initial = 16;
size_t NuclDeck::mult = 2;
//    Out-of-class definition: std::min binds it by reference.
const size_t NuclBuffer::STN;

//    Positions in [left, right) satisfy zero - capacity <= pos < zero + capacity, and zero <= right.
size_t NuclDeck::physical(int pos) const {
    return pos >= zero ? size_t(pos - zero) : data.size() + pos - zero;
}

NuclBuffer::ST NuclDeck::readBits(int pos, size_t cnt) const {
    size_t p = physical(pos);
    size_t avail = data.size() - p;
    if (cnt <= avail)
        return data.getBits(p, cnt);
    return data.getBits(p, avail) | (data.getBits(0, cnt - avail) << (avail << 1u));
}

void NuclDeck::writeBits(int pos, size_t cnt, NuclBuffer::ST bits) {
    size_t p = physical(pos);
    size_t avail = data.size() - p;
    if (cnt <= avail) {
        data.setBits(p, cnt, bits);
    } else {
        data.setBits(p, avail, bits);
        data.setBits(0, cnt - avail, bits >> (avail << 1u));
    }
}

void NuclDeck::relinearize(size_t new_capacity) {
    VERIFY(new_capacity >= size());
    NuclBuffer new_data(std::max<size_t>(new_capacity, 1));
    size_t n = size();
    for (size_t i = 0; i < n; i += NuclBuffer::STN) {
        size_t cnt = std::min(NuclBuffer::STN, n - i);
        new_data.setBits(i, cnt, readBits(left + int(i), cnt));
    }
    data = std::move(new_data);
    zero = left;
}

void NuclDeck::reserveBack(size_t n) {
    if (size() + n > data.size()) {
        relinearize(std::max(data.size() * mult, size() + n));
    } else if (right + int(n) > zero + int(data.size())) {
//        The new positions wrap around: they go to the beginning of the buffer, so left > zero and the old content
//        keeps its buffer indices relative to the moved zero.
        zero += int(data.size());
    }
}

NuclDeck::NuclDeck(NuclDeck &&other) noexcept : data(std::move(other.data)), left(other.left), right(other.right),
        zero(other.zero) {
    other.left = other.right = other.zero = 0;
}

NuclDeck &NuclDeck::operator=(NuclDeck &&other) noexcept {
    data = std::move(other.data);
    left = other.left;
    right = other.right;
    zero = other.zero;
    other.left = other.right = other.zero = 0;
    return *this;
}

NuclDeck::NuclDeck(Iterator begin, Iterator end) {
    VERIFY(&begin.getDeck() == &end.getDeck() && begin <= end);
    size_t n = end - begin;
    if (n > 0) {
        data = NuclBuffer(n);
        for (size_t i = 0; i < n; i += NuclBuffer::STN) {
            size_t cnt = std::min(NuclBuffer::STN, n - i);
            data.setBits(i, cnt, begin.getDeck().readBits(begin.getCur() + int(i), cnt));
        }
        right = int(n);
    }
}

void NuclDeck::reserveFront(size_t n) {
    if (size() + n > data.size())
        relinearize(std::max(data.size() * mult, size() + n));
}

void NuclDeck::push_back(const NuclDeck &s) {
    if (&s == this) {
        NuclDeck copy(s);
        push_back(copy);
        return;
    }
    size_t n = s.size();
    reserveBack(n);
    for (size_t i = 0; i < n; i += NuclBuffer::STN) {
        size_t cnt = std::min(NuclBuffer::STN, n - i);
        writeBits(right + int(i), cnt, s.readBits(s.left + int(i), cnt));
    }
    right += int(n);
}

void NuclDeck::push_front(const NuclDeck &s) {
    if (&s == this) {
        NuclDeck copy(s);
        push_front(copy);
        return;
    }
    size_t n = s.size();
    reserveFront(n);
    for (size_t i = 0; i < n; i += NuclBuffer::STN) {
        size_t cnt = std::min(NuclBuffer::STN, n - i);
        writeBits(left - int(n) + int(i), cnt, s.readBits(s.left + int(i), cnt));
    }
    left -= int(n);
}

void NuclDeck::push_back(unsigned char c) {
    if(right == left + data.size()) {
        resize();
    }
//    reserveBack may have moved zero ahead of right already; physical handles both cases.
    if(right == zero + int(data.size()))
        zero = right;
    data.setNucl(physical(right), c);
    right++;
}

void NuclDeck::pop_back() {
    VERIFY(size() > 0);
    if(right == zero)
        zero -= int(data.size());
    right--;
}

void NuclDeck::pop_back(size_t val) {
    right -= int(val);
    if(right < zero)
        zero -= int(data.size());
}

void NuclDeck::push_front(unsigned char c) {
    if(right == left + data.size()) {
        resize();
    }
    left--;
    data.setNucl(getIndex(left), c);
}

void NuclDeck::pop_front() {
    VERIFY(size() > 0);
    left++;
}

void NuclDeck::resize() {
    relinearize(std::max(data.size() * mult, initial));
}

size_t NuclDeck::getIndex(int pos) const {
    VERIFY_MSG(pos >= left && pos < right, "Out of bounds: " << pos << " " << left << " " << right << " " << zero);
    if(pos >= zero)
        return pos - zero;
    else
        return data.size() + pos - zero;
}

NuclBuffer::ST NuclBuffer::getBits(size_t pos, size_t cnt) const {
    size_t w = pos >> STNBits;
    size_t off = (pos & (STN - 1u)) << 1u;
    ST res = _data[w] >> off;
    if (off != 0 && off + (cnt << 1u) > STBits)
        res |= _data[w + 1] << (STBits - off);
    return cnt == STN ? res : res & ((ST(1) << (cnt << 1u)) - 1u);
}

void NuclBuffer::setBits(size_t pos, size_t cnt, ST bits) {
    size_t w = pos >> STNBits;
    size_t off = (pos & (STN - 1u)) << 1u;
    ST mask = cnt == STN ? ~ST(0) : (ST(1) << (cnt << 1u)) - 1u;
    bits &= mask;
    _data[w] = (_data[w] & ~(mask << off)) | (bits << off);
    if (off != 0 && off + (cnt << 1u) > STBits) {
        size_t sh = STBits - off;
        _data[w + 1] = (_data[w + 1] & ~(mask >> sh)) | (bits >> sh);
    }
}

void NuclBuffer::setNucl(size_t pos, unsigned char value) {
    size_t x = dignucl(value)^getNucl(pos);
    size_t sh = x << ((pos & (STN - 1u)) << 1u);
    size_t old = _data[pos >> STNBits];
    _data[pos >> STNBits] = old^sh;
    VERIFY(getNucl(pos) == value);
}

std::string NuclDeck::str() const {
    std::stringstream ss;
    for(unsigned char c : *this)
        ss << nucl(char(c));
    return ss.str();
}

bool NuclDeck::operator==(const NuclDeck &other) const {
    if(size() != other.size()) return false;
    for (size_t i = 0; i < size(); i += NuclBuffer::STN) {
        size_t cnt = std::min(NuclBuffer::STN, size() - i);
        if (readBits(left + int(i), cnt) != other.readBits(other.left + int(i), cnt))
            return false;
    }
    return true;
}