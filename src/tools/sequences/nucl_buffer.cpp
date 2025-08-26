#include <common/verify.hpp>
#include <vector>
#include <common/oneline_utils.hpp>
#include "nucl_buffer.hpp"
size_t NuclDeck::initial = 16;
size_t NuclDeck::mult = 2;

void NuclDeck::push_back(const NuclDeck &s) {
    auto e = s.end();
    for(auto b = s.begin(); b != e; ++b)
        push_back((unsigned char) dignucl(char(*b)));
}

void NuclDeck::push_front(const NuclDeck &s) {
    auto b = s.begin();
    for(auto e = s.end(); e != b;) {
        --e;
        push_front((unsigned char)dignucl(char(*e)));
    }
}

void NuclDeck::push_back(unsigned char c) {
    if(right == left + data.size()) {
        resize();
    }
    size_t pos = right - zero;
    if(pos == data.size()) {
        zero = right;
        pos = 0;
    }
    right++;
    data.setNucl(pos, c);
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
    std::vector<char> tmp = oneline::initialize<char>(begin(), end());
    NuclBuffer new_data(data.size() * mult);
    if(left < zero) {
        for (int i = left; i < zero; i++)
            new_data.setNucl(new_data.size() + i - zero, data.getNucl(data.size() + i - zero));
        for (int i = zero; i < right; i++)
            new_data.setNucl(i - zero, data.getNucl(i - zero));
    } else {
        for (int i = left; i < right; i++)
            new_data.setNucl(i - zero, data.getNucl(i - zero));
    }
    data = std::move(new_data);
    std::vector<char> tmp1 = oneline::initialize<char>(begin(), end());
    VERIFY(tmp == tmp1);
}

size_t NuclDeck::getIndex(int pos) const {
    VERIFY_MSG(pos >= left && pos < right, "Out of bounds: " << pos << " " << left << " " << right << " " << zero);
    if(pos >= zero)
        return pos - zero;
    else
        return data.size() + pos - zero;
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
    for(size_t i = 0; i < size(); i++)
        if((*this)[i] != other[i])
            return false;
    return true;
}