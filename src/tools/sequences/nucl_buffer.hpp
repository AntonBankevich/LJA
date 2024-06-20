#pragma once
#include "nucl.hpp"
#include "IntrusiveRefCntPtr.h"
#include <cstring>
#include <memory>
#include "sstream"


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

class NuclBuffer {
public:
    typedef u_int64_t ST;
// Number of bits in ST
    const static size_t STBits = sizeof(ST) << 3u;
// Number of nucleotides in ST
    const static size_t STN = (STBits >> 1u);
// Number of bits in STN (for faster div and mod)
    const static size_t STNBits = log_<STN, 2>::value;

    static size_t DataSize(size_t size) {return (size + STN - 1) >> STNBits;}

    NuclBuffer(const NuclBuffer &other) : _data(new ST[DataSize(other._size)]), _size(other._size) {
        memcpy(_data, other._data, (DataSize(other._size) << STNBits) >> 2u);
    }
    NuclBuffer(NuclBuffer &&other)  noexcept : _data(other._data), _size(other._size) {other._data = nullptr;}
    NuclBuffer &operator=(NuclBuffer &&other)  noexcept {std::swap(_data, other._data); std::swap(_size, other._size); return *this;};
    NuclBuffer &operator=(const NuclBuffer &other) {*this = NuclBuffer(other); return *this;}

    NuclBuffer() : _data(nullptr), _size(0) {};
    explicit NuclBuffer(size_t nucls) : _data(new ST[DataSize(nucls)]), _size(nucls) {
    }
    NuclBuffer(size_t nucls, ST *buf) : _data(new ST[DataSize(nucls)]), _size(nucls) {std::uninitialized_copy(buf, buf + DataSize(nucls), _data);
    }


    void InitZero() {memset(_data, 0, size());}
    template<typename S>
    void InitFromNucls(const S &s, size_t sz, bool rc = false);
    template<typename S>
    void InitFromNucls(const S &s, bool rc = false) {InitFromNucls(s, size(), rc);}
    template<typename T>
    void InitFromNucls(T begin, T end, size_t sz, bool rc);


private:
    ST *_data;
    size_t _size;
public:
    unsigned char getNucl(size_t pos) const {return (_data[pos >> STNBits] >> ((pos & (STN - 1u)) << 1u)) & 3u;}
    void setNucl(size_t pos, unsigned char value);
    const ST *data() const { return _data; }
    size_t size() const {return _size;}
    virtual ~NuclBuffer() {delete[] _data;}
};

template<typename S>
void NuclBuffer::InitFromNucls(const S &s, size_t sz, bool rc) {
    VERIFY(sz <= size());
    size_t bytes_size = DataSize(size());

    // data -- one temporary variable corresponding to the i-th array element
    // and some counters
    ST data = 0;
    size_t cnt = 0;
    size_t cur = 0;

    if (rc) {
        for (int i = (int) sz - 1; i >= 0; --i) {
            VERIFY_MSG(is_dignucl(s[i]) || is_nucl(s[i]), "Non-ACGTacgt symbols in the input sequences");
            char c = complement(dignucl(s[(unsigned) i]));

            data = data | (ST(c) << cnt);
            cnt += 2;

            if (cnt == STBits) {
                _data[cur++] = data;
                cnt = 0;
                data = 0;
            }
        }
    } else {
        for (size_t i = 0; i < sz; ++i) {
            VERIFY_MSG(is_dignucl(s[i]) || is_nucl(s[i]), "Non-ACGTacgt symbols in the input sequences");
            char c = dignucl(s[i]);

            data = data | (ST(c) << cnt);
            cnt += 2;

            if (cnt == STBits) {
                _data[cur++] = data;
                cnt = 0;
                data = 0;
            }
        }
    }

    if (cnt != 0)
        _data[cur++] = data;

    for (; cur < bytes_size; ++cur)
        _data[cur] = 0;
}

template<typename T>
void NuclBuffer::InitFromNucls(T begin, T end, size_t sz, bool rc) {
    VERIFY(sz <= size());
    size_t bytes_size = DataSize(size());

    // data -- one temporary variable corresponding to the i-th array element
    // and some counters
    ST data = 0;
    size_t cnt = 0;
    size_t cur = 0;

    if (rc) {
        for (int i = (int) sz - 1; i >= 0; --i) {
            --end;
            VERIFY_MSG(is_dignucl(*end) || is_nucl(*end), "Non-ACGTacgt symbols in the input sequences");
            char c = complement(dignucl(*end));

            data = data | (ST(c) << cnt);
            cnt += 2;

            if (cnt == STBits) {
                _data[cur++] = data;
                cnt = 0;
                data = 0;
            }
        }
    } else {
        for (size_t i = 0; i < sz; ++i) {
            VERIFY_MSG(is_dignucl(*begin) || is_nucl(*begin), "Non-ACGTacgt symbols in the input sequences");
            char c = dignucl(*begin);

            data = data | (ST(c) << cnt);
            cnt += 2;

            if (cnt == STBits) {
                _data[cur++] = data;
                cnt = 0;
                data = 0;
            }
            ++begin;
        }
    }

    if (cnt != 0)
        _data[cur++] = data;

    for (; cur < bytes_size; ++cur)
        _data[cur] = 0;
}

class NuclDeck {
public:
    class Iterator;
    friend class Iterator;
private:
    NuclBuffer data;
    int left = 0;
    int right = 0;
    int zero = 0;
    static size_t initial;
    static size_t mult;

    size_t getIndex(int pos) const;

//    TODO: remake it with memcpy. Make borders always good.
    void resize();
public:
    class IteratorInfo {
    protected:
        int cur;
    public:
        IteratorInfo(int cur) : cur(cur) {} // NOLINT(google-explicit-constructor)
        int getCur() const {return cur;}
    };
    class Iterator : public IteratorInfo {
    private:
        const NuclDeck *deck;
    public:
        typedef char value_type;
        typedef char reference;
        Iterator(const NuclDeck &deck, int cur) : IteratorInfo(cur), deck(&deck) {}
        Iterator(const NuclDeck &deck, IteratorInfo cur) : IteratorInfo(cur), deck(&deck) {}
        unsigned char operator*() const {return deck->data.getNucl(deck->getIndex(cur));}
        Iterator &operator++() {cur++; return *this;}
        Iterator operator++(int) const {return {*deck, cur + 1};} // NOLINT(cert-dcl21-cpp)
        Iterator &operator--() {cur--; return *this;}
        Iterator operator--(int) const {return {*deck, cur - 1};} // NOLINT(cert-dcl21-cpp)
        int operator-(const Iterator &other) const {return cur - other.cur;}
        Iterator &operator +=(int val) {cur += val; return *this;}
        Iterator &operator -=(int val) {cur -= val; return *this;}
        Iterator operator +(int val) {Iterator res = *this; res += val; return res;}
        Iterator operator -(int val) {Iterator res = *this; res -= val; return res;}
        bool operator==(const Iterator &other) const {return deck == other.deck && cur == other.cur;}
        bool operator!=(const Iterator &other) const {return deck != other.deck || cur != other.cur;}
        bool operator<=(const Iterator &other) const {return cur <= other.cur;}
        bool operator<(const Iterator &other) const {return cur < other.cur;}
        bool operator>=(const Iterator &other) const {return cur >= other.cur;}
        bool operator>(const Iterator &other) const {return cur > other.cur;}
        size_t getPos() const {return cur - deck->left;}
        const NuclDeck &getDeck() const {return *deck;}
    };
    explicit NuclDeck(size_t size = initial) : data(size) {VERIFY(size > 0);}
    template<class I>
    NuclDeck(I begin, I end);
    template<class S>
    explicit NuclDeck(const S&s) : data(s.size()) {push_back(s);}
    NuclDeck(const NuclDeck &) = default;
    NuclDeck(NuclDeck &&) = default;
    NuclDeck &operator=(const NuclDeck &other) = default;
    NuclDeck &operator=(NuclDeck &&) = default;
    size_t size() const {return right - left;}
    void push_back(unsigned char c);
    void push_back(const NuclDeck &s);
    template<class S>
    void push_back(const S &s);
    void pop_back();
    void pop_back(size_t val);
    void push_front(unsigned char c);
    void push_front(const NuclDeck &s);
    template<class S>
    void push_front(const S &s);
    void pop_front();
    void pop_front(size_t val) {left += int(val);}
    template<class S>
    void replaceFront(size_t sz, const S &s);
    template<class S>
    void replaceBack(size_t sz, const S &s);
    Iterator begin() const {return {*this, left};}
    Iterator end() const {return {*this, right};}
    unsigned char operator[](size_t ind) const {return data.getNucl(this->getIndex(left + int(ind)));}
    bool empty() const {return left == right;}
    std::string str() const;
    bool operator==(const NuclDeck &other) const;

    bool startsWith(const NuclDeck &other) const;
    bool endsWith(const NuclDeck &other) const;
};

template<class I>
NuclDeck::NuclDeck(I begin, I end) : NuclDeck() {
    while(begin != end) {
        push_back(*begin);
        ++begin;
    }
}

template<class S>
void NuclDeck::push_back(const S &s) {
    for(size_t i = 0; i < s.size(); i++)
        push_back((unsigned char)(dignucl(s[i])));
}

template<class S>
void NuclDeck::push_front(const S &s) {
    for(size_t i = s.size(); i > 0; i--)
        push_front(s[i - 1]);
}

template<class S>
void NuclDeck::replaceFront(size_t sz, const S &s) {
    VERIFY(sz <= size());
    for(size_t i = 0; i < sz && i < s.size(); i++) {
        data.setNucl(getIndex(left + int(sz - 1 - i)), s[s.size() - 1 - i]);
    }
//        If replacement is shorter than initial, we adjust left and shorten storage
    left = std::max<int>(left, left + sz - s.size());
//        If adjustment is longer than initial, we add the remaining information to front and extend storage
    for(size_t i = sz; i < s.size(); i++) {
        push_front(s[s.size() - 1 - i]);
    }
}

template<class S>
void NuclDeck::replaceBack(size_t sz, const S &s) {
    VERIFY(sz <= size());
    for(size_t i = 0; i < sz && i < s.size(); i++) {
        data.setNucl(getIndex(right - sz + i), s[i]);
    }
//        If replacement is shorter than initial, we adjust right and shorten storage
    right = std::min<int>(right, right + s.size() - sz);
//        If adjustment is longer than initial, we add the remaining information to back and extend storage
    for(size_t i = sz; i < s.size(); i++) {
        push_back(s[i]);
    }
}

inline std::ostream &operator<<(std::ostream &os, const NuclDeck &nucls) {
    for(unsigned char c : nucls) {
        os << nucl(char(c));
    }
    return os;
}
