#pragma once

#include "verify.hpp"
#include <functional>

template<class Iterator>
class SkippingIterator {
public:
    typedef typename Iterator::value_type value_type;
private:
    Iterator iterator;
    Iterator end;
    std::function<bool(const value_type &)> use;

    void seek() {
        while(iterator != end && !use(*iterator)) {
            ++iterator;
        }
    }
public:
    SkippingIterator(Iterator iterator, Iterator end, const std::function<bool(const value_type &)> &use) :
                                        iterator(iterator), end(end), use(use) {
        seek();
    }

    value_type &operator*() const {
        return *iterator;
    }

    SkippingIterator& operator++() {
        ++iterator;
        seek();
        return *this;
    }

    SkippingIterator operator++(int) const {
        SkippingIterator other = *this;
        ++other;
        return other;
    }

    bool operator==(const SkippingIterator &other) const {
        return iterator== other.iterator && end == other.end;
    }

    bool operator!=(const SkippingIterator &other) const {
        return !operator==(other);
    }
};

template<class Iterator, typename V, size_t MAX_SIZE = 8>
class ApplyingIterator {
public:
    typedef V value_type;
private:
    typedef typename Iterator::reference old_value_type;

    Iterator iterator;
    Iterator end;
    std::function<std::array<value_type*, MAX_SIZE>(old_value_type)> apply;
    std::array<value_type*, MAX_SIZE> values;
    size_t cur;

    void fill() {
        if(iterator == end){
            values = {};
        } else {
            values = apply(*iterator);
        }
    }

    void seek() {
        while(iterator != end && (cur == MAX_SIZE || values[cur] == nullptr)) {
            cur = 0;
            ++iterator;
            fill();
        }
    }

public:
    ApplyingIterator(Iterator iterator, Iterator end,
                     const std::function<std::array<value_type*, MAX_SIZE>(old_value_type)> &apply) :
                        iterator(iterator), end(end), apply(apply), values(), cur(0) {
        fill();
        seek();
    }

    value_type &operator*() const {
        VERIFY(cur < MAX_SIZE);
        VERIFY(values[cur] != nullptr);
        return *values[cur];
    }

    ApplyingIterator& operator++() {
        cur++;
        seek();
        return *this;
    }

    ApplyingIterator operator++(int) const {
        ApplyingIterator other = *this;
        ++other;
        return other;
    }

    bool operator==(const ApplyingIterator &other) const {
        return iterator== other.iterator && end == other.end && cur == other.cur;
    }

    bool operator!=(const ApplyingIterator &other) const {
        return !operator==(other);
    }
};

template<class Iterator>
class IterableStorage {
private:
    Iterator _begin;
    Iterator _end;
public:
    IterableStorage(const Iterator &_begin, const Iterator &_end) : _begin(_begin), _end(_end) {
    }

    Iterator begin() const {
        return _begin;
    }

    Iterator end() const {
        return _end;
    }
};