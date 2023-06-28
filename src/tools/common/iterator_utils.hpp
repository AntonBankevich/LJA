#pragma once

#include "verify.hpp"
#include <functional>

template<class Iterator>
class SkippingIterator {
public:
    typedef typename Iterator::value_type value_type;
    typedef value_type &reference;
    typedef value_type &pointer;
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

    pointer operator->() const {
        return &(operator*());
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
    typedef V &reference;
    typedef V &pointer;
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

    reference operator*() const {
        VERIFY(cur < MAX_SIZE);
        VERIFY(values[cur] != nullptr);
        return *values[cur];
    }

    pointer operator->() const {
        return &(operator*());
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

template<class T>
T& Dereference(T*&pointer) {return *pointer;}

template<class T>
T& ConstDereference(T*const &pointer) {return *pointer;}


template<class Iterator, typename V>
class TransformingIterator {
public:
    typedef V value_type;
    typedef V &reference;
    typedef V *pointer;
private:
    typedef typename Iterator::reference old_value_type;

    Iterator iterator;
    Iterator end;
    std::function<reference(old_value_type)> transform;

public:
    TransformingIterator(Iterator iterator, Iterator end,
                     const std::function<reference(old_value_type)> &transform) :
            iterator(iterator), end(end), transform(transform) {
    }

    static TransformingIterator<Iterator, V> DereferencingIterator(Iterator iterator, Iterator end) {
        return {iterator, end, std::function<V&(V*&)>(&Dereference<V>)};
    }

    static TransformingIterator<Iterator, V> DereferencingConstIterator(Iterator iterator, Iterator end) {
        return {iterator, end, std::function<V&(V*const&)>(&ConstDereference<V>)};
    }

    reference operator*() const {
        return transform(*iterator);
    }

    pointer operator->() const {
        return &(operator*());
    }

    TransformingIterator& operator++() {
        ++iterator;
        return *this;
    }

    TransformingIterator operator++(int) const {
        TransformingIterator other = *this;
        ++other;
        return other;
    }

    bool operator==(const TransformingIterator &other) const {
        return iterator== other.iterator && end == other.end;
    }

    bool operator!=(const TransformingIterator &other) const {
        return !operator==(other);
    }
};

template<class Iterator, typename V>
class TransformingGenerator {
public:
    typedef V value_type;
    typedef V reference;
private:
    typedef typename Iterator::reference old_value_type;

    Iterator iterator;
    Iterator end;
    std::function<reference(old_value_type)> transform;

public:
    TransformingGenerator(Iterator iterator, Iterator end,
                         const std::function<reference(old_value_type)> &transform) :
            iterator(iterator), end(end), transform(transform) {
    }

    reference operator*() const {
        return transform(*iterator);
    }

    TransformingGenerator& operator++() {
        return ++iterator;
    }

    TransformingGenerator operator++(int) const {
        TransformingGenerator other = *this;
        ++other;
        return other;
    }

    bool operator==(const TransformingGenerator &other) const {
        return iterator== other.iterator && end == other.end;
    }

    bool operator!=(const TransformingGenerator &other) const {
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