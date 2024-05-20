#pragma once

#include "verify.hpp"
#include <functional>
#include <vector>

template<class Iterator>
class SkippingIterator {
public:
    typedef typename Iterator::value_type value_type;
    typedef typename Iterator::reference reference;
    typedef value_type &pointer;
private:
    Iterator iterator;
    Iterator end;
    std::function<bool(reference)> use;

    void seek() {
        while(iterator != end && !use(*iterator)) {
            ++iterator;
        }
    }
public:
    SkippingIterator(Iterator iterator, Iterator end, const std::function<bool(reference)> &use) :
                                        iterator(iterator), end(end), use(use) {
        seek();
    }

    reference operator*() const {
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

template<typename I>
class CountingIterator {
public:
    typedef I value_type;
    typedef I reference;
    typedef I pointer;
private:
    I pos;

public:
    CountingIterator(I pos) : pos(pos){
    }

    size_t operator*() const {
        return pos;
    }

    CountingIterator& operator++() {
        pos++;
        return *this;
    }

    CountingIterator operator++(int) const {
        CountingIterator other = *this;
        ++other;
        return other;
    }

    bool operator==(const CountingIterator &other) const {
        return pos == other.pos;
    }

    bool operator!=(const CountingIterator &other) const {
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

template<class T, class U = T*>
T& Dereference(U&pointer) {return *pointer;}

template<class T, class U = T*>
T& ConstDereference(U const &pointer) {return *pointer;}


template<class Iterator, typename V>
class TransformingIterator {
public:
    typedef V value_type;
    typedef V &reference;
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
        return {iterator, end, std::function<V&(old_value_type &)>(&Dereference<V, old_value_type>)};
    }

    static TransformingIterator<Iterator, V> DereferencingConstIterator(Iterator iterator, Iterator end) {
        return {iterator, end, std::function<V&(old_value_type const&)>(&ConstDereference<V, old_value_type>)};
    }

    reference operator*() const {
        old_value_type v = *iterator;
        return transform(v);
    }

    V* operator->() const {
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

    TransformingIterator& operator--() {
        --iterator;
        return *this;
    }

    TransformingIterator operator--(int) const {
        TransformingIterator other = *this;
        --other;
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
class Generator {
public:
    typedef V value_type;
    typedef V reference;
private:
    typedef typename Iterator::reference old_value_type;

    Iterator iterator;
    Iterator end;
    std::function<bool(old_value_type)> use;
    std::function<reference(old_value_type)> generate;

    void seek() {
        while(iterator != end && !use(*iterator)) {
            ++iterator;
        }
    }

public:
    Generator(Iterator iterator, Iterator end, const std::function<reference(old_value_type)> &generate, 
                      const std::function<bool(old_value_type)> &use = [](old_value_type){return true;}) :
            iterator(iterator), end(end), use(use), generate(generate) {
        seek();
    }

    reference operator*() const {
        return generate(*iterator);
    }

    Generator& operator++() {
        ++iterator;
        seek();
        return *this;
    }

    Generator operator++(int) const {
        Generator other = *this;
        ++other;
        return other;
    }

    bool operator==(const Generator &other) const {
        return iterator== other.iterator && end == other.end;
    }

    bool operator!=(const Generator &other) const {
        return !operator==(other);
    }
};


template<class Iterator>
class IterableStorage {
private:
    Iterator _begin;
    Iterator _end;
public:
    typedef Iterator iterator;
    IterableStorage(const Iterator &_begin, const Iterator &_end) : _begin(_begin), _end(_end) {
    }

    Iterator begin() const {
        return _begin;
    }

    Iterator end() const {
        return _end;
    }
};

template<class Iterator>
class ComplexIterator {
private:
    typedef typename std::vector<IterableStorage<Iterator>>::const_iterator storage_iterator;
    storage_iterator cur_storage;
    storage_iterator end_storage;
    Iterator *cur_item;
    void seek() {
        VERIFY(cur_storage == end_storage || cur_item != nullptr);
        while(cur_storage != end_storage && *cur_item == cur_storage->end()) {
            ++cur_storage;
            if(cur_storage != end_storage) {
                *cur_item = cur_storage->begin();
            } else {
                delete cur_item;
                cur_item = nullptr;
            }
        }
    }
public:
    typedef typename Iterator::reference reference;
    typedef typename Iterator::value_type value_type;

    ComplexIterator(storage_iterator cur_storage, storage_iterator end_storage, Iterator cur_item) : cur_storage(
            cur_storage), end_storage(end_storage), cur_item(new Iterator(cur_item)) {
        seek();
    }

    ComplexIterator(storage_iterator end_storage) : cur_storage(end_storage), end_storage(end_storage), cur_item(nullptr) {
    }
    ComplexIterator(const ComplexIterator<Iterator> &other) : cur_storage(other.cur_storage), end_storage(other.end_storage),
                    cur_item(new Iterator(*other.cur_item)) {}

    ComplexIterator(ComplexIterator<Iterator> &&other) : cur_storage(std::move(other.cur_storage)), end_storage(std::move(other.end_storage)),
                                                              cur_item(new Iterator(*other.cur_item)) {
        other.cur_item = nullptr;
    }

    ~ComplexIterator() {
        delete cur_item;
        cur_item = nullptr;
    }

    reference operator*() const {
        VERIFY(cur_storage != end_storage);
        VERIFY(cur_item != nullptr);
        return **cur_item;
    }

    ComplexIterator& operator++() {
        VERIFY(cur_item != nullptr);
        Iterator &item = *cur_item;
        ++item;
        seek();
        return *this;
    }

    ComplexIterator operator++(int) const {
        ComplexIterator other = *this;
        ++other;
        return other;
    }

    bool operator==(const ComplexIterator &other) const {
        return cur_storage == other.cur_storage && (cur_item == other.cur_item ||
                    (cur_item != nullptr && other.cur_item != nullptr && *cur_item == *other.cur_item));
    }

    bool operator!=(const ComplexIterator &other) const {
        return !operator==(other);
    }

};

template<class Iterator>
class ComplexIterableStorage {
private:
    std::vector<IterableStorage<Iterator>> storages;
public:
    ComplexIterableStorage() {}

    bool empty() const {
        return storages.empty();
    }
    void operator+=(IterableStorage<Iterator> storage) & {
        storages.template emplace_back(std::move(storage));
    }
    ComplexIterator<Iterator> begin() const & {
        if(storages.empty()) {
            return {storages.end()};
        } else {
            return {storages.begin(), storages.end(), storages.begin()->begin()};
        }
    }
    ComplexIterator<Iterator> end() const & {
        return {storages.end()};
    }
};