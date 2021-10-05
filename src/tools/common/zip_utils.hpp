//
// Created by anton on 9/1/20.
//
#pragma once

template<class V1, class V2>
class RecordPair {
public:
    V1 first;
    V2 second;
    RecordPair(V1 v1, V2 v2) : first(std::move(v1)), second(std::move(v2)) {
    }

    size_t size() const {
        return first.size() + second.size();
    }
};


template<class I1, class I2>
class Zip {
public:
    typedef RecordPair<typename I1::value_type, typename I2::value_type> value_type;
    I1 iterator1;
    I2 iterator2;

    Zip(I1 it1, I2 it2) : iterator1(it1), iterator2(it2){}

    value_type operator*() {
        return {std::move(*iterator1), std::move(*iterator2)};
    }

    bool operator==(const Zip & other) const {
        return iterator1 == other.iterator1 || iterator2 == other.iterator2;
    }

    bool operator!=(const Zip & other) const {
        return !(*this == other);
    }

    void operator++() {
        ++iterator1;
        ++iterator2;
    }
};