//
// Created by anton on 23.01.2020.
//

#pragma once

#include <cstdlib>
#include <vector>

template<class T>
T CountSum(T const * from, const T * const to) {
    T res(0);
    while(from != to) {
        res += *from;
        from += 1;
    }
    return res;
}

template<class T>
size_t total_size(const std::vector<T> &data) {
    size_t res = 0;
    for(const T & val : data) {
        res += val.size();
    }
    return res;
}

template<class Iterator>
std::vector<size_t> histogram(Iterator begin, Iterator end, size_t max_value, size_t bucket_size = 1) {
    std::vector<size_t> res((max_value + bucket_size - 1)  / bucket_size + 1);
    while(begin != end) {
        res[std::min((*begin + bucket_size - 1) / bucket_size, res.size() - 1)] += 1;
        ++begin;
    }
    return res;
}