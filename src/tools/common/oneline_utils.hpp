//
// Created by anton on 08.07.2020.
//

#pragma once
#include <vector>
#include <functional>
#include <algorithm>

namespace oneline {
    template<class T>
    std::vector<T> Concat(const std::vector<T> &v1, const std::vector<T> &v2) {
        std::vector<T> res = v1;
        res.insert(res.end(), v2.begin(), v2.end());
        return std::move(res);
    }

    template<class U, class V, class I>
    std::vector<V> map(I begin, I end, const std::function<V(U &)> &f) {
        std::vector<V> result;
        std::for_each(begin, end, [&](U& param){ result.push_back(f(param));});
        return std::move(result);
    }

    template<class V, class I>
    std::vector<V> filter(I begin, I end, const std::function<bool(const V&)> &f) {
        std::vector<V> result;
        std::for_each(begin, end, [&](typename I::reference param){if(f(param)) result.emplace_back(std::move(param));});
        return std::move(result);
    }

    template<class V, class C>
    C filter(const C&container, const std::function<bool(const V&)> &f) {
        return std::move(filter<V, typename C::const_iterator>(container.begin(), container.getFinish(), f));
    }

    template<class V, class I>
    std::vector<V> initialize(I begin, const I &end) {
        std::vector<V> result;
        std::for_each(begin, end, [&](const typename I::value_type & param){ result.emplace_back(param);});
        return std::move(result);
    }

    template<typename V, class C>
    std::vector<V> initialize(const C &container) {
        return std::move(initialize<V, typename C::const_iterator>(container.begin(), container.end()));
    }

    template<class V, class U, class I>
    std::vector<V> initialize(I begin, const I &end) {
        std::vector<V> result;
        std::for_each(begin, end, [&](const U & param){ result.emplace_back(param);});
        return std::move(result);
    }

}