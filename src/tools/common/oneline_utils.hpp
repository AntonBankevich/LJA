//
// Created by anton on 08.07.2020.
//

#pragma once
#include <vector>
#include <functional>
#include <algorithm>

namespace oneline {
    template<class U, class V, class I>
    std::vector<V> map(I begin, I end, std::function<V(U &)> f) {
        std::vector<V> result;
        std::for_each(begin, end, [&](U& param){ result.push_back(f(param));});
        return std::move(result);
    }

    template<class V, class I>
    std::vector<V> filter(I begin, I end, std::function<bool(V&)> f) {
        std::vector<V> result;
        std::for_each(begin, end, [&](V& param){if(f(param)) result.emplace_back(std::move(param));});
        return std::move(result);
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