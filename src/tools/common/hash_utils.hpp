#pragma once
#include <iostream>
#include <vector>

namespace hashing {
    typedef unsigned __int128 htype;

    template<class Key>
    struct alt_hasher {
        size_t operator()(const Key &k) const;
    };

    template<>
    struct alt_hasher<htype> {
        size_t operator()(const htype &x) const {
            return (size_t(x) * 31) ^ size_t(x >> 64u);
        }
    };
}

inline std::ostream &operator<<(std::ostream &os, hashing::htype val) {
    std::vector<size_t> res;
    while (val != 0) {
        res.push_back(val % 10);
        val /= 10;
    }
    for (auto it = res.rbegin(); it != res.rend(); ++it) {
        os << *it;
    }
    return os;
}

inline std::istream &operator>>(std::istream &is, hashing::htype &val) {
    val = 0;
    std::string tmp;
    is >> tmp;
    for (char c : tmp) {
        val = val * 10 + c - '0';
    }
    return is;
}
