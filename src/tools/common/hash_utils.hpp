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
    if(val == 0) {
        os << '0';
        return os;
    }
    while (val != 0) {
        res.push_back(val % 10);
        val /= 10;
    }
    for (auto it = res.rbegin(); it != res.rend(); ++it) {
        os << *it;
    }
    return os;
}

inline hashing::htype stohtype(const std::string &s) {
    hashing::htype val = 0;
    for (char c : s) {
        val = val * 10 + c - '0';
    }
    return val;
}

inline hashing::htype stohtype(const std::string &s, size_t start_pos) {
    hashing::htype val = 0;
    while(start_pos < s.size() && s[start_pos] >= '0' && s[start_pos] <= '9') {
        val = val * 10 + s[start_pos] - '0';
        start_pos++;
    }
    return val;
}

inline std::istream &operator>>(std::istream &is, hashing::htype &val) {
    std::string tmp;
    is >> tmp;
    val = stohtype(tmp);
    return is;
}
