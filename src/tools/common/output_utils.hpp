//
// Created by anton on 23.01.2020.
//

#pragma once
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

template<class U, class V>
std::ostream& operator<<(std::ostream& out, const std::pair<U, V>& item) {
    return out << "(" << item.first << ", " << item.second << ")";
}

//inline std::ostream& operator<<(std::ostream& out, const unsigned __int128& item) {
//    std::vector<char> res;
//    unsigned __int128 tmp = item;
//    while(tmp != 0) {
//        res.push_back(char((tmp % 10) + '0'));
//        tmp /= 10;
//    }
//    return out << std::string(res.rbegin(), res.rend());
//}


template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& tree) {
    if(tree.size() == 0) {
        return out << "[]" << std::endl;
    }
    out << "[";
    for(size_t i = 0; i + 1 < tree.size(); i += 1) {
        out << tree[i] << ", ";
    }
    return out << tree[tree.size() - 1] << "]";
}