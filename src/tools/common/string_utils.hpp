//
// Created by anton on 02.04.2020.
//

#pragma once
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include "verify.hpp"

inline std::string itos(size_t val, size_t min_size = 0) {
    std::stringstream ss;
    ss << val;
    std::string res = ss.str();
    while(res.size() < min_size) {
        res = "0" + res;
    }
    return res;
}

inline std::string itos(int val, size_t min_size = 0) {
    std::stringstream ss;
    ss << val;
    std::string res = ss.str();
    while(res.size() < min_size) {
        res = "0" + res;
    }
    return res;
}

template<class T>
inline T Parse(const std::string &s);

template<>
inline int Parse<int>(const std::string &s) {return std::stoi(s);}

template<>
inline std::string Parse<std::string>(const std::string &s) {return s;}

inline int ParseInt(const std::string &s, size_t pos = 0) {
    bool neg = false;
    if(pos < s.size() && s[pos] == '-') {
        neg = true;
        pos++;
    }
    VERIFY_MSG(pos < s.size() && s[pos] >= '0' && s[pos] <= '9', itos(pos) + " " + itos(s.size()) + " " + (pos < s.size() ? s[pos] : ' '));
    int res = 0;
    while(pos < s.size() && s[pos] >= '0' && s[pos] <= '9') {
        res = res * 10 + s[pos] - '0';
        pos++;
    }
    return neg ? -res : res;
}

static bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

static bool startsWith(const std::string& str, const std::string& prefix)
{
    return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
}

static inline void ltrim_inplace(std::string &s) {
    s.erase(s.begin(), std::find_if_not(s.begin(), s.end(),
                                    [](auto c){return std::isspace(c);}));
}

static inline void rtrim_inplace(std::string &s) {
    s.erase(std::find_if_not(s.rbegin(), s.rend(),
                         [](auto c){return std::isspace(c);}).base(), s.end());
}

static inline std::string trim(std::string s) {
    ltrim_inplace(s);
    rtrim_inplace(s);
    return s;
}

static inline std::string & compress_inplace(std::string &s) {
    s.erase(std::unique(s.begin(), s.end()), s.end());
    return s;
}

static inline std::string mask(const std::string &s, const std::string &pattern = "/\\\"", char value = '_') {
    std::string res = s;
    for(char &c : res) {
        if(pattern.find(c) != size_t(-1))
            c = value;
    }
    return std::move(res);
}

inline std::string join(const std::string &s, const std::vector<std::string> &arr) {
    if(arr.empty())
        return "";
    std::stringstream ss;
    ss << arr[0];
    for(size_t i = 1; i < arr.size(); i++) {
        ss << s << arr[i];
    }
    return ss.str();
}

template<class I>
inline std::string join(const std::string &s, I begin, I end) {
    if(begin == end)
        return "";
    std::stringstream ss;
    ss << *begin;
    ++begin;
    while(begin != end) {
        ss << s << *begin;
        ++begin;
    }
    return ss.str();
}


inline std::vector<std::string> split(const std::string &s) {
    std::vector<std::string> res;
    size_t cur = 0;
    std::string bad = " \n\t";
    while(cur < s.size()) {
        size_t next = cur;
        while(next < s.size() && bad.find(s[next]) == size_t(-1)) {
//            std::cout << s[cur] << " " << size_t(s[next]) << " " << size_t('\t') << std::endl;
            next += 1;
        }
        if (next > cur) {
            res.push_back(s.substr(cur, next - cur));
        }
        cur = next + 1;
    }
    if(res.empty())
        res.emplace_back("");
    return res;
}

inline std::vector<std::string> split(const std::string &s, const std::string &delim) {
    std::vector<std::string> res;
    size_t cur = 0;
    while(cur < s.size()) {
        size_t next = cur;
        while(next < s.size() && delim.find(s[next]) == size_t(-1)) {
            next += 1;
        }
        if (next > cur) {
            res.push_back(s.substr(cur, next - cur));
        }
        cur = next + 1;
    }
    return res;
}