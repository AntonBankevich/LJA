#pragma once
#include "deque"
#include <cstddef>
//TODO: add namespace
template<class T>
class MinQueue {
    std::deque<std::pair<T, std::size_t>> q;
    size_t max_pos;
public:
    MinQueue() = default;

    bool empty() const {return q.empty();}
    T get() const {return q.front().first;}
    size_t size() const {return q.size();}

    void push(const T &item) {
        while (!q.empty() && q.back().first > item) {
            q.pop_back();
        }
        q.emplace_back(item, max_pos);
        max_pos++;
    }
    void pop(size_t pos) {
        if (!q.empty() && q.front().second < pos) {
            q.pop_front();
        }
    }
};