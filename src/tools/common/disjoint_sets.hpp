#pragma once
#include <unordered_map>
#include <vector>
#include <iostream>

template<typename T>
class DisjointSet {
private:
    std::unordered_map<T, T> parent;
    bool cnt = false;
public:
    T get(T obj) {
        auto it = parent.find(obj);
        if(it == parent.end() || obj == it->second) {
            return obj;
        } else {
            parent[obj] = get(parent[obj]);
            return parent[obj];
        }
    }

    void link(T obj1, T obj2) {
        if(obj1 == obj2)
            return;
        obj1 = get(obj1);
        obj2 = get(obj2);
        if(obj1 == obj2)
            return;
        T obj1_group = get(obj1);
        T obj2_group = get(obj2);
        if(cnt) {
            parent[get(obj1)] = get(obj2);
            parent[obj2_group] = obj2_group;
        } else {
            parent[get(obj2)] = get(obj1);
            parent[obj1_group] = obj1_group;
        }
        cnt = !cnt;
    }

    template<class S>
    std::unordered_map<T, std::vector<T>> subsets(S all) {
        std::unordered_map<T, std::vector<T>> res;
        for(T obj : all) {
            res[get(obj)].emplace_back(obj);
        }
        return std::move(res);
    }

    std::unordered_map<T, std::vector<T>> nontrivialSubsets() {
        std::unordered_map<T, std::vector<T>> res;
        for(auto it : parent) {
            res[get(it.first)].emplace_back(it.first);
        }
        return std::move(res);
    }

    typename std::unordered_map<T, T>::iterator begin() {
        return parent.begin();
    }
    typename std::unordered_map<T, T>::iterator end() {
        return parent.end();
    }
};