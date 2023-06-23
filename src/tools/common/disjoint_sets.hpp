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
        if(it == parent.end()) {
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
        if(cnt)
            parent[get(obj1)] = get(obj2);
        else
            parent[get(obj2)] = get(obj1);
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
        for(auto it : res) {
            it.second.emplace_back(it.first);
        }
        return std::move(res);
    }
};