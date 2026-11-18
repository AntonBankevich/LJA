#include "graph_coloring.hpp"
#include "verify.hpp"
#include <atomic>
#include <cstdint>
#include <numeric>
#include <utility>
#include <omp.h>

namespace {
//    A fixed pseudo-random priority (splitmix64), so that the coloring does not depend on the number of threads.
    uint64_t Priority(size_t index) {
        uint64_t z = index + 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    bool HigherPriority(size_t a, size_t b) {
        return std::make_pair(Priority(a), a) > std::make_pair(Priority(b), b);
    }
}

//    A vertex is colored once all its neighbours of higher priority are, with the smallest color none of them has, so
//    each color is a function of the graph alone, whatever the order of the parallel computation. The uncolored vertex
//    of highest priority is always ready, so every round colors at least one vertex.
std::vector<std::vector<size_t>> ParallelColorClasses(size_t threads, const std::vector<std::vector<size_t>> &conflicts) {
    omp_set_num_threads(threads);
    const size_t n = conflicts.size();
    std::vector<std::atomic<int>> colors(n);
    for (std::atomic<int> &c : colors)
        c.store(-1, std::memory_order_relaxed);
    std::vector<size_t> remaining(n);
    std::iota(remaining.begin(), remaining.end(), 0);
    while (!remaining.empty()) {
        std::vector<char> colored(remaining.size(), 0);
#pragma omp parallel for schedule(dynamic, 100)
        for (size_t r = 0; r < remaining.size(); r++) {
            const size_t i = remaining[r];
            uint64_t used = 0;
            bool ready = true;
            for (size_t j : conflicts[i]) {
                if (!HigherPriority(j, i))
                    continue;
                const int c = colors[j].load(std::memory_order_acquire);
                if (c < 0) {
                    ready = false;
                    break;
                }
                used |= uint64_t(1) << c;
            }
            if (!ready)
                continue;
            VERIFY(~used != 0);
            colors[i].store(__builtin_ctzll(~used), std::memory_order_release);
            colored[r] = 1;
        }
        std::vector<size_t> next;
        for (size_t r = 0; r < remaining.size(); r++)
            if (!colored[r])
                next.push_back(remaining[r]);
        remaining = std::move(next);
    }
    std::vector<std::vector<size_t>> classes;
    for (size_t i = 0; i < n; i++) {
        const size_t color = colors[i].load(std::memory_order_relaxed);
        if (classes.size() <= color)
            classes.resize(color + 1);
        classes[color].push_back(i);
    }
    return classes;
}
