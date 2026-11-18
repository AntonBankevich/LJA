#pragma once
#include <cstddef>
#include <vector>

//    Colors the vertices 0..n-1 of an undirected graph, given by adjacency lists (conflicts[i] lists the neighbours of i,
//    each edge in both lists; repetitions are allowed), so that neighbours get different colors. Returns the color
//    classes, each listing its vertices in increasing order. Parallel (Jones-Plassmann with fixed pseudo-random
//    priorities): the result depends only on the graph, not on the number of threads. Greedy, so it uses at most
//    (max degree + 1) colors; the degree must be below 64.
std::vector<std::vector<size_t>> ParallelColorClasses(size_t threads, const std::vector<std::vector<size_t>> &conflicts);
