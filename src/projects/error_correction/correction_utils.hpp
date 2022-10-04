#pragma once
#include "dbg/sparse_dbg.hpp"
#include "dbg/paths.hpp"

std::unordered_map<dbg::Vertex *, size_t> findReachable(dbg::Vertex &start, double min_cov, size_t max_dist);
std::vector<dbg::GraphAlignment> FindPlausibleBulgeAlternatives(const dbg::GraphAlignment &path,
                                                                       size_t max_diff, double min_cov);
dbg::GraphAlignment FindReliableExtension(dbg::Vertex &start, size_t len, double min_cov);
std::vector<dbg::GraphAlignment> FindPlausibleTipAlternatives(const dbg::GraphAlignment &path,
                                                                size_t max_diff, double min_cov);

dbg::GraphAlignment FindLongestCoveredForwardExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov);
dbg::GraphAlignment FindLongestCoveredExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov);

