#pragma once
#include "dbg/sparse_dbg.hpp"
#include "dbg/paths.hpp"

std::unordered_map<dbg::Vertex *, size_t> findReachable(dbg::Vertex &start, double min_cov, size_t max_dist);
std::vector<DBGGraphPath> FindPlausibleBulgeAlternatives(const DBGGraphPath &path,
                                                         size_t max_diff, double min_cov);
DBGGraphPath FindReliableExtension(dbg::Vertex &start, size_t len, double min_cov);
std::vector<DBGGraphPath> FindPlausibleTipAlternatives(const DBGGraphPath &path,
                                                       size_t max_diff, double min_cov);

DBGGraphPath FindLongestCoveredForwardExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov);
DBGGraphPath FindLongestCoveredExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov);

