#pragma once

#include <assembly_graph/data_structures/suffix_tracker.hpp>
#include "assembly_graph/random_access_paths.hpp"
#include "dbg/sparse_dbg.hpp"

namespace dbg {
    ag::GraphPath FindReliableExtension(Vertex &start, size_t len, double min_cov);
    std::vector<ag::GraphPath> FindPlausibleBulgeAlternatives(const ag::GraphPath &path, size_t max_diff, double min_cov);
    std::vector<ag::GraphPath> FindPlausibleTipAlternatives(const ag::GraphPath &path, size_t max_diff, double min_cov);

    ag::GraphPath FindLongestCoveredForwardExtension(Edge &start, size_t max_len, double min_rel_cov, double max_err_cov);
    ag::GraphPath FindLongestCoveredExtension(Edge &start, size_t max_len, double min_rel_cov, double max_err_cov);

    std::vector<ag::GraphPath> SuffixSupportedBulgeAlternatives(const ag::SuffixTracker &tracker, const ag::GraphPath &bulge, size_t threshold);
    EdgeId SuffixSupportedExtension(const ag::SuffixRecord &record, const ag::GraphPath &start, size_t min_good, size_t max_bad);
    ag::GraphPath FullSuffixSupportedExtension(const ag::SuffixRecord &record, ag::GraphPath start, size_t min_good_cov,
                                        size_t max_bad_cov, size_t max_size = size_t(-1));
    std::vector<ag::GraphPath> SuffixSupportedTipAlternatives(const ag::SuffixTracker &record, const ag::GraphPath &tip, double threshold);
}