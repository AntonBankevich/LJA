#pragma once

#include <assembly_graph/data_structures/suffix_tracker.hpp>
#include "assembly_graph/random_access_paths.hpp"
#include "dbg/sparse_dbg.hpp"

namespace dbg {
    GraphPath FindReliableExtension(Vertex &start, size_t len, double min_cov);
    std::vector<GraphPath> FindPlausibleBulgeAlternatives(const GraphPath &path, size_t max_diff, double min_cov);
    std::vector<GraphPath> FindPlausibleTipAlternatives(const GraphPath &path, size_t max_diff, double min_cov);
    
    GraphPath FindLongestCoveredForwardExtension(Edge &start, size_t max_len, double min_rel_cov, double max_err_cov);
    GraphPath FindLongestCoveredExtension(Edge &start, size_t max_len, double min_rel_cov, double max_err_cov);

    std::vector<GraphPath> SuffixSupportedBulgeAlternatives(const ag::SuffixTracker<DBGTraits> &tracker, const GraphPath &bulge, size_t threshold);
    EdgeId SuffixSupportedExtension(const ag::SuffixRecord<DBGTraits> &record, const GraphPath &start, size_t min_good, size_t max_bad);
    GraphPath FullSuffixSupportedExtension(const ag::SuffixRecord<DBGTraits> &record, GraphPath start, size_t min_good_cov,
                                        size_t max_bad_cov, size_t max_size = size_t(-1));
    std::vector<GraphPath> SuffixSupportedTipAlternatives(const ag::SuffixTracker<DBGTraits> &record, const GraphPath &tip, double threshold);
}