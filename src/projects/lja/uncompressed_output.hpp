#pragma once

#include <common/logging.hpp>
#include "multi_graph.hpp"

std::vector<Contig> printUncompressedResults(logging::Logger &logger, size_t threads, multigraph::MultiGraph &graph,
                              const std::vector<Contig> &uncompressed, const std::experimental::filesystem::path &out_dir, bool debug);