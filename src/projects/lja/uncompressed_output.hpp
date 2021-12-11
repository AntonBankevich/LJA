#pragma once

#include <common/logging.hpp>
#include "multi_graph.hpp"

void printUncompressedResults(logging::Logger &logger, size_t threads, multigraph::MultiGraph &graph,
                              const std::vector<Contig> &uncompressed, const std::experimental::filesystem::path &out_dir);