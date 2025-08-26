#pragma once

#include <common/logging.hpp>
#include "dbg/multi_graph.hpp"

void printUncompressedResults(logging::Logger &logger, size_t threads, multigraph::MultiGraph &graph,
                              const std::unordered_map<ag::VertexId, Segment<ag::Vertex>> &segs,
                              const std::unordered_map<ag::VertexId , Sequence> &uncompressed_results,
                              const std::experimental::filesystem::path &out_dir, bool debug);