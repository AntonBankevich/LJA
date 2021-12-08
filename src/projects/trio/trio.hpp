#pragma once


#include <array>
#include <vector>
#include <unordered_map>
#include <sequences/seqio.hpp>
#include "common/logging.hpp"
#include <common/string_utils.hpp>
#include "lja/multi_graph.hpp"
#include "lja/subdataset_processing.hpp"



std::experimental::filesystem::path simplifyHaplo(logging::Logger &logger, size_t threads,
                                           const std::experimental::filesystem::path &output_file,
                                           const std::experimental::filesystem::path &diplo_graph,
                                           const std::experimental::filesystem::path &haployak,
                                           const char haplotype);

