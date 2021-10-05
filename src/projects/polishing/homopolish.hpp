#pragma once

#include <sequences/seqio.hpp>
#include "common/logging.hpp"
std::experimental::filesystem::path Polish(logging::Logger &logger, size_t threads,
                                           const std::experimental::filesystem::path &output_file,
                                           const std::experimental::filesystem::path &contigs_file,
                                           const std::experimental::filesystem::path &alignments,
                                           const io::Library &reads, size_t dicompress);
