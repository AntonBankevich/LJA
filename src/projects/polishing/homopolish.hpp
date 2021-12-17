#pragma once

#include <sequences/seqio.hpp>
#include "common/logging.hpp"
std::vector<Contig> Polish(logging::Logger &logger, size_t threads,
                                           const std::vector<Contig> &contigs_file,
                                           const std::experimental::filesystem::path &alignments,
                                           const io::Library &reads, size_t dicompress);
