#pragma once

#include "common/string_utils.hpp"
#include "contigs.hpp"
#include <experimental/filesystem>
#include <vector>

namespace io{
//    Files of reads, read by io::ReadLoader.
    typedef std::vector<std::experimental::filesystem::path> Library;
}

inline io::Library operator+(const io::Library &lib1, const io::Library &lib2) {
    io::Library res = lib1;
    res.insert(res.end(), lib2.begin(), lib2.end());
    return res;
}
