//
// Created by anton on 09.07.2020.
//

#pragma once

#include <string>
#include <sys/stat.h>
#include <experimental/filesystem>

//TODO: throw exception if this is file
inline void ensure_dir_existance(const std::experimental::filesystem::path & path) {
    struct stat statbuf{};
    if (not std::experimental::filesystem::is_directory(path)) {
        std::experimental::filesystem::create_directories(path);
    }
}

//TODO: throw exception if this is file
inline void recreate_dir(const std::experimental::filesystem::path & path) {
    struct stat statbuf{};
    if (std::experimental::filesystem::is_directory(path)) {
        std::experimental::filesystem::remove_all(path);
    }
    std::experimental::filesystem::create_directories(path);
}