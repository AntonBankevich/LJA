//
// Created by Anton Zamyatin on 5/20/24.
//
#pragma once
#include "sparse_dbg.hpp"
#include "sequences/read_loader.hpp"
#include "dbg_read_alignment_storage.hpp"
#include <map>

namespace dbg {
//    io::ReadLoader that also reads .paths files: reads aligned to a graph, read as the sequences of their paths and
//    not compressed, as the graph already is. A .paths file names a graph (GFA) on its first line and the read
//    alignments to it (.aln) on its second, and gives k on its third. Each graph is loaded once.
    class ReadLoader : public io::ReadLoader {
    private:
        std::map<std::experimental::filesystem::path, SparseDBG> graphs;
    protected:
        std::unique_ptr<io::ReadSource> open(const std::experimental::filesystem::path &path) override;
    public:
        using io::ReadLoader::ReadLoader;
    };
}
