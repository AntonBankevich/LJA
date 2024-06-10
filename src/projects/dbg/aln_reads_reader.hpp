//
// Created by Anton Zamyatin on 5/20/24.
//
#pragma once
#include "sparse_dbg.hpp"
#include "sequences/seqio.hpp"
#include "graph_alignment_storage.hpp"
#include <map>

namespace dbg {

    class SeqReader : public io::ISeqReader {
    private:
        size_t threads;
        logging::Logger &logger;
        std::map<std::experimental::filesystem::path, SparseDBG> loaded_graphs;

    public:
        void initReader(const std::experimental::filesystem::path &file_name) override;
        explicit SeqReader(io::Library _lib, logging::Logger &_logger, size_t threads,
                           size_t _min_read_size = size_t(-1) / 2, size_t _overlap = size_t(-1) / 8);
        explicit SeqReader(const std::experimental::filesystem::path & file_name,
                           logging::Logger &_logger, size_t threads,
                           size_t _min_read_size = size_t(-1) / 2, size_t _overlap = size_t(-1) / 8);
    };

    class DBGAlignedReadsReader final : public io::IContigFromFileReader {
    private:
        SparseDBG * sparse_dbg = nullptr;
        IdIndex<Vertex> id_index;
        size_t number_of_paths;
    public:
        explicit DBGAlignedReadsReader(SparseDBG &_sparse_dbg, const std::experimental::filesystem::path& _file_name);
        void inner_read() override;
    };
};