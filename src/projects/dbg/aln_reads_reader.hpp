//
// Created by Anton Zamyatin on 5/20/24.
//
#pragma once
#include "sparse_dbg.hpp"
#include "sequences/seqio.hpp"
#include "graph_alignment_storage.hpp"

namespace dbg {
    class DBGAlignedReadsReader : public io::IContigFromFileReader {
    private:
        SparseDBG * sparse_dbg = nullptr;
        IdIndex<Vertex> id_index;
    public:
        explicit DBGAlignedReadsReader(SparseDBG & _sparse_dbg, const std::experimental::filesystem::path& _file_name);
        void inner_read() override;
    };
};