#pragma once

#include "dbg_read_alignment_storage.hpp"
#include <common/logging.hpp>

namespace dbg {
    void SimpleRemoveUncovered(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg);
    void SplitUncovered(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                        const std::vector<ag::AlignedReadStorage *> &storages);
    void RemoveUncovered(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                         const std::vector<ag::AlignedReadStorage *> &storages);

    class Connection {
    public:
        VertexId tip1;
        VertexId tip2;
        AlignmentForm al;

        Connection(Edge &tip1, Edge &tip2, AlignmentForm al) : tip1(tip1.getFinish().getId()), tip2(tip2.getFinish().getId()), al(std::move(al)) {}
    };
}