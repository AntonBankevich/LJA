#pragma once

#include "assembly_graph/compact_path.hpp"
#include "graph_alignment_storage.hpp"
#include <common/logging.hpp>

namespace dbg {
    void SimpleRemoveUncovered(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                               const std::vector<dbg::ReadAlignmentStorage *> &storages, size_t new_extension_size = 0);

    void RemoveUncovered(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                         const std::vector<ReadAlignmentStorage *> &storages, size_t new_extension_size = 0);

    class Connection {
    public:
        dbg::EdgePosition pos1;
        dbg::EdgePosition pos2;
        Sequence connection;

        Connection(dbg::EdgePosition pos1, dbg::EdgePosition pos2, Sequence connection);

        Connection shrink() const;
    };

    dbg::SparseDBG AddConnections(logging::Logger &logger, size_t threads, const dbg::SparseDBG &dbg,
                                  const std::vector<ReadAlignmentStorage*> &storages,
                                  const std::vector<Connection> &connections);
}