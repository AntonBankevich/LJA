#pragma once

#include "compact_path.hpp"
#include "graph_alignment_storage.hpp"
#include <common/logging.hpp>

void SimpleRemoveUncovered(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                     const std::vector<RecordStorage*> &storages, size_t new_extension_size = 0);
void RemoveUncovered(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                            const std::vector<RecordStorage*> &storages, size_t new_extension_size = 0);

class Connection {
public:
    dbg::EdgePosition pos1;
    dbg::EdgePosition pos2;
    Sequence connection;
    Connection(dbg::EdgePosition pos1, dbg::EdgePosition pos2, Sequence connection);
    Connection shrink() const;
};

void AddConnections(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                    const std::vector<RecordStorage*> &storages, const std::vector<Connection> &connections);