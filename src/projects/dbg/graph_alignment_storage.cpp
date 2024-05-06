#include <common/logging.hpp>
#include <sequences/sequence.hpp>
#include "graph_alignment_storage.hpp"

using namespace dbg;

void ReadAlignmentStorage::checkCoverage(const SparseDBG &dbg) const {
    VERIFY_MSG(track_cov, "Error: checking coverage using read storage that does not contribute to coverage");
    std::unordered_map<ConstEdgeId, size_t> map;
    for (const Edge &edge: dbg.edges()) {
        map[edge.getId()] = 0;
    }
    for (const ag::AlignedRead<DBGTraits> &read: *this) {
        if (!read.valid())
            continue;
        GraphPath path = read.path.unpack();
        for (auto seg: path) {
            map[seg.contig().getId()] += seg.size();
            map[seg.contig().rc().getId()] += seg.size();
        }
    }
    for (const Edge &edge: dbg.edges()) {
        VERIFY_MSG(edge.intCov() == map[edge.getId()], "Coverage check failed");
    }
}
