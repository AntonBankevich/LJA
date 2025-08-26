#include "sparse_dbg.hpp"

void dbg::HashListener::fireAddVertex(ag::Vertex &v) {
    if(v.getHash() == ag::Vertex::default_hash && !v.getSeq().empty()) {
        getFire<SparseDBG>().setHash(v, hashing::MovingKWH(hasher, v.getSeq(), 0).hash());
    }
}

dbg::HashListener::HashListener(dbg::SparseDBG &dbg, const hashing::RollingHash &hasher) :
        ag::ResolutionListener(dbg, "HashListener"), hasher(hasher) {}
