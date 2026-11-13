#include "sparse_dbg.hpp"

void dbg::HashListener::fireAddVertex(ag::Vertex &v) {
    if(v.getHash() == ag::Vertex::default_hash && !v.getSeq().empty()) {
        SparseDBG &dbg = getFire<SparseDBG>();
        VERIFY(v.getSeq().size() == dbg.getK());
        dbg.setHash(v, hashing::MovingKWH(dbg.hasher(), v.getSeq(), 0).hash());
    }
}

bool dbg::HashListener::checkConsistency(std::ostream &os) {
    SparseDBG &dbg = getFire<SparseDBG>();
    bool ok = true;
    std::unordered_set<hashing::htype> seen;
    for(ag::Vertex &v : dbg.vertices()) {
        if(v.getHash() != v.rc().getHash()) {
            os << "HashListener: vertex " << v.getId() << " and its rc have different hashes" << std::endl;
            ok = false;
        }
        if(v.getSeq().empty()) {
            if(v.getHash() == ag::Vertex::default_hash) {
                os << "HashListener: vertex " << v.getId() << " has neither sequence nor hash" << std::endl;
                ok = false;
            }
        } else if(v.getSeq().size() != dbg.getK()) {
            os << "HashListener: vertex " << v.getId() << " has size " << v.getSeq().size() << " != k" << std::endl;
            ok = false;
        } else if(hashing::MovingKWH(dbg.hasher(), v.getSeq(), 0).hash() != v.getHash()) {
            os << "HashListener: vertex " << v.getId() << " hash does not match its sequence" << std::endl;
            ok = false;
        }
        if(v.isCanonical() && !seen.insert(v.getHash()).second) {
            os << "HashListener: hash of vertex " << v.getId() << " is not unique" << std::endl;
            ok = false;
        }
    }
    return ok;
}

dbg::HashListener::HashListener(dbg::SparseDBG &dbg) : DBGListener(dbg, "HashListener") {}
