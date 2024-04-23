#include "supregraph_base.hpp"

bool spg::SPGVertex::isCore() {
    for(Edge &edge: *this) {
        if(edge.isSuffix())
            return false;
    }
    for(Edge &edge: rc()) {
        if(edge.isSuffix())
            return false;
    }
    return true;
}

spg::SPGEdge &spg::SPGVertex::addSPEdgeLockFree(spg::SPGVertex &end, ag::BaseEdge<spg::SPGTraits>::id_type eid,
                                                ag::BaseEdge<spg::SPGTraits>::id_type rcid) {
    if (end == *this) {
        return addEdgeLockFree(end, getSeq() + getSeq(), EdgeData(), eid, rcid);
    } else {
        VERIFY(size() != end.size());
        if (size() < end.size()) {
            VERIFY(end.getSeq().startsWith(getSeq()));
            return addEdgeLockFree(end, end.getSeq(), EdgeData(), eid, rcid);
        } else {
            VERIFY(size() > end.size())
            VERIFY(getSeq().endsWith(end.getSeq()));
            return addEdgeLockFree(end, getSeq(), EdgeData(), eid, rcid);
        }
    }
}

spg::SPGEdge &spg::SPGVertex::addSPEdge(spg::SPGVertex &end, ag::BaseEdge<spg::SPGTraits>::id_type eid,
                                        ag::BaseEdge<spg::SPGTraits>::id_type rcid) {
    ag::Locker <BaseVertex<SPGTraits>> locker({this, &end.rc()});
    return addSPEdgeLockFree(end, eid, rcid);
}
