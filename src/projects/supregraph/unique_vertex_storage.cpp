#include "unique_vertex_storage.hpp"

spg::VertexId spg::UniqueVertexStorage::nextOutAfterDelete(spg::Vertex &cur, spg::Vertex &deleted_core) {
    VertexId res;
    for(Edge &edge : cur) {
        if(edge.getFinish() != deleted_core && edge.getFinish() != deleted_core.rc()) {
            if(res.valid())
                return {};
            res = edge.getFinish().getId();
        }
    }
    return res;
}

spg::VertexId spg::UniqueVertexStorage::nextInAfterDelete(spg::Vertex &cur, spg::Vertex &deleted_core) {
    return nextOutAfterDelete(cur.rc(), deleted_core);
}

void spg::UniqueVertexStorage::propagateUniquenessForward(spg::Vertex &uv, spg::Vertex &deleted_core) {
    VERIFY(isUnique(uv));
    VertexId cur = nextOutAfterDelete(uv, deleted_core);
    while(cur.valid() && nextInAfterDelete(*cur, deleted_core).valid() && !isUnique(*cur)) {
        add(*cur);
        cur = nextOutAfterDelete(*cur, deleted_core);
    }
}

void spg::UniqueVertexStorage::propagateUniqueness(spg::Vertex &uv, spg::Vertex &deleted_core) {
    VERIFY(isUnique(uv));
    propagateUniquenessForward(uv, deleted_core);
    propagateUniquenessForward(uv.rc(), deleted_core);
}

void spg::UniqueVertexStorage::add(spg::Vertex &v) {
    unique.emplace(v.getId());
    unique.emplace(v.rc().getId());
}

void spg::UniqueVertexStorage::remove(spg::Vertex &v) {
    unique.erase(v.getId());
    unique.erase(v.rc().getId());
}

bool spg::UniqueVertexStorage::isUnique(spg::Vertex &v) const {
    return unique.find(v.getId()) != unique.end();
}

void spg::UniqueVertexStorage::fireResolveVertex(spg::Vertex &core, const spg::VertexResolutionResult &resolution) {
    for(Vertex &v : resolution.newVertices()) {
        VERIFY(v.inDeg() == 1);
        VERIFY(v.outDeg() == 1);
        if(isUnique(v.front().getFinish()))
            propagateUniqueness(v.front().getFinish(), core);
        if(isUnique(v.rc().front().getFinish()))
            propagateUniqueness(v.rc().front().getFinish(), core);
    }
}
