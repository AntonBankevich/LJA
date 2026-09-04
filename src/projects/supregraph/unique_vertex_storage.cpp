#include "unique_vertex_storage.hpp"
#include "assembly_graph/assembly_graph.hpp"

ag::VertexId spg::UniqueVertexStorage::nextOutAfterDelete(Vertex &cur, Vertex &deleted_core) {
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

ag::VertexId spg::UniqueVertexStorage::nextInAfterDelete(Vertex &cur, Vertex &deleted_core) {
    return nextOutAfterDelete(cur.rc(), deleted_core);
}

void spg::UniqueVertexStorage::propagateUniquenessForward(Vertex &uv, Vertex &deleted_core) {
    VERIFY(isUnique(uv));
    VertexId cur = nextOutAfterDelete(uv, deleted_core);
    while(cur.valid() && nextInAfterDelete(*cur, deleted_core).valid() && !isUnique(*cur)) {
        add(*cur);
        cur = nextOutAfterDelete(*cur, deleted_core);
    }
}

void spg::UniqueVertexStorage::propagateUniqueness(Vertex &uv, Vertex &deleted_core) {
    VERIFY(isUnique(uv));
    propagateUniquenessForward(uv, deleted_core);
    propagateUniquenessForward(uv.rc(), deleted_core);
}

void spg::UniqueVertexStorage::add(const Vertex &v) {
    unique.emplace(v.getId());
    unique.emplace(v.rc().getId());
}

void spg::UniqueVertexStorage::remove(const Vertex &v) {
    unique.erase(v.getId());
    unique.erase(v.rc().getId());
}

bool spg::UniqueVertexStorage::isUnique(const Vertex &v) const {
    return unique.find(v.getId()) != unique.end();
}

void spg::UniqueVertexStorage::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
    for(Vertex &v : resolution.newVertices()) {
        VERIFY(v.inDeg() == 1);
        VERIFY(v.outDeg() == 1);
        if(isUnique(v.front().getFinish()))
            propagateUniqueness(v.front().getFinish(), core);
        if(isUnique(v.rc().front().getFinish()))
            propagateUniqueness(v.rc().front().getFinish(), core);
    }
}

void spg::UniqueVertexStorage::fireMergePath(const ag::RAGraphPath &path, Vertex &vertex) {
    for(Edge &e : path.edges()) {
        if(isUnique(e.getStart())) {
            add(vertex);
            return;
        }
    }
    if(isUnique(path.backEdge().getFinish())) {
        add(vertex);
    }
}
void spg::UniqueVertexStorage::fireMergeLoop(const ag::GraphPath &path, Vertex &vertex) {
    fireMergePath(path.asRAPath(), vertex);
}
spg::UniqueVertexStorage::UniqueVertexStorage(ag::AssemblyGraph &spg, const std::function<bool(Vertex &)> &is_unique, size_t unique_threshold)
            : ag::ResolutionListener(spg, "UniqueVertexStorage"), unique_threshold(unique_threshold){
    for (Vertex &vertex : spg.vertices()) {fireAddVertex(vertex);}
    for(Vertex &vertex : spg.vertices()) {
        if(is_unique(vertex))
            add(vertex);
    }
}

spg::UniqueVertexStorage::UniqueVertexStorage(ag::AssemblyGraph &spg, size_t unique_threshold): ag::ResolutionListener(spg, "UniqueVertexStorage"), unique_threshold(unique_threshold) {
    for (Vertex &vertex : spg.vertices()) {fireAddVertex(vertex);}
}

