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
    unique.insert(v.getId(), true);
    unique.insert(v.rc().getId(), true);
}

void spg::UniqueVertexStorage::remove(const Vertex &v) {
    unique.erase(v.getId());
    unique.erase(v.rc().getId());
}

bool spg::UniqueVertexStorage::isUnique(const Vertex &v) const {
    return unique.contains(v.getId());
}

void spg::UniqueVertexStorage::fireMergePath(const ag::RAGraphPath &path, Vertex &vertex) {
    for(Vertex &v : path.vertices()) {
        if(isUnique(v)) {
            add(vertex);
            return;
        }
    }
}

void spg::UniqueVertexStorage::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
//    v spells e, so it contains both of e's ends.
    if (isUnique(e.getStart()) || isUnique(e.getFinish()))
        add(v);
}

void spg::UniqueVertexStorage::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
//    Every new vertex is a + core + b.
    if (isUnique(core))
        for (Vertex &newv : resolution.newVertices())
            add(newv);
}

void spg::UniqueVertexStorage::fireMergeLoop(const ag::GraphPath &path, Vertex &vertex) {
    fireMergePath(path.asRAPath(), vertex);
}
spg::UniqueVertexStorage::UniqueVertexStorage(ag::AssemblyGraph &spg, const std::function<bool(Vertex &)> &is_unique, size_t unique_threshold)
            : SPGListener(spg, "UniqueVertexStorage"), unique_threshold(unique_threshold){
    for (Vertex &vertex : spg.vertices()) {fireAddVertex(vertex);}
    for(Vertex &vertex : spg.vertices()) {
        if(is_unique(vertex))
            add(vertex);
    }
}

spg::UniqueVertexStorage::UniqueVertexStorage(ag::AssemblyGraph &spg, size_t unique_threshold): SPGListener(spg, "UniqueVertexStorage"), unique_threshold(unique_threshold) {
    for (Vertex &vertex : spg.vertices()) {fireAddVertex(vertex);}
}

