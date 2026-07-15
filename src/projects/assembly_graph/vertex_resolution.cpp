#include "vertex_resolution.hpp"

using namespace ag;

std::ostream &ag::operator<<(std::ostream &stream, const InOutEdgePair &pair) {
    return stream << "(" << pair.first << "|" << pair.second << ")";
}

InOutEdgePair::InOutEdgePair(Edge &first, Edge &second) : first(first.getId()), second(second.getId()) {
    VERIFY(this->first->getFinish() == this->second->getStart());
}

void VertexResolutionResult::innerAdd(Vertex &new_vertex, const InOutEdgePair &edgePair) {
    VERIFY(new_vertices.find(new_vertex.getId()) == new_vertices.end());
    new_vertices.emplace(new_vertex.getId(), edgePair);
    edge_mapping[edgePair.incoming().getId()][edgePair.outgoing().getId()] = new_vertex.getId();
}

VertexResolutionResult VertexResolutionResult::RC() const {
    VertexResolutionResult res(core->rc());
    for(auto it : new_vertices) {
        res.innerAdd(it.first->rc(), it.second.RC());
    }
    return std::move(res);
}

bool VertexResolutionResult::contains(Edge &edge1, Edge &edge2) const {
    VERIFY(edge1.getFinish() == *core);
    VERIFY(edge2.getStart() == *core);
    return edge_mapping.find(edge1.getId()) != edge_mapping.end() &&
           edge_mapping.at(edge1.getId()).find(edge2.getId()) != edge_mapping.at(edge1.getId()).end();
}

Vertex &VertexResolutionResult::get(Edge &edge1, Edge &edge2) const {
    return *edge_mapping.at(edge1.getId()).at(edge2.getId());
}

const InOutEdgePair &VertexResolutionResult::get(Vertex &new_vertex) const {
    return new_vertices.at(new_vertex.getId());
}

void VertexResolutionResult::add(Vertex &new_vertex, const InOutEdgePair &edgePair) {
    innerAdd(new_vertex, edgePair);
    if(*core == core->rc() && new_vertex != new_vertex.rc()) {
        innerAdd(new_vertex.rc(), edgePair.RC());
    }
}

void VertexResolutionResult::add(Vertex &new_vertex, Edge &edge1, Edge &edge2) {add(new_vertex, {edge1, edge2});}

IterableStorage<TransformingIterator<typename std::unordered_map<VertexId, InOutEdgePair>::const_iterator, Vertex>>
VertexResolutionResult::newVertices() const {
    std::function<Vertex &(const std::pair<VertexId, InOutEdgePair> &)> transform = [](const std::pair<VertexId, InOutEdgePair> &val) ->Vertex& {
        return *val.first;
    };
    return {{new_vertices.begin(), new_vertices.end(), transform},
            {new_vertices.end(),   new_vertices.end(), transform}};
}

std::ostream &ag::operator<<(std::ostream &stream, const VertexResolutionResult &vr) {
    stream << "VRResult." << vr.getCore().getId() << ":";
    for(Vertex & it : vr.newVertices()) {
        stream << it.getId() << vr.get(it);
    }
    return stream;
}

void VertexResolutionPlan::add(const InOutEdgePair &edgePair) {
    VERIFY(edgePair.incoming().getFinish().getId() == v && edgePair.outgoing().getStart().getId() == v);
    for(const InOutEdgePair &ep : edge_pairs)
        if(edgePair == ep)
            return;
    edge_pairs.emplace_back(edgePair);
    InOutEdgePair rc = edgePair.RC();
    if(*v == v->rc() && rc != edgePair)
        edge_pairs.emplace_back(rc);
    sorted = false;
}

IterableStorage<std::vector<InOutEdgePair>::const_iterator> VertexResolutionPlan::connections() const {
    sort();
    return {edge_pairs.begin(), edge_pairs.end()};
}

IterableStorage<SkippingIterator<std::vector<InOutEdgePair>::const_iterator>>
VertexResolutionPlan::connectionsUnique() const {
    sort();
    std::function<bool(const InOutEdgePair &)> use = [](const InOutEdgePair &ep)->bool {
        return ep.middle() != ep.middle().rc() || ep.incoming().getId() < ep.outgoing().rc().getId() ||
               (ep.incoming() == ep.outgoing().rc() && ep.outgoing().getId() < ep.incoming().rc().getId());
    };
    return {{edge_pairs.begin(), edge_pairs.end(), use}, {edge_pairs.end(), edge_pairs.end(), use}};
}

void VertexResolutionPlan::sort() const {
    if(!sorted)
        std::sort(edge_pairs.begin(), edge_pairs.end());
    sorted = true;
}

bool VertexResolutionPlan::incConnected(Edge &edge) const {
    for(const auto &it : edge_pairs)
        if(edge == it.incoming())
            return true;
    return false;
}

bool VertexResolutionPlan::outConnected(Edge &edge) const {
    for(const auto &it : edge_pairs)
        if(edge == it.outgoing())
            return true;
    return false;
}

bool VertexResolutionPlan::incConnected() const {
    for(Edge &edge : v->incoming())
        if(!incConnected(edge))
            return false;
    return true;

}

bool VertexResolutionPlan::outConnected() const {
    for(Edge &edge : *v)
        if(!outConnected(edge))
            return false;
    return true;

}

bool VertexResolutionPlan::allConnected() const {
    for(Edge &edge : v->incoming())
        if(!incConnected(edge))
            return false;
    for(Edge &edge : *v)
        if(!outConnected(edge))
            return false;
    return true;
}

VertexResolutionPlan VertexResolutionPlan::RC() const {
    VertexResolutionPlan res(v->rc());
    for(const InOutEdgePair &ep : edge_pairs) {
        res.add(ep.RC());
    }
    return std::move(res);
}

std::ostream &ag::operator<<(std::ostream &stream, const VertexResolutionPlan &vr) {
    stream << "VRResult." << vr.getCore().getId() << ":";
    for(const InOutEdgePair & it : vr.connections()) {
        stream << "(" << it.incoming().getId() << "|" << it.outgoing().getId() << ")";
    }
    return stream;
}