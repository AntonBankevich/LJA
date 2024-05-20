#include "vertex_resolution.hpp"

void spg::VertexResolutionPlan::add(const spg::EdgePair &edgePair) {
    VERIFY(edgePair.first->getFinish().getId() == v && edgePair.second->getStart().getId() == v);
    for(const EdgePair &ep : edge_pairs)
        if(edgePair == ep)
            return;
    edge_pairs.emplace_back(edgePair);
    EdgePair rc = edgePair.RC();
    if(*v == v->rc() && rc != edgePair)
        edge_pairs.emplace_back(rc);
    sorted = false;
}

IterableStorage<std::vector<spg::EdgePair>::const_iterator> spg::VertexResolutionPlan::connections() const {
    sort();
    return {edge_pairs.begin(), edge_pairs.end()};
}

IterableStorage<SkippingIterator<std::vector<spg::EdgePair>::const_iterator>>
spg::VertexResolutionPlan::connectionsUnique() const {
    sort();
    std::function<bool(const EdgePair &)> use = [](const EdgePair &ep)->bool {
        return ep.middle() != ep.middle().rc() || ep.first < ep.second->rc().getId() || ep.second < ep.first->rc().getId();
    };
    return {{edge_pairs.begin(), edge_pairs.end(), use}, {edge_pairs.end(), edge_pairs.end(), use}};
}

void spg::VertexResolutionPlan::sort() const {
    if(!sorted)
        std::sort(edge_pairs.begin(), edge_pairs.end());
    sorted = true;
}

bool spg::VertexResolutionPlan::incConnected(spg::Edge &edge) const {
    for(const auto &it : edge_pairs)
        if(edge == *it.first)
            return true;
    return false;
}

bool spg::VertexResolutionPlan::outConnected(spg::Edge &edge) const {
    for(const auto &it : edge_pairs)
        if(edge == *it.second)
            return true;
    return false;
}

bool spg::VertexResolutionPlan::allConnected() const {
    for(Edge &edge : v->incoming())
        if(!incConnected(edge))
            return false;
    for(Edge &edge : *v)
        if(!outConnected(edge))
            return false;
    return true;
}
