#include "vertex_resolution.hpp"

void spg::VertexResolutionPlan::add(const spg::InOutEdgePair &edgePair) {
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

IterableStorage<std::vector<spg::InOutEdgePair>::const_iterator> spg::VertexResolutionPlan::connections() const {
    sort();
    return {edge_pairs.begin(), edge_pairs.end()};
}

IterableStorage<SkippingIterator<std::vector<spg::InOutEdgePair>::const_iterator>>
spg::VertexResolutionPlan::connectionsUnique() const {
    sort();
    std::function<bool(const InOutEdgePair &)> use = [](const InOutEdgePair &ep)->bool {
        return ep.middle() != ep.middle().rc() || ep.incoming().getId() < ep.outgoing().rc().getId() ||
                    (ep.incoming() == ep.outgoing().rc() && ep.outgoing().getId() < ep.incoming().rc().getId());
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
        if(edge == it.incoming())
            return true;
    return false;
}

bool spg::VertexResolutionPlan::outConnected(spg::Edge &edge) const {
    for(const auto &it : edge_pairs)
        if(edge == it.outgoing())
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
