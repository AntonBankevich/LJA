#include "decision_rules.hpp"

size_t spg::ChainRule::getDiveSize(spg::Edge &edge) {
    size_t res = 0;
    for(auto it : storage->getReadPositions(edge)) {
        ReadDirection dir = it.first;
        PathIterator pos = it.second;
        if(pos == dir.begin())
            res = std::max(res, edge.getStart().size() - dir.cutLeft());
    }
    return res;
}

spg::VertexResolutionPlan spg::ChainRule::judge(spg::Vertex &v) {
    VertexResolutionPlan res(v);
    bool has_covering = false;
    for(auto it : storage->getPassing(v)) {
        res.add(*it.edges.first, *it.edges.second);
        has_covering = true;
    }
    std::vector<Segment<Vertex>> segs = storage->getInnerReads(v);
    std::sort(segs.begin(), segs.end());
    bool has_unpassable = false;
    for(Edge &inc : v.incoming()) {
        size_t max = getDiveSize(inc.rc());
        for(Segment<Vertex> &p : segs) {
            if(max >= p.left + k) {
                max = std::max(max, p.right);
            }
        }
        for(Edge &out : v) {
            size_t min = v.size() - getDiveSize(out);
            if(min + k <= max) {
                res.add(inc, out);
            } else {
                has_unpassable = true;
            }
        }
    }
    if(has_covering || has_unpassable)
        return std::move(res);
    else
        return {v};
}

void spg::AndreyRule::loopHeuristic(spg::VertexResolutionPlan &res) const {
    Vertex &core = res.getCore();
    size_t loop_cnt = 0;
    EdgeId loop_start;
    EdgeId loop_end;
    for(Edge &e : core) {
        ag::GraphPath<SPGTraits> path = ag::GraphPath<SPGTraits>::WalkForward(e);
        if(path.finish() == core) {
            loop_start = e.getId();
            loop_end = path.backEdge().getId();
            break;
        }
    }
    if(loop_start.valid() && loop_end.valid() && unique_storage->isUnique(loop_start->getFinish())) {
        if(core.inDeg() == 2)
            for(Edge &e : core.incoming())
                if(e != *loop_end)
                    res.add(e, *loop_start);
        if(core.outDeg() == 2)
            for(Edge &e : core)
                if(e != *loop_start)
                    res.add(*loop_end, e);
    }
}

spg::EdgeId spg::AndreyRule::getUniqueDisconnectedInc(const spg::VertexResolutionPlan &plan) {
    Vertex &v = plan.getCore();
    EdgeId res;
    for(Edge &e: v.incoming()) {
        if(plan.incConnected(e)) {
            if (!unique_storage->isUnique(e.getFinish()))
                return {};
        } else {
            if(res.valid())
                return {};
            else
                res = e.getId();
        }
    }
    return res;
}

void spg::AndreyRule::uniqueHeuristic(spg::VertexResolutionPlan &res) {
    EdgeId disconnected_inc = getUniqueDisconnectedInc(res);
    EdgeId disconnected_out = getUniqueDisconnectedInc(res.RC());
    if(disconnected_inc.valid() && disconnected_out.valid()) {
        res.add(*disconnected_inc, disconnected_out->rc());
    }
}

void spg::AndreyRule::noChoiceHeuristic(spg::VertexResolutionPlan &res) {
    Vertex &core = res.getCore();
    if(core.inDeg() == 1 && !core.isInfLeft()) {
        for(Edge &edge : core) {
            res.add(*core.incoming().begin(), edge);
        }
    }
    if(core.outDeg() == 1 && !core.isInfRight()) {
        for(Edge &edge : core.incoming()) {
            res.add(edge, core.front());
        }
    }
}

spg::VertexResolutionPlan spg::AndreyRule::judge(spg::Vertex &v) {
    VertexResolutionPlan res(v);
    const std::vector<std::pair<spg::ReadDirection, PathIterator>> &tmp = storage->getReadPositions(*v.incoming().begin());
    for(auto it : storage->getPassing(v)) {
        res.add(*it.edges.first, *it.edges.second);
    }
    std::cout << "Passing: " << res << std::endl;
    loopHeuristic(res);
    std::cout << "Loop: " << res << std::endl;
    uniqueHeuristic(res);
    std::cout << "Unique: " << res << std::endl;
    noChoiceHeuristic(res);
    std::cout << "NoChoice: " << res << std::endl;
    if(res.allConnected())
        return std::move(res);
    return {v};
}
