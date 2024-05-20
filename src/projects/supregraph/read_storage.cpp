#include "read_storage.hpp"

using namespace spg;

spg::PathStorage::PathStorage(spg::SupreGraph &spg) : spg(&spg){
}

const std::vector<std::pair<spg::ReadDirection, PathIterator>> &
spg::PathIndex::getReadPositions(spg::Edge &edge) const {
    return read_index.at(edge.getId());
}

const std::vector<spg::ReadDirection> &spg::PathIndex::getOutgoingReads(spg::Vertex &v) const {
    return outgoing_index.at(v.getId());
}

spg::ReadRecord &spg::PathStorage::addRead(std::string name, spg::GraphPath path) {
    shrink(path);
    reads.emplace_back(std::move(name), std::move(path));
    return reads.back();
}

void spg::PathIndex::processPassing(Edge &edge, const spg::VertexResolutionResult &resolution) {
    std::vector<std::pair<size_t, PathIterator>> res;
    for(auto &p : read_index[edge.getId()]) {
        ReadDirection dir = p.first;
        auto it = p.second;
        auto nit = it + 1;
        if(nit == dir.end())
            continue;
        if(*it != edge) {
            VERIFY(it->getStart() == edge.getStart());
        }
        Edge &nedge = (*it == edge) ? *nit : *resolution.get(it->getFinish()).second;
        ++nit;
        VERIFY(resolution.contains(edge, nedge));
        Vertex &new_vertex = resolution.get(edge ,nedge);
        Edge &new_edge1 = new_vertex.rc().front().rc();
        Edge &new_edge2 = new_vertex.front();
        dir.rerouteSameSize(it, nit, {new_edge1.getId(), new_edge2.getId()});
        read_index[new_edge1.getId()].emplace_back(dir, it);
        read_index[new_edge2.getId()].emplace_back(dir, it + 1);
    }
}

void spg::PathIndex::fireAddVertex(spg::Vertex &vertex) {
    outgoing_index[vertex.getId()] = {};
    inner_index[vertex.getId()] = {};
}

void spg::PathIndex::fireAddEdge(spg::Edge &edge) {
    read_index[edge.getId()] = {};
}

void spg::PathIndex::fireDeleteVertex(spg::Vertex &vertex) {
    for(ReadDirection dir : getOutgoingReads(vertex)) {
        if(!dir.empty()) {
            VERIFY(dir.getgetStart() == vertex);
            dir.pop_front();
        }
    }
    outgoing_index.erase(vertex.getId());
//    if(vertex != vertex.rc()) {
//        outgoing_index.erase(vertex.rc().getId());
//    }
}

void spg::PathIndex::fireDeleteEdge(spg::Edge &edge) {
    std::cout << "Fire delete " << edge.getId() << std::endl;
    read_index.erase(edge.getId());
}

void spg::PathIndex::fireResolveVertex(spg::Vertex &core, const spg::VertexResolutionResult &resolution) {
    std::cout << "FireResolveVertex " << resolution << std::endl;
    for(Vertex &new_vertex : resolution.newVertices()) {
        fireAddVertex(new_vertex);
        VERIFY(new_vertex.outDeg() == 1);
        VERIFY(new_vertex.inDeg() == 1);
        fireAddEdge(new_vertex.front());
    }
    std::vector<std::pair<size_t, PathIterator>> to_reroute;
    for(Edge &edge : core.incoming()) {
        processPassing(edge, resolution);
    }
    for(ReadDirection dir : getOutgoingReads(core)) {
        if(dir.empty()) {
            inner_index[dir.getgetStart().getId()].emplace_back(dir.getgetStart(), dir.cutLeft(), dir.getgetStart().size() - dir.cutRight());
        } else {
            outgoing_index[dir.begin()->getFinish().getId()].emplace_back(dir);
        }
    }
    for(Vertex &new_vertex : resolution.newVertices()) {
        EdgePair edges = resolution.get(new_vertex);
        for(Segment<Vertex> &seg : inner_index[core.getId()]) {
            size_t left = seg.left + edges.first->rc().truncSize();
            size_t right = left + seg.size();
            inner_index[new_vertex.getId()].emplace_back(new_vertex, left, right);
        }
    }
//    TODO: When real fiddreDeleteVertex exists, need to invalidate all reads with empty paths. Also need a way to preserve their information
}

const std::vector<Segment<Vertex>> &PathIndex::getInnerReads(Vertex &v) const {
    return inner_index.at(v.getId());
}

ComplexIterableStorage<Generator<std::vector<std::pair<ReadDirection, PathIterator>>::iterator, PathIndex::PassingRead>>
PathIndex::getPassing(Vertex &v) &{
    std::function<PassingRead(std::pair<ReadDirection, PathIterator> &)> generate = [](std::pair<ReadDirection, PathIterator> &p) -> PassingRead {
        return {p.first, p.second};
    };
    std::function<bool(std::pair<ReadDirection, PathIterator> &)> use = [](std::pair<ReadDirection, PathIterator> &p) {
        return p.second + 1 != p.first.end();
    };
    ComplexIterableStorage<Generator<std::vector<std::pair<ReadDirection, PathIterator>>::iterator, PassingRead>> res;
    for(Edge &edge : v.incoming()) {
        res += {{read_index.at(edge.getId()).begin(), read_index.at(edge.getId()).end(), generate, use},
                {read_index.at(edge.getId()).end(), read_index.at(edge.getId()).end(), generate, use}};
    }
    return std::move(res);
}

PathIndex::PathIndex(SupreGraph &spg, PathStorage &storage) : ResolutionListener(spg){
    for(Edge &edge : spg.edges()) {
        read_index[edge.getId()] = {};
    }
    for(Vertex &vertex : spg.vertices()) {
        outgoing_index[vertex.getId()] = {};
    }
    for(ReadRecord &rr : storage) {
        if (!rr.path.empty()) {
            outgoing_index[rr.path.start().getId()].emplace_back(rr.forward());
            outgoing_index[rr.path.finish().rc().getId()].emplace_back(rr.backward());
            for (auto it = rr.forward().begin(); it != rr.forward().end(); ++it) {
                Edge &edge = *it;
                read_index[it->getId()].emplace_back(rr.forward(), it);
            }
            for (auto it = rr.backward().begin(); it != rr.backward().end(); ++it) {
                Edge &edge = *it;
                read_index[it->getId()].emplace_back(rr.backward(), it);
            }
        } else if (rr.path.valid()) {
            inner_index[rr.path.start().getId()].emplace_back(rr.path.start(), rr.path.cutLeft(),
                                                                 rr.path.start().size() - rr.path.cutRight());
            inner_index[rr.path.start().rc().getId()].emplace_back(rr.path.start().rc(),
                                                                      rr.path.cutRight(), rr.path.start().size() - rr.path.cutLeft());
        }
    }
}

bool PathIndex::checkReadIndexConsistency() const {
    for(auto &it1 : read_index) {
        for(auto p : it1.second) {
            bool found = false;
            size_t cnt = 0;
            for(auto it = p.first.begin(); it != p.first.end(); ++it) {
                if (it == p.second)
                    found = true;
                cnt++;
            }
            if(!found) {
                std::cout << "Found a problem in read index " << it1.first << " " << p.first << " " << p.second.str() << std::endl;
                return false;
            }
        }
    }
    return true;
}

ReadDirection ReadRecord::forward() {return {*this, false};}

ReadDirection ReadRecord::backward() {return {*this, true};}
