#include "read_storage.hpp"
#include "listeners.hpp"

using namespace spg;

spg::PathStorage::PathStorage(spg::SupreGraph &spg) : spg(&spg){
}

const std::vector<std::pair<spg::ReadDirection, PathIterator>> &
spg::PathIndex::getReadPositions(spg::Edge &edge) const {
    prepare();
    return read_index.at(edge.getId());
}

spg::ReadRecord &spg::PathStorage::addRead(std::string name, spg::GraphPath path) {
    shrink(path);
    reads.emplace_back(std::move(name), std::move(path));
    return reads.back();
}

ReadRecord &PathStorage::Oppa() {return reads[4078];}

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
        size_t cut_left = it == dir.begin() ? dir.cutLeft() : 0;
        size_t cut_right = nit == dir.end() ? dir.cutRight() : 0;
        dir.rerouteSameSize(it, nit, GraphPath(new_edge1.getStart(), std::vector<EdgeId>({new_edge1.getId(), new_edge2.getId()}), cut_left, cut_right));
        read_index[new_edge1.getId()].emplace_back(dir, it);
        read_index[new_edge2.getId()].emplace_back(dir, it + 1);
    }
}

void spg::PathIndex::fireAddVertex(spg::Vertex &vertex) {
    std::cout << "Fire Add Vertex " << vertex.getId();
    if(vertex.size() < 10)
        std::cout << " " << vertex.getSeq();
    std::cout << std::endl;
    inner_index[vertex.getId()] = {};
}

void spg::PathIndex::fireAddEdge(spg::Edge &edge) {
    std::cout << "Fire Add Edge " << edge.getId() << " " << edge.truncSize() << std::endl;
    read_index[edge.getId()] = {};
}

void spg::PathIndex::fireDeleteVertex(spg::Vertex &vertex) {
//    if(vertex != vertex.rc()) {
//        outgoing_index.erase(vertex.rc().getId());
//    }
}

void spg::PathIndex::fireDeleteEdge(spg::Edge &edge) {
    std::cout << "Fire delete " << edge.getId() << std::endl;
//    if(edge.isPrefix())
//        for(auto it : read_index[edge.getId()]) {
//            ReadDirection dir = it.first;
//            if(it.second == dir.begin()) {
//                dir.pop_front();
//                if(dir.empty()) {
//                    inner_index[dir.getgetStart().getId()].emplace_back(dir.getgetStart(), dir.cutLeft(), dir.getgetStart().size() - dir.cutRight());
//                    dir.invalidate();
//                }
//            }
//        }
    read_index.erase(edge.getId());
}

void spg::PathIndex::fireResolveVertex(spg::Vertex &core, const spg::VertexResolutionResult &resolution) {
    std::cout << "FireResolveVertex " << resolution << std::endl;
    prepare();
    std::vector<std::pair<size_t, PathIterator>> to_reroute;
    for(Edge &edge : core) {
        for(auto &p : read_index[edge.getId()]) {
            ReadDirection dir = p.first;
            auto it = p.second;
            if(*it != edge) {
                VERIFY(it->getFinish() == edge.getFinish());
                continue;
            }
            if(it == dir.begin()) {
                dir.pop_front();
                if(dir.empty()) {
                    inner_index[dir.getgetStart().getId()].emplace_back(dir.getgetStart(), dir.cutLeft(), dir.getgetStart().size() - dir.cutRight());
                    dir.invalidate();
                }
            } else {
                Edge &pedge = *(it - 1);
                VERIFY(resolution.contains(pedge, edge));
                Vertex &new_vertex = resolution.get(pedge, edge);
                Edge &new_edge1 = *new_vertex.incoming().begin();
                Edge &new_edge2 = new_vertex.front();
                size_t cut_left = it - 1 == dir.begin() ? dir.cutLeft() : 0;
                size_t cut_right = it + 1 == dir.end() ? dir.cutRight() : 0;
                std::cout << "Reroute " << dir.getRead().name << std::endl;
                std::cout << dir << std::endl;
                dir.rerouteSameSize(it - 1, it + 1, GraphPath(new_edge1.getStart(),
                                                       std::vector<EdgeId>({new_edge1.getId(), new_edge2.getId()}),
                                                       cut_left, cut_right));
                std::cout << dir << std::endl;
                read_index[new_edge1.getId()].emplace_back(dir, it - 1);
                PathIterator it1 = (it - 1).SameElementRC();
                read_index[new_edge1.rc().getId()].emplace_back(dir.RC(), (it - 1).SameElementRC());
                read_index[new_edge2.getId()].emplace_back(dir, it);
                read_index[new_edge2.rc().getId()].emplace_back(dir.RC(), it.SameElementRC());
            }
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
    for(ReadRecord &rec : *storage) {
        rec.path.lenStr();
    }
//    TODO: When real fiddreDeleteVertex exists, need to invalidate all reads with empty paths. Also need a way to preserve their information
}

const std::vector<Segment<Vertex>> &PathIndex::getInnerReads(Vertex &v) const {
    return inner_index.at(v.getId());
}

ComplexIterableStorage<Generator<std::vector<std::pair<ReadDirection, PathIterator>>::iterator, PathIndex::PassingRead>>
PathIndex::getPassing(Vertex &v) const &{
    prepare();
    std::function<PassingRead(std::pair<ReadDirection, PathIterator> &)> generate = [&v](std::pair<ReadDirection, PathIterator> &p) -> PassingRead {
        std::cout << "Passing: " << v.getId() << " " << p.first.getRead().name << " " << p.second.str() << std::endl;
        std::cout << p.second->getId() << " " << (p.second + 1)->getId() << std::endl;
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

PathIndex::PathIndex(SupreGraph &spg, PathStorage &storage) : ResolutionListener(spg), storage(&storage) {
    for(Edge &edge : spg.edges()) {
        read_index[edge.getId()] = {};
    }
    for(Vertex &vertex : spg.vertices()) {
        inner_index[vertex.getId()] = {};
    }
    prepareIndex();
    reads_ready = true;
}

bool PathIndex::checkReadIndexConsistency() const {
    prepare();
    for(auto &it1 : read_index) {
        for(auto p : it1.second) {
            bool found = false;
            for(auto it = p.first.begin(); it != p.first.end(); ++it) {
                if (it == p.second)
                    found = true;
            }
            if(!found) {
                std::cout << "Found a problem in read index " << it1.first << " " << p.first << " " << p.second.str() << std::endl;
                return false;
            }
        }
    }
    return true;
}

void PathIndex::prepareIndex() const {
    for(ReadRecord &rr : *this->storage) {
        if (!rr.path.empty()) {
            for (auto it = rr.forward().begin(); it != rr.forward().end(); ++it)
                read_index[it->getId()].emplace_back(rr.forward(), it);
            for (auto it = rr.backward().begin(); it != rr.backward().end(); ++it)
                read_index[it->getId()].emplace_back(rr.backward(), it);
        } else if (rr.path.valid()) {
            inner_index[rr.path.start().getId()].emplace_back(rr.path.start(), rr.path.cutLeft(),
                                                              rr.path.start().size() - rr.path.cutRight());
            inner_index[rr.path.start().rc().getId()].emplace_back(rr.path.start().rc(),
                                                                   rr.path.cutRight(), rr.path.start().size() - rr.path.cutLeft());
            rr.path.invalidate();
        }
    }
}

void PathIndex::fireMergePath(const GraphPath &path, Vertex &new_vertex) {
    unprepare();
    std::cout << "FIreMergePath " << path.lenStr() << std::endl;
    embedding.remap(path, new_vertex);
    size_t right = path.start().size();
    size_t left = 0;
    for(Edge &edge : path.edges()) {
        if(edge == path.backEdge())
            break;
        Vertex &vertex = edge.getFinish();
        right += edge.truncSize();
        left = right - vertex.size();
        Segment<Vertex> vertex_embedding(new_vertex, left, right);
        for(const Segment<Vertex> &seg : this->getInnerReads(vertex)) {
            inner_index[new_vertex.getId()].emplace_back(seg.nest(vertex_embedding));
        }
    }
}

void PathIndex::fireMergeLoop(const GraphPath &path, Vertex &new_vertex) {
    VERIFY_MSG(false, "Not supported yet");
}

ReadDirection ReadRecord::forward() {return {*this, false};}

ReadDirection ReadRecord::backward() {return {*this, true};}

Segment<Edge> Embedding::embed(EdgeId eid, size_t cut_left, size_t cut_right) const {
    IdSegment<EdgeId> res(eid, cut_left, cut_right);
    while(edge_embedding.find(res.id) != edge_embedding.end()) {
        res = res.nest(edge_embedding.at(res.id));
    }
    return res.asSegment();
}

Segment<Vertex> Embedding::embed(VertexId vid, size_t cut_left, size_t cut_right) const {
    IdSegment<VertexId> res(vid, cut_left, cut_right);
    while(vertex_embedding.find(res.id) != vertex_embedding.end()) {
        res = res.nest(vertex_embedding.at(res.id));
    }
    return res.asSegment();
}

void Embedding::remap(const GraphPath &path, Vertex &v) {
    Edge &new_edge = v == path.start() ? v.front() : *v.incoming().begin();
    size_t start = 0;
    size_t sz = path.truncLen();
    for(Edge &edge : path.edges()) {
        edge_embedding[edge.getId()] = {new_edge.getId(), start, sz - start - edge.truncSize()};
        VERIFY(edge.truncSize() == edge_embedding[edge.getId()].truncSize());
        start += edge.truncSize();
        if(edge != path.backEdge()) {
            vertex_embedding[edge.getFinish().getId()] = {v.getId(),
                                                          start + path.start().size() - edge.getFinish().size(),
                                                          sz - start};
            VERIFY(edge.getFinish().truncSize() == vertex_embedding[edge.getFinish().getId()].truncSize());
        }
    }
}

GraphPath Embedding::embed(const GraphPath &path) {
    if(!path.valid())
        return {};
    if(path.size() == 0) {
        Segment<Vertex> res(embed(path.startId(), path.cutLeft(), path.cutRight()));
        return {res.contig(), res.left, res.contig().size() - res.right};
    }
    GraphPath new_path;
    for(EdgeId eid : path.edgeIds()) {
        new_path += embed(eid);
    }
    size_t extra_cut = new_path.start().size() - embed(path.startId()).size();
    new_path.setCutLeft(path.cutLeft() + extra_cut + new_path.cutLeft());
    new_path.setCutRight(path.cutRight() + new_path.cutRight());
    return std::move(new_path);
}
