#include "component.hpp"

size_t dbg::Component::countBorderEdges() const {
    size_t res = 0;
    for (Vertex &vert : vertices()) {
        for(Edge &edge : vert) {
            if(!contains(*edge.end()))
                res++;
        }
    }
    return res;
}

size_t dbg::Component::countTips() const {
    size_t res = 0;
    for (Vertex &vert : vertices()) {
        if(vert.inDeg() == 0)
            res++;
    }
    return res;
}

size_t dbg::Component::isAcyclic() const {
    std::unordered_set<Vertex *> visited;
    for(hashing::htype hash : v) {
        Vertex & compStart = graph().getVertex(hash);
        for(Vertex *vit : {&compStart, &compStart.rc()}) {
            if(visited.find(vit) != visited.end())
                continue;
            std::vector<Vertex *> queue;
            queue.emplace_back(vit);
            while(!queue.empty()) {
                Vertex &vert = *queue.back();
                queue.pop_back();
                if(visited.find(&vert) != visited.end())
                    continue;
                bool ok = true;
                for(Edge &edge : vert.rc()) {
                    Vertex &prev = edge.end()->rc();
                    if(contains(prev) && visited.find(&prev) == visited.end())
                        ok = false;
                }
                if(!ok)
                    continue;
                visited.emplace(&vert);
                for(Edge &edge : vert) {
                    if(contains(*edge.end()))
                        queue.emplace_back(edge.end());
                }
            }
        }
    }
    return visited.size() == v.size() * 2;
}

size_t dbg::Component::realCC() const {
    std::unordered_set<Vertex *> visited;
    size_t cnt = 0;
    for(hashing::htype hash : v) {
        Vertex & compStart = graph().getVertex(hash);
        for(Vertex *vit : {&compStart, &compStart.rc()}) {
            if(visited.find(vit) != visited.end())
                continue;
            cnt += 1;
            std::vector<Vertex *> queue;
            queue.emplace_back(vit);
            while(!queue.empty()) {
                Vertex &vert = *queue.back();
                queue.pop_back();
                if(visited.find(&vert) != visited.end())
                    continue;
                visited.emplace(&vert);
                for(Edge &edge : vert) {
                    if(contains(*edge.end()))
                        queue.emplace_back(edge.end());
                }
                for(Edge &edge : vert.rc()) {
                    if(contains(*edge.end()))
                        queue.emplace_back(&edge.end()->rc());
                }
            }
        }
    }
    return cnt;
}

IterableStorage<ApplyingIterator<dbg::Component::iterator, dbg::Vertex, 2>> dbg::Component::vertices(bool unique) const {
    std::function<std::array<Vertex*, 2>(const hashing::htype &)> apply = [this, unique](const hashing::htype &hash) -> std::array<Vertex*, 2> {
        size_t cur = 0;
        Vertex &vertex = graph().getVertex(hash);
        if(unique)
            return {&graph().getVertex(hash)};
        else
            return graph().getVertices(hash);
    };
    ApplyingIterator<iterator, Vertex, 2> begin(v.begin(), v.end(), apply);
    ApplyingIterator<iterator, Vertex, 2> end(v.end(), v.end(), apply);
    return {begin, end};
}

IterableStorage<ApplyingIterator<dbg::Component::iterator, dbg::Vertex, 2>> dbg::Component::verticesUnique(bool unique) const {
    return vertices(true);
}

IterableStorage<ApplyingIterator<dbg::Component::iterator, dbg::Edge, 16>> dbg::Component::edges(bool inner, bool unique) const {
    std::function<std::array<Edge*, 16>(const hashing::htype &)> apply = [this, inner, unique](const hashing::htype &hash) {
        std::array<Edge*, 16> res = {};
        size_t cur = 0;
        for(Vertex * vertex : graph().getVertices(hash)) {
            for (Edge &edge : *vertex) {
                bool isInner = contains(*edge.end());
                if(inner && !isInner)
                    continue;
                if(!isInner || edge <= edge.rc()) {
                    res[cur] = &edge;
                    cur++;
                    if(!unique) {
                        res[cur] = &edge.rc();
                        cur++;
                    }
                }
            }
        }
        return res;
    };
    ApplyingIterator<iterator, Edge, 16> begin(v.begin(), v.end(), apply);
    ApplyingIterator<iterator, Edge, 16> end(v.end(), v.end(), apply);
    return {begin, end};
}

IterableStorage<ApplyingIterator<dbg::Component::iterator, dbg::Edge, 16>> dbg::Component::edgesInner() const {
    return edges(true, false);
}

IterableStorage<ApplyingIterator<dbg::Component::iterator, dbg::Edge, 16>> dbg::Component::edgesUnique() const {
    return edges(false, true);
}

IterableStorage<ApplyingIterator<dbg::Component::iterator, dbg::Edge, 16>> dbg::Component::edgesInnerUnique() const {
    return edges(true, true);
}

bool dbg::Component::covers(const dbg::Vertex &vert) const {
    if(contains(vert))
        return true;
    for(Edge &edge :vert)
        if(!contains(*edge.end()))
            return false;
    for(Edge &edge :vert.rc())
        if(!contains(*edge.end()))
            return false;
    return true;
}

dbg::Component::Component(dbg::SparseDBG &_graph) : _graph(&_graph) {
    for (auto &vert : graph())
        v.emplace(vert.second.hash());
}

dbg::Component dbg::Component::neighbourhood(dbg::SparseDBG &graph, Contig &contig, size_t radius) {
    std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
    typedef std::pair<size_t, hashing::htype> StoredValue;
    std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
    std::vector<PerfectAlignment<Contig, Edge>> als1 = GraphAligner(graph).carefulAlign(contig);
    Contig rc_contig = contig.RC();
//        TODO Every edge must have a full sequence stored as a composite sequence
    for (PerfectAlignment<Contig, Edge> &al : als1) {
        if(al.seg_to.left < radius)
            queue.emplace(0, al.seg_to.contig().start()->hash());
        VERIFY(al.seg_to.right <= al.seg_to.contig().size());
        if(al.seg_to.contig().size() < radius + al.seg_to.right)
            queue.emplace(0, al.seg_to.contig().end()->hash());
    }
    if(queue.empty()) {
        for (PerfectAlignment<Contig, Edge> &al : als1) {
            queue.emplace(0, al.seg_to.contig().start()->hash());
            queue.emplace(0, al.seg_to.contig().end()->hash());
        }
    }
    while (!queue.empty()) {
        std::pair<size_t, hashing::htype> val = queue.top();
        queue.pop();
        if (v.find(val.second) != v.end() || val.first > radius)
            continue;
        v.insert(val.second);
        for(Vertex *vit : graph.getVertices(val.second))
            for (Edge &edge : *vit) {
                queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
    }
    return Component(graph, v.begin(), v.end());
}

dbg::Component dbg::Component::longEdgeNeighbourhood(dbg::SparseDBG &graph, Contig &contig, size_t long_edge_threshold) {
    std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
    typedef std::pair<size_t, hashing::htype> StoredValue;
    std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
    std::vector<PerfectAlignment<Contig, Edge>> als1 = GraphAligner(graph).carefulAlign(contig);
    Contig rc_contig = contig.RC();
//        TODO Every edge must have a full sequence stored as a composite sequence
    for (PerfectAlignment<Contig, Edge> &al : als1) {
        queue.emplace(0, al.seg_to.contig().start()->hash());
        queue.emplace(0, al.seg_to.contig().end()->hash());
    }
    while (!queue.empty()) {
        std::pair<size_t, hashing::htype> val = queue.top();
        queue.pop();
        if (v.find(val.second) != v.end())
            continue;
        v.insert(val.second);
        for(Vertex *vit : graph.getVertices(val.second))
            for (Edge &edge : *vit) {
                if (edge.size() < long_edge_threshold)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
    }
    return Component(graph, v.begin(), v.end());
}

std::vector<dbg::Component> dbg::ConditionSplitter::split(const dbg::Component &comp) const {
    SparseDBG &dbg = comp.graph();
    std::vector<Component> res;
    std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> visited;
    size_t size = 0;
    for (Vertex &v : comp.verticesUnique()) {
        std::vector<hashing::htype> queue;
        if (visited.find(v.hash()) != visited.end())
            continue;
        queue.push_back(v.hash());
        std::vector<hashing::htype> component;
        while (!queue.empty()) {
            hashing::htype val = queue.back();
            queue.pop_back();
            if (visited.find(val) != visited.end())
                continue;
            visited.insert(val);
            component.emplace_back(val);
            for(Vertex *vert : dbg.getVertices(val)) {
                for (Edge &edge : *vert) {
                    if (!splitEdge(edge) && comp.contains(*edge.end())) {
                        queue.emplace_back(edge.end()->hash());
                    }
                }
            }
        }
        res.emplace_back(dbg, component.begin(), component.end());
        size += res.back().size();
    }
    VERIFY(size == comp.size());
    return std::move(res);
}
