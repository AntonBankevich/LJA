#include "component.hpp"
#include "dbg_graph_aligner.hpp"

std::vector<dbg::Vertex *> dbg::Component::borderVertices() const {
    std::vector<dbg::Vertex *> res;
    for (Vertex &vert : vertices()) {
        for(Edge &edge : vert) {
            if(!contains(edge.getFinish())) {
                res.emplace_back(&vert);
            }
        }
    }
    return std::move(res);
}

size_t dbg::Component::countBorderEdges() const {
    size_t res = 0;
    for (Vertex &vert : vertices()) {
        for(Edge &edge : vert) {
            if(!contains(edge.getFinish()))
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
    std::unordered_set<VertexId> visited;
    for(Vertex &start: vertices()) {
        if(visited.find(start.getId()) != visited.end())
            continue;
        std::vector<VertexId> queue;
        queue.emplace_back(start.getId());
        while(!queue.empty()) {
            Vertex &vert = *queue.back();
            queue.pop_back();
            if(visited.find(vert.getId()) != visited.end())
                continue;
            bool ok = true;
            for(Edge &edge : vert.rc()) {
                Vertex &prev = edge.getFinish().rc();
                if(contains(prev) && visited.find(prev.getId()) == visited.end())
                    ok = false;
            }
            if(!ok)
                continue;
            visited.emplace(vert.getId());
            for(Edge &edge : vert) {
                if(contains(edge.getFinish()))
                    queue.emplace_back(edge.getFinish().getId());
            }
        }
    }
    return visited.size() == v.size() * 2;
}

size_t dbg::Component::realCC() const {
    std::unordered_set<VertexId> visited;
    size_t cnt = 0;
    for(Vertex &start : vertices()) {
        if(visited.find(start.getId()) != visited.end())
            continue;
        cnt += 1;
        std::vector<VertexId> queue;
        queue.emplace_back(start.getId());
        while(!queue.empty()) {
            Vertex &vert = *queue.back();
            queue.pop_back();
            if(visited.find(vert.getId()) != visited.end())
                continue;
            visited.emplace(vert.getId());
            for(Edge &edge : vert) {
                if(contains(edge.getFinish()))
                    queue.emplace_back(edge.getFinish().getId());
            }
            for(Edge &edge : vert.rc()) {
                if(contains(edge.getFinish()))
                    queue.emplace_back(edge.getFinish().rc().getId());
            }
        }
    }
    return cnt;
}

IterableStorage<TransformingIterator<dbg::Component::iterator, dbg::Vertex>> dbg::Component::vertices() const {
    TransformingIterator<iterator, Vertex> begin = TransformingIterator<iterator, Vertex>::DereferencingIterator(v.begin(), v.end());
    TransformingIterator<iterator, Vertex> end = TransformingIterator<iterator, Vertex>::DereferencingIterator(v.end(), v.end());
    return {begin, end};
}

IterableStorage<SkippingIterator<TransformingIterator<dbg::Component::iterator, dbg::Vertex>>> dbg::Component::verticesUnique() const {
    std::function<bool(const Vertex &)> use = [](const Vertex &vert) {return vert.isCanonical();};
    TransformingIterator<iterator, Vertex> begin0 = TransformingIterator<iterator, Vertex>::DereferencingIterator(v.begin(), v.end());
    TransformingIterator<iterator, Vertex> end0 = TransformingIterator<iterator, Vertex>::DereferencingIterator(v.end(), v.end());
    SkippingIterator<TransformingIterator<iterator, Vertex>> begin(begin0, end0, use);
    SkippingIterator<TransformingIterator<iterator, Vertex>> end(end0, end0, use);
    return {begin, end};
}

IterableStorage<ApplyingIterator<dbg::Component::iterator, dbg::Edge, 16>> dbg::Component::edges(bool inner, bool unique) const {
    std::function<std::array<Edge*, 16>(const VertexId &)> apply = [this, inner, unique](const VertexId &vid) {
        std::array<Edge*, 16> res = {};
        size_t cur = 0;
        for (Edge &edge : *vid) {
            bool endInner = contains(edge.getFinish());
            if(!inner && !endInner) {
                res[cur] = &edge.rc();
                cur++;
            }
            if((!inner || endInner) && (edge.isCanonical() || !unique)) {
                res[cur] = &edge;
                cur++;
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

bool dbg::Component::covers(dbg::Vertex &vert) const {
    if(contains(vert))
        return true;
    for(Edge &edge :vert)
        if(!contains(edge.getFinish()))
            return false;
    for(Edge &edge :vert.rc())
        if(!contains(edge.getFinish()))
            return false;
    return true;
}

dbg::Component::Component(dbg::SparseDBG &_graph) : _graph(&_graph) {
    for (auto &vert : graph().verticesUnique())
        addVertex(vert.getId());
}

dbg::Component dbg::Component::neighbourhood(dbg::SparseDBG &graph, const std::vector<PerfectAlignment<Contig, Edge>> &als1, size_t radius) {
    std::unordered_set<VertexId> v;
    typedef std::pair<size_t, VertexId> StoredValue;
    std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
//        TODO Every edge must have a full sequence stored as a composite sequence
    for (const PerfectAlignment<Contig, Edge> &al : als1) {
        if(al.seg_to.left < radius)
            queue.emplace(0, al.seg_to.contig().getStart().getId());
        VERIFY(al.seg_to.right <= al.seg_to.contig().truncSize());
        if(al.seg_to.contig().truncSize() < radius + al.seg_to.right)
            queue.emplace(0, al.seg_to.contig().getFinish().getId());
    }
    if(queue.empty()) {
        for (const PerfectAlignment<Contig, Edge> &al : als1) {
            queue.emplace(0, al.seg_to.contig().getStart().getId());
            queue.emplace(0, al.seg_to.contig().getFinish().getId());
        }
    }
    while (!queue.empty()) {
        std::pair<size_t, VertexId> val = queue.top();
        queue.pop();
        if (v.find(val.second) != v.end() || val.first > radius)
            continue;
        v.insert(val.second);
        for(VertexId vit : {val.second, val.second->rc().getId()})
            for (Edge &edge : *vit) {
                queue.emplace(val.first + edge.truncSize(), edge.getFinish().getId());
            }
    }
    return Component(graph, v.begin(), v.end());
}

dbg::Component dbg::Component::longEdgeNeighbourhood(dbg::SparseDBG &graph, const std::vector<PerfectAlignment<Contig, Edge>> &als1, size_t long_edge_threshold) {
    std::unordered_set<VertexId> v;
    typedef std::pair<size_t, VertexId> StoredValue;
    std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
//        TODO Every edge must have a full sequence stored as a composite sequence
    for (const PerfectAlignment<Contig, Edge> &al : als1) {
        queue.emplace(0, al.seg_to.contig().getStart().getId());
        queue.emplace(0, al.seg_to.contig().getFinish().getId());
    }
    while (!queue.empty()) {
        std::pair<size_t, VertexId> val = queue.top();
        queue.pop();
        if (v.find(val.second) != v.end())
            continue;
        v.insert(val.second);
        for(VertexId vit : {val.second, val.second->rc().getId()})
            for (Edge &edge : *vit) {
                if (edge.truncSize() < long_edge_threshold)
                    queue.emplace(val.first + edge.truncSize(), edge.getFinish().getId());
            }
    }
    return Component(graph, v.begin(), v.end());
}

std::vector<dbg::Component> dbg::ConditionSplitter::split(const dbg::Component &comp) const {
    SparseDBG &dbg = comp.graph();
    std::vector<Component> res;
    std::unordered_set<VertexId> visited;
    size_t size = 0;
    for (Vertex &v : comp.verticesUnique()) {
        std::vector<VertexId> queue;
        if (visited.find(v.getId()) != visited.end())
            continue;
        queue.push_back(v.getId());
        queue.push_back(v.rc().getId());
        std::vector<VertexId> component;
        while (!queue.empty()) {
            VertexId vert = queue.back();
            queue.pop_back();
            if (visited.find(vert) != visited.end())
                continue;
            visited.insert(vert);
            component.emplace_back(vert);
            for (Edge &edge : *vert) {
                if (!splitEdge(edge) && comp.contains(edge.getFinish())) {
                    queue.emplace_back(edge.getFinish().getId());
                    queue.emplace_back(edge.getFinish().rc().getId());
                }
            }
        }
        res.emplace_back(dbg, component.begin(), component.end());
        size += res.back().uniqueSize();
    }
    VERIFY(size == comp.uniqueSize());
    return std::move(res);
}

std::vector<dbg::Vertex *> dbg::Component::topSort() const {
//            Undefined behavior if component contains cycles.
    std::vector<Vertex *> stack;
    std::vector<Vertex *> res;
    std::unordered_set<Vertex *> visited;
    for(Vertex &v : vertices()) {
        bool ok = true;
        for(Edge &edge : v.rc()) {
            if(contains(edge.getFinish())) {
                ok = false;
                break;
            }
        }
        if(ok) {
            stack.emplace_back(&v);
        }
    }
    VERIFY(!stack.empty());
    while(!stack.empty()) {
        Vertex *cur = stack.back();
        if(visited.find(cur) != visited.end()) {
            stack.pop_back();
            continue;
        }
        bool ok = true;
        for(Edge &e : *cur) {
            if(contains(e.getFinish()) && visited.find(&e.getFinish()) == visited.end())
                ok = false;
        }
        if(ok) {
            visited.emplace(cur);
            stack.pop_back();
            res.emplace_back(&cur->rc());
        } else {
            for(Edge &e : *cur) {
                if(contains(e.getFinish()) && visited.find(&e.getFinish()) == visited.end())
                    stack.emplace_back(&e.getFinish());
            }
        }
    }
    return std::move(res);
}

void dbg::Component::addVertex(dbg::VertexId vid) {
    size_t old_size = v.size();
    v.insert(vid);
    v.insert(vid->rc().getId());
    if(v.size() > old_size) {
        sz++;
    }
}
