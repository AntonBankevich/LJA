#pragma once

#include <queue>
#include <unordered_map>
#include <common/verify.hpp>
#include <unordered_set>
#include "vector"
class Network {
public:
    struct Vertex;
    struct Edge {
        int id;
        size_t start;
        size_t end;
        size_t capacity;
        size_t min_flow;
        Edge(int id, size_t start, size_t end, size_t capacity, size_t min_capacity = 0) :
                id(id), start(start), end(end), capacity(capacity), min_flow(min_capacity) {}
    };

    struct Vertex {
        explicit Vertex(size_t _id) : id(_id) {
        }

        size_t id;
        std::vector<int> out;
        std::vector<int> inc;
    };

private:
    int innerAddEdge(size_t from, size_t to, size_t max_capacity, size_t flow = 0) {
        int eid = edges.size() + 1;
        edges.emplace_back(eid, from, to, max_capacity, flow);
        back_edges.emplace_back(-eid, to, from, 0);
        vertices[from].out.push_back(eid);
        vertices[to].inc.push_back(eid);
        vertices[to].out.push_back(-eid);
        vertices[from].inc.push_back(-eid);
        return eid;
    }

protected:
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Edge> back_edges;
    size_t source;
    size_t sink;

public:
    Edge &getEdge(int id) {
        if (id > 0)
            return edges[id - 1];
        else
            return back_edges[-id - 1];
    }

    size_t getFlow(int id) {
        return getEdge(id).min_flow + getEdge(-id).capacity;
    }

    Network() {
        source = addVertex();
        sink = addVertex();
    }

    size_t addVertex() {
        vertices.emplace_back(vertices.size());
        return vertices.size() - 1;
    }

    int addEdge(size_t from, size_t to, size_t min_capacity, size_t max_capacity) {
        if(min_capacity > 0) {
            addSource(to, min_capacity);
            addSink(from, min_capacity);
            max_capacity -= min_capacity;
        }
        return innerAddEdge(from, to, max_capacity, min_capacity);
    }

    void addSource(size_t id, size_t capacity) {
        innerAddEdge(source, id, capacity);
    }

    void addSink(size_t id, size_t capacity) {
        innerAddEdge(id, sink, capacity);
    }

private:
    std::vector<int> bfs(int startId, int endId, int avoidEdge = 0) {
        std::queue<size_t> queue;
        std::unordered_map<size_t, int> prev;
        prev[startId] = 0;
        queue.push(startId);
        while(!queue.empty()) {
            size_t next = queue.front();
            queue.pop();
            if (next == endId) {
                std::vector<int> res;
                while(next != startId) {
                    res.emplace_back(prev[next]);
                    next = getEdge(prev[next]).start;
                }
                return {res.rbegin(), res.rend()};
            }
            for(int eid : vertices[next].out) {
                Edge &edge = getEdge(eid);
                size_t end = edge.end;
                if(edge.id != avoidEdge && edge.capacity > 0 && prev.find(end) == prev.end()) {
                    prev[end] = eid;
                    queue.emplace(end);
                }
            }
        }
        return {};
    }

    void pushFlow(int edgeId, size_t val = 1) {
        VERIFY(val <= getEdge(edgeId).capacity);
        getEdge(edgeId).capacity -= val;
        getEdge(-edgeId).capacity += val;
    }

    void pushFlow(const std::vector<int> &path, size_t val = 1) {
        for(int eid : path)
            pushFlow(eid, val);
    }

    size_t outCapasity(size_t vId) {
        size_t outdeg = 0;
        for(int eid : vertices[vId].out) {
            outdeg += getEdge(eid).capacity;
        }
        return outdeg;
    }

    size_t inCapasity(size_t vId) {
        size_t indeg = 0;
        for(int eid : vertices[vId].inc) {
            indeg += getEdge(eid).capacity;
        }
        return indeg;
    }

public:
    bool fillNetwork() {
        size_t outdeg = outCapasity(source);
        for(size_t i = 0; i < outdeg; i++) {
            std::vector<int> path = bfs(source, sink);
            if(path.empty())
                return false;
            pushFlow(path, 1);
        }
        return true;
    }

    bool isInLoop(int edgeId) {
        return getEdge(edgeId).capacity > 0 && !findLoop(edgeId).empty();
    }

    std::vector<int> findLoop(int edgeId) {
        Edge edge = getEdge(edgeId);
        if(edge.capacity == 0)
            return {};
        Edge &redge = getEdge(-edgeId);
        return bfs(redge.start, redge.end, -edgeId);
    }

    size_t maxFlow(int edgeId) {
        if(isInLoop(edgeId)) {
            return getEdge(edgeId).capacity + getFlow(edgeId);
        } else {
            return getFlow(edgeId);
        }
    }

    size_t minFlow(int edgeId) {
        if(isInLoop(-edgeId)) {
            return getEdge(edgeId).min_flow;
        } else {
            return getFlow(edgeId);
        }
    }

    std::unordered_map<int, std::pair<size_t, size_t>> findBounds() {
        std::unordered_map<int, std::pair<size_t, size_t>> res;
        for(Edge &e : edges) {
            size_t minflow = minFlow(e.id);
            size_t maxflow = maxFlow(e.id);
            res[e.id] = {minflow, maxflow};
        }
        return std::move(res);
    }

    std::unordered_map<int, size_t> findFixedMultiplicities() {
        std::unordered_map<int, size_t> res;
        for(Edge &e : edges) {
            size_t minflow = minFlow(e.id);
            size_t maxflow = maxFlow(e.id);
            if(e.start != source && e.end != sink && minflow == maxflow) {
                res[e.id] = minflow;
            }
        }
        return std::move(res);
    }

    std::vector<int> findUnique() {
        std::vector<int> res;
        for(Edge &e : edges) {
            if(e.start != source && e.end != sink && maxFlow(e.id) == 1 && minFlow(e.id) == 1) {
                res.emplace_back(e.id);
            }
        }
        return std::move(res);
    }

    std::vector<size_t> topSort() {
        std::unordered_set<size_t> visited;
        std::vector<size_t> stack;
        std::vector<size_t> result;
        for(Vertex &v : vertices)
            stack.push_back(v.id);
        while(!stack.empty()) {
            size_t vid = stack.back();
            stack.pop_back();
            if(vid >= vertices.size()) {
                result.push_back(vid - vertices.size());
                continue;
            }
            if(visited.count(vid) > 0)
                continue;
            visited.emplace(vid);
            stack.push_back(vid + vertices.size());
            for(int eid : vertices[vid].out) {
                Edge & edge = getEdge(eid);
                if(edge.capacity > 0 && visited.count(edge.end) == 0) {
                    stack.push_back(edge.end);
                }
            }
        }
        return {result.rbegin(), result.rend()};
    }

    std::unordered_map<size_t, size_t> strongComponents() {
        std::vector<size_t> order = topSort();
        std::unordered_map<size_t, size_t> res;
        std::vector<std::pair<size_t, size_t>> stack;
        for(size_t vid : order)
            stack.emplace_back(vid, vid);
        while(!stack.empty()) {
            size_t vid = stack.back().first;
            size_t color = stack.back().second;
            stack.pop_back();
            if(res.find(vid) != res.end())
                continue;
            res[vid] = color;
            for(int eid : vertices[vid].inc) {
                Edge & edge = getEdge(eid);
                if(edge.capacity > 0 && res.find(edge.start) == res.end()) {
                    stack.emplace_back(edge.start, color);
                }
            }
        }
        return std::move(res);
    }
};