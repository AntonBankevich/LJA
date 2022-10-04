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
    public:
        int id;
        size_t start;
        size_t end;
        Edge(int id, size_t start, size_t end) :
                id(id), start(start), end(end) {}
        Edge(const Edge &) = delete;
        virtual bool isReal() const = 0;
        virtual size_t getMinFlow() const = 0;
        virtual size_t getFlow() const = 0;
        virtual void addFlow(int val) = 0;
        virtual size_t getRemainingCapacity() const = 0;
        virtual Edge &back() = 0;
        bool operator==(const Edge &other) const {return id == other.id;}
    };

    struct VirtualEdge;

    struct RealEdge : public Edge {
    private:
        friend class VirtualEdge;
        size_t capacity;
        size_t min_flow;
        size_t extra_flow;
        VirtualEdge *backEdge;
    public:
        RealEdge(int id, size_t start, size_t end, size_t capacity, size_t minFlow)
                : Edge(id, start, end), capacity(capacity), min_flow(minFlow), extra_flow(0), backEdge(nullptr) {}
        RealEdge(RealEdge &&other) noexcept;

        Edge &back() override {
            return *backEdge;
        }
        bool isReal() const override {
            return true;
        }
        size_t getFlow() const override {
            return min_flow + extra_flow;
        }
        size_t getMinFlow() const override {
            return min_flow;
        }
        void addFlow(int val) override {
            VERIFY(val + extra_flow >= 0);
            VERIFY(val + min_flow + extra_flow <= capacity);
            extra_flow += val;
        }
        size_t getRemainingCapacity() const override {
            return capacity - min_flow - extra_flow;
        }
    };

    struct VirtualEdge : public Edge {
    private:
        friend class RealEdge;
        RealEdge *backEdge;
    public:
        explicit VirtualEdge(RealEdge &edge) : Edge(-edge.id, edge.end, edge.start), backEdge(&edge) {
            backEdge->backEdge = this;
        }
        VirtualEdge(VirtualEdge &&other) noexcept : Edge(other.id, other.start, other.end), backEdge(other.backEdge) {
            backEdge->backEdge = this;
        }
        Edge &back() override {
            return *backEdge;
        }
        bool isReal() const override {
            return false;
        }
        size_t getFlow() const override {
            return 0;
        }
        size_t getMinFlow() const override {
            return 0;
        }
        void addFlow(int val) override {
            backEdge->addFlow(-val);
        }
        size_t getRemainingCapacity() const override {
            return backEdge->extra_flow;
        }
    };


    struct Vertex {
        explicit Vertex(size_t _id) : id(_id) {
        }

        size_t id;
        std::vector<int> out;
        std::vector<int> inc;
    };

private:
    int innerAddEdge(size_t from, size_t to, size_t max_capacity, size_t min_flow = 0) {
        int eid = edges.size() + 1;
        edges.emplace_back(eid, from, to, max_capacity, min_flow);
        back_edges.emplace_back(edges.back());
        vertices[from].out.push_back(eid);
        vertices[to].inc.push_back(eid);
        vertices[to].out.push_back(-eid);
        vertices[from].inc.push_back(-eid);
        return eid;
    }

protected:
    std::vector<Vertex> vertices;
    std::vector<RealEdge> edges;
    std::vector<VirtualEdge> back_edges;
    size_t source;
    size_t sink;

public:
    Edge &getEdge(int id) {
        if (id > 0)
            return edges[id - 1];
        else
            return back_edges[-id - 1];
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
    std::vector<int> bfs(size_t startId, size_t endId, int avoidEdge = 0) {
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
                if(edge.id != avoidEdge && edge.getRemainingCapacity() > 0 && prev.find(end) == prev.end()) {
                    prev[end] = eid;
                    queue.emplace(end);
                }
            }
        }
        return {};
    }

    void pushFlow(int edgeId, int val = 1) {
        getEdge(edgeId).addFlow(val);
    }

    void pushFlow(const std::vector<int> &path, int val = 1) {
        for(int eid : path)
            pushFlow(eid, val);
    }

    size_t outCapasity(size_t vId) {
        size_t outdeg = 0;
        for(int eid : vertices[vId].out) {
            outdeg += getEdge(eid).getRemainingCapacity();
        }
        return outdeg;
    }

    size_t inCapasity(size_t vId) {
        size_t indeg = 0;
        for(int eid : vertices[vId].inc) {
            indeg += getEdge(eid).getRemainingCapacity();
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
        return getEdge(edgeId).getRemainingCapacity() > 0 && !findLoop(edgeId).empty();
    }

    std::vector<int> findLoop(int edgeId) {
        Edge &edge = getEdge(edgeId);
        if(edge.getRemainingCapacity() == 0)
            return {};
        return bfs(edge.end, edge.start, edge.back().id);
    }

    size_t maxFlow(int edgeId) {
        Edge &edge = getEdge(edgeId);
        if(isInLoop(edgeId)) {
            return edge.getRemainingCapacity() + edge.getFlow();
        } else {
            return edge.getFlow();
        }
    }

    size_t minFlow(int edgeId) {
        Edge &edge = getEdge(edgeId);
        if(isInLoop(edge.back().id)) {
            return edge.getMinFlow();
        } else {
            return edge.getFlow();
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
};

inline Network::RealEdge::RealEdge(Network::RealEdge &&other) noexcept: Edge(other.id, other.start, other.end), capacity(other.capacity), min_flow(other.min_flow),
                                                                 extra_flow(other.extra_flow), backEdge(other.backEdge) {
    backEdge->backEdge = this;
}
