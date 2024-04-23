#pragma once

#include "paths.hpp"
#include "alignment_chain.hpp"
#include "common/iterator_utils.hpp"
#include <utility>
#include <unordered_set>
#include <queue>

namespace ag {
    //TODO: move to namespace ag and make it universal for all graphs!!!
    template<class Traits>
    class Component {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex::VertexId VertexId;
        typedef typename Traits::Edge::EdgeId EdgeId;
    private:
        AssemblyGraph<Traits> *_graph;
        std::unordered_set<VertexId> v;
        size_t sz = 0;

        void addVertex(VertexId vid);

    public:
        template<class I>
        Component(AssemblyGraph<Traits> &_graph, I begin, I end) {
            for (; begin != end; ++begin) {
                addVertex(*begin);
            }
        }

        explicit Component(AssemblyGraph<Traits> &_graph);

        template<class I>
        static Component
        neighbourhood(AssemblyGraph<Traits> &graph, I begin, I end, size_t radius, size_t min_coverage = 0) {
            std::unordered_set<VertexId> v;
            typedef std::pair<size_t, VertexId> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            while (begin != end) {
                queue.emplace(0, *begin);
                ++begin;
            }
            while (!queue.empty()) {
                StoredValue val = queue.top();
                queue.pop();
                if (v.find(val.second) != v.end())
                    continue;
                v.insert(val.second);
                if (val.first > radius)
                    continue;
                Vertex &vert = *val.second;
                for (Edge &edge: vert) {
                    if (edge.getCoverage() >= min_coverage)
                        queue.emplace(val.first + edge.truncSize(), edge.getFinish().getHash());
                }
                for (Edge &edge: vert.rc()) {
                    if (edge.getCoverage() >= min_coverage)
                        queue.emplace(val.first + edge.truncSize(), edge.getFinish().getHash());
                }
            }
            return {graph, v.begin(), v.end()};
        }

        static Component neighbourhood(AssemblyGraph<Traits> &graph, const std::vector<AlignmentChain<Contig, Edge>> &als1,
                      size_t radius);

        static Component longEdgeNeighbourhood(AssemblyGraph<Traits> &graph, const std::vector<AlignmentChain<Contig, Edge>> &als1,
                              size_t long_edge_threshold);

        AssemblyGraph<Traits> &graph() const { return *_graph; }

        bool contains(Vertex &vert) const { return v.find(vert.getId()) != v.end(); }

        bool covers(Vertex &vert) const;

        size_t uniqueSize() const { return sz; }

        typedef typename std::unordered_set<VertexId>::const_iterator iterator;

        IterableStorage<TransformingIterator<iterator, Vertex>> vertices() const;

        IterableStorage<SkippingIterator<TransformingIterator<iterator, Vertex>>> verticesUnique() const;

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edges(bool inner = false, bool unique = false) const;

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInner() const;

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesUnique() const;

        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInnerUnique() const;

        std::vector<Vertex *> borderVertices() const;

        size_t countBorderEdges() const;

        size_t countTips() const;

        size_t isAcyclic() const;

        size_t realCC() const;

        std::vector<Vertex *> topSort() const;
    };

    template<class Traits>
    std::vector<typename Traits::Vertex *> Component<Traits>::borderVertices() const {
        std::vector<Vertex *> res;
        for (Vertex &vert: vertices()) {
            for (Edge &edge: vert) {
                if (!contains(edge.getFinish())) {
                    res.emplace_back(&vert);
                }
            }
        }
        return std::move(res);
    }

    template<class Traits>
    size_t Component<Traits>::countBorderEdges() const {
        size_t res = 0;
        for (Vertex &vert: vertices()) {
            for (Edge &edge: vert) {
                if (!contains(edge.getFinish()))
                    res++;
            }
        }
        return res;
    }

    template<class Traits>
    size_t Component<Traits>::countTips() const {
        size_t res = 0;
        for (Vertex &vert: vertices()) {
            if (vert.inDeg() == 0)
                res++;
        }
        return res;
    }

    template<class Traits>
    size_t Component<Traits>::isAcyclic() const {
        std::unordered_set<VertexId> visited;
        for (Vertex &start: vertices()) {
            if (visited.find(start.getId()) != visited.end())
                continue;
            std::vector<VertexId> queue;
            queue.emplace_back(start.getId());
            while (!queue.empty()) {
                Vertex &vert = *queue.back();
                queue.pop_back();
                if (visited.find(vert.getId()) != visited.end())
                    continue;
                bool ok = true;
                for (Edge &edge: vert.rc()) {
                    Vertex &prev = edge.getFinish().rc();
                    if (contains(prev) && visited.find(prev.getId()) == visited.end())
                        ok = false;
                }
                if (!ok)
                    continue;
                visited.emplace(vert.getId());
                for (Edge &edge: vert) {
                    if (contains(edge.getFinish()))
                        queue.emplace_back(edge.getFinish().getId());
                }
            }
        }
        return visited.size() == v.size() * 2;
    }

    template<class Traits>
    size_t Component<Traits>::realCC() const {
        std::unordered_set<VertexId> visited;
        size_t cnt = 0;
        for (Vertex &start: vertices()) {
            if (visited.find(start.getId()) != visited.end())
                continue;
            cnt += 1;
            std::vector<VertexId> queue;
            queue.emplace_back(start.getId());
            while (!queue.empty()) {
                Vertex &vert = *queue.back();
                queue.pop_back();
                if (visited.find(vert.getId()) != visited.end())
                    continue;
                visited.emplace(vert.getId());
                for (Edge &edge: vert) {
                    if (contains(edge.getFinish()))
                        queue.emplace_back(edge.getFinish().getId());
                }
                for (Edge &edge: vert.rc()) {
                    if (contains(edge.getFinish()))
                        queue.emplace_back(edge.getFinish().rc().getId());
                }
            }
        }
        return cnt;
    }

    template<class Traits>
    IterableStorage<TransformingIterator<typename Component<Traits>::iterator, typename Traits::Vertex>>
    Component<Traits>::vertices() const {
        TransformingIterator<iterator, Vertex> begin = TransformingIterator<iterator, Vertex>::DereferencingIterator(
                v.begin(), v.end());
        TransformingIterator<iterator, Vertex> end = TransformingIterator<iterator, Vertex>::DereferencingIterator(
                v.end(), v.end());
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<SkippingIterator<TransformingIterator<typename Component<Traits>::iterator, typename Traits::Vertex>>>
    Component<Traits>::verticesUnique() const {
        std::function<bool(const Vertex &)> use = [](
                const Vertex &vert) { return vert.isCanonical(); };
        TransformingIterator<iterator, Vertex> begin0 = TransformingIterator<iterator, Vertex>::DereferencingIterator(
                v.begin(), v.end());
        TransformingIterator<iterator, Vertex> end0 = TransformingIterator<iterator, Vertex>::DereferencingIterator(
                v.end(), v.end());
        SkippingIterator<TransformingIterator<iterator, Vertex>> begin(begin0, end0, use);
        SkippingIterator<TransformingIterator<iterator, Vertex>> end(end0, end0, use);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename Component<Traits>::iterator, typename Traits::Edge, 16>>
    Component<Traits>::edges(bool inner, bool unique) const {
        std::function<std::array<Edge *, 16>(const VertexId &)> apply = [this, inner, unique](
                const VertexId &vid) {
            std::array<Edge *, 16> res = {};
            size_t cur = 0;
            for (Edge &edge: *vid) {
                bool endInner = contains(edge.getFinish());
                if (!inner && !endInner) {
                    res[cur] = &edge.rc();
                    cur++;
                }
                if ((!inner || endInner) && (edge.isCanonical() || !unique)) {
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

    template<class Traits>
    IterableStorage<ApplyingIterator<typename Component<Traits>::iterator, typename Traits::Edge, 16>>
    Component<Traits>::edgesInner() const {
        return edges(true, false);
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename Component<Traits>::iterator, typename Traits::Edge, 16>>
    Component<Traits>::edgesUnique() const {
        return edges(false, true);
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename Component<Traits>::iterator, typename Traits::Edge, 16>>
    Component<Traits>::edgesInnerUnique() const {
        return edges(true, true);
    }

    template<class Traits>
    bool Component<Traits>::covers(typename Traits::Vertex &vert) const {
        if (contains(vert))
            return true;
        for (Edge &edge: vert)
            if (!contains(edge.getFinish()))
                return false;
        for (Edge &edge: vert.rc())
            if (!contains(edge.getFinish()))
                return false;
        return true;
    }

    template<class Traits>
    Component<Traits>::Component(AssemblyGraph<Traits> &_graph) : _graph(&_graph) {
        for (auto &vert: graph().verticesUnique())
            addVertex(vert.getId());
    }

    template<class Traits>
    Component<Traits> Component<Traits>::neighbourhood(AssemblyGraph<Traits> &graph,
                                                       const std::vector<AlignmentChain<Contig, Edge>> &als1,
                                                       size_t radius) {
        std::unordered_set<VertexId> v;
        typedef std::pair<size_t, VertexId> StoredValue;
        std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
//        TODO Every edge must have a full sequence stored as a composite sequence
        for (const AlignmentChain<Contig, Edge> &al: als1) {
            if (al.seg_to.left < radius)
                queue.emplace(0, al.seg_to.contig().getStart().getId());
            VERIFY(al.seg_to.right <= al.seg_to.contig().truncSize());
            if (al.seg_to.contig().truncSize() < radius + al.seg_to.right)
                queue.emplace(0, al.seg_to.contig().getFinish().getId());
        }
        if (queue.empty()) {
            for (const AlignmentChain<Contig, Edge> &al: als1) {
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
            for (VertexId vit: {val.second, val.second->rc().getId()})
                for (Edge &edge: *vit) {
                    queue.emplace(val.first + edge.truncSize(), edge.getFinish().getId());
                }
        }
        return Component<Traits>(graph, v.begin(), v.end());
    }

    template<class Traits>
    Component<Traits> Component<Traits>::longEdgeNeighbourhood(AssemblyGraph<Traits> &graph,
                                             const std::vector<AlignmentChain<Contig, Edge>> &als1,
                                             size_t long_edge_threshold) {
        std::unordered_set<VertexId> v;
        typedef std::pair<size_t, VertexId> StoredValue;
        std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
//        TODO Every edge must have a full sequence stored as a composite sequence
        for (const AlignmentChain<Contig, Edge> &al: als1) {
            queue.emplace(0, al.seg_to.contig().getStart().getId());
            queue.emplace(0, al.seg_to.contig().getFinish().getId());
        }
        while (!queue.empty()) {
            std::pair<size_t, VertexId> val = queue.top();
            queue.pop();
            if (v.find(val.second) != v.end())
                continue;
            v.insert(val.second);
            for (VertexId vit: {val.second, val.second->rc().getId()})
                for (Edge &edge: *vit) {
                    if (edge.truncSize() < long_edge_threshold)
                        queue.emplace(val.first + edge.truncSize(), edge.getFinish().getId());
                }
        }
        return Component<Traits>(graph, v.begin(), v.end());
    }


    template<class Traits>
    std::vector<typename Traits::Vertex *> Component<Traits>::topSort() const {
//            Undefined behavior if component contains cycles.
        std::vector<Vertex *> stack;
        std::vector<Vertex *> res;
        std::unordered_set<Vertex *> visited;
        for (Vertex &vertex: vertices()) {
            bool ok = true;
            for (Edge &edge: vertex.rc()) {
                if (contains(edge.getFinish())) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                stack.emplace_back(&vertex);
            }
        }
        VERIFY(!stack.empty());
        while (!stack.empty()) {
            typename Traits::Vertex *cur = stack.back();
            if (visited.find(cur) != visited.end()) {
                stack.pop_back();
                continue;
            }
            bool ok = true;
            for (Edge &e: *cur) {
                if (contains(e.getFinish()) && visited.find(&e.getFinish()) == visited.end())
                    ok = false;
            }
            if (ok) {
                visited.emplace(cur);
                stack.pop_back();
                res.emplace_back(&cur->rc());
            } else {
                for (Edge &e: *cur) {
                    if (contains(e.getFinish()) && visited.find(&e.getFinish()) == visited.end())
                        stack.emplace_back(&e.getFinish());
                }
            }
        }
        return std::move(res);
    }

    template<class Traits>
    void Component<Traits>::addVertex(VertexId vid) {
        size_t old_size = v.size();
        v.insert(vid);
        v.insert(vid->rc().getId());
        if (v.size() > old_size) {
            sz++;
        }
    }

}