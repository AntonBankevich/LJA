#pragma once

#include "assembly_graph/alignment_chain.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "common/iterator_utils.hpp"
#include <utility>
#include <unordered_set>
#include <queue>

namespace ag {
    //TODO: move to namespace ag and make it universal for all graphs!!!
    class Component {
    public:
        typedef typename std::unordered_set<VertexId>::const_iterator iterator;
    private:
        AssemblyGraph *_graph;
        std::unordered_set<VertexId> v;
        size_t sz = 0;
        mutable Locker<VertexId> locker = {};

        void addVertex(VertexId vid);

    public:
        Component(Component &&) = default;
        // Component(const Component &component) : _graph(component._graph), v(component.v), sz(component.sz) {}
        template<class I>
        Component(AssemblyGraph &_graph, I begin, I end);
        explicit Component(AssemblyGraph &_graph);
        template<class I>
        static Component neighbourhood(AssemblyGraph &graph, I begin, I end, size_t radius, size_t min_coverage = 0);
        static Component neighbourhood(AssemblyGraph &graph, const std::vector<AlignmentChain<Contig, Edge>> &als1,
                      size_t radius, size_t max_size = size_t(-1));
        static Component neighbourhood(AssemblyGraph &graph, const std::vector<VertexId> &vertices, size_t radius, size_t max_size);
        static Component neighbourhood(AssemblyGraph &graph, const ag::GraphPath &path, size_t radius, size_t max_size);
//        Thread safe version of neighbothood
        static Component neighbourhoodTT(AssemblyGraph &graph, const std::vector<VertexId> &vertices, size_t radius, size_t max_size);
        static Component longEdgeNeighbourhood(AssemblyGraph &graph, const std::vector<AlignmentChain<Contig, Edge>> &als1,
                              size_t long_edge_threshold, size_t max_size = size_t(-1));

        void lock() {locker = Locker<VertexId>(v.begin(), v.end());}
        void unlock() {locker = {};}


        AssemblyGraph &getGraph() const { return *_graph; }

        bool contains(Vertex &vert) const { return v.find(vert.getId()) != v.end(); }
        bool covers(Vertex &vert) const;

        size_t uniqueSize() const { return sz; }
        size_t size() const {return v.size();}
        size_t countBorderEdges() const;
        size_t countTips() const;
        size_t isAcyclic() const;
        size_t realCC() const;

        IterableStorage<TransformingIterator<iterator, Vertex>> vertices() const;
        IterableStorage<SkippingIterator<TransformingIterator<iterator, Vertex>>> verticesUnique() const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edges(bool inner = false, bool unique = false) const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInner() const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesUnique() const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInnerUnique() const;

        std::vector<VertexId> borderVertices() const;
        std::vector<VertexId> topSort() const;
    };

    template<class I>
    Component::Component(AssemblyGraph &_graph, I begin, I end) {
        for (; begin != end; ++begin) {
            addVertex(*begin);
        }
    }

    template<class I>
    Component Component::neighbourhood(AssemblyGraph &graph, I begin, I end, size_t radius, size_t min_coverage) {
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
                    queue.emplace(val.first + edge.truncSize(), edge.getFinish().getId());
            }
            for (Edge &edge: vert.rc()) {
                if (edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.truncSize(), edge.getFinish().getId());
            }
        }
        return {graph, v.begin(), v.end()};
    }
}