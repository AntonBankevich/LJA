#pragma once

#include "dbg_graph_aligner.hpp"
#include "paths.hpp"
#include "sparse_dbg.hpp"
#include "common/iterator_utils.hpp"
#include <utility>

namespace dbg {
    //TODO: move to namespace ag and make it universal for all graphs!!!
    class Component {
    private:
        SparseDBG *_graph;
        std::unordered_set<VertexId> v;
        size_t sz = 0;
        void addVertex(VertexId vid);
    public:
        template<class I>
        Component(SparseDBG &_graph, I begin, I end);
        explicit Component(SparseDBG &_graph);

        template<class I>
        static Component neighbourhood(SparseDBG &graph, I begin, I end, size_t radius, size_t min_coverage = 0);
        static Component neighbourhood(dbg::SparseDBG &graph, const std::vector<PerfectAlignment<Contig, Edge>> &als1, size_t radius);
        static Component longEdgeNeighbourhood(SparseDBG &graph, const std::vector<PerfectAlignment<Contig, Edge>> &als1, size_t long_edge_threshold);

        SparseDBG &graph() const {return *_graph;}
        bool contains(Vertex &vert) const {return v.find(vert.getId()) != v.end();}
        bool covers(Vertex &vert) const;
        size_t uniqueSize() const {return sz;}

        typedef std::unordered_set<VertexId>::const_iterator iterator;
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

    class AbstractSplitter {
    public:
        virtual std::vector<Component> split(const Component &component) const = 0;
        std::vector<Component> splitGraph(SparseDBG &dbg) const {return split(Component(dbg));}
    };

    class ConditionSplitter : public AbstractSplitter {
    private:
        std::function<bool(const Edge &)> splitEdge;
    public:
        explicit ConditionSplitter(std::function<bool(const Edge &)> splitEdge) : splitEdge(std::move(splitEdge)) {}
        std::vector<Component> split(const Component &comp) const override;
    };

    class CCSplitter : public ConditionSplitter {
    private:
        std::function<bool(const Edge &)> splitEdge;
    public:
        explicit CCSplitter() : ConditionSplitter([](const Edge &){return false;}) {}
    };

    class LengthSplitter : public ConditionSplitter {
    public:
        explicit LengthSplitter(size_t min_len) :
                    ConditionSplitter([min_len](const Edge& edge){return edge.fullSize() > min_len;}){
        }
    };

    template<class I>
    dbg::Component dbg::Component::neighbourhood(dbg::SparseDBG &graph, I begin, I end, size_t radius, size_t min_coverage) {
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
            for (Edge &edge : vert) {
                if (edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.truncSize(), edge.getFinish().getHash());
            }
            for (Edge &edge : vert.rc()) {
                if (edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.truncSize(), edge.getFinish().getHash());
            }
        }
        return {graph, v.begin(), v.end()};
    }

    template<class I>
    Component::Component(SparseDBG &_graph, I begin, I end) : _graph(&_graph) {
        for(;begin != end; ++begin) {
            addVertex(*begin);
        }
    }

}