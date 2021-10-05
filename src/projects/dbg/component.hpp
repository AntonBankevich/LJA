#pragma once

#include "paths.hpp"
#include "sparse_dbg.hpp"
#include "common/iterator_utils.hpp"
#include <utility>

namespace dbg {
    class Component {
    private:
        SparseDBG *_graph;
        std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
    public:
        template<class I>
        Component(SparseDBG &_graph, I begin, I end) : _graph(&_graph), v(begin, end) {}
        explicit Component(SparseDBG &_graph);
        template<class I>
        static Component neighbourhood(SparseDBG &graph, I begin, I end, size_t radius, size_t min_coverage = 0);
        static Component neighbourhood(SparseDBG &graph, Contig &contig, size_t radius);
        static Component longEdgeNeighbourhood(SparseDBG &graph, Contig &contig, size_t long_edge_threshold);

        SparseDBG &graph() const {return *_graph;}
        bool contains(const Vertex &vert) const {return v.find(vert.hash()) != v.end();}
        bool covers(const Vertex &vert) const;
        size_t size() const {return v.size();}

        typedef std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>>::const_iterator iterator;
        IterableStorage<ApplyingIterator<iterator, Vertex, 2>> vertices(bool unique = false) const;
        IterableStorage<ApplyingIterator<iterator, Vertex, 2>> verticesUnique(bool unique = false) const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edges(bool inner = false, bool unique = false) const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInner() const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesUnique() const;
        IterableStorage<ApplyingIterator<iterator, Edge, 16>> edgesInnerUnique() const;

        size_t countBorderEdges() const;
        size_t countTips() const;
        size_t isAcyclic() const;
        size_t realCC() const;
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
                    ConditionSplitter([min_len](const Edge& edge){return edge.size() > min_len;}){
        }
    };

    template<class I>
    dbg::Component dbg::Component::neighbourhood(dbg::SparseDBG &graph, I begin, I end, size_t radius, size_t min_coverage) {
        std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> v;
        typedef std::pair<size_t, hashing::htype> StoredValue;
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
            Vertex &vert = graph.getVertex(val.second);
            for (Edge &edge : vert) {
                if (edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
            for (Edge &edge : vert.rc()) {
                if (edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
        }
        return Component(graph, v.begin(), v.end());
    }
}