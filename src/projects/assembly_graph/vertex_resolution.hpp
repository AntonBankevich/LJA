#pragma once
#include "assembly_graph_base.hpp"
#include "unordered_map"

namespace ag {
    class InOutEdgePair;
    inline std::ostream &operator<<(std::ostream &stream, const InOutEdgePair &pair);

    struct InOutEdgePair {
        friend std::ostream &operator<<(std::ostream &, const InOutEdgePair &);
    private:
        EdgeId first;
        EdgeId second;
    public:
        InOutEdgePair(Edge & first, Edge &second);
        Edge & incoming() const {return *first;};
        Edge & outgoing() const {return *second;};
        InOutEdgePair RC() const {return {second->rc(), first->rc()};}
        Sequence getSeq() const {return first->getStart().getSeq() + second->truncSeq();}
        Vertex &middle() const {return first->getFinish();}

        bool operator==(const InOutEdgePair &other) const {return first == other.first && second == other.second;}
        bool operator!=(const InOutEdgePair &other) const {return first != other.first || second != other.second;}
        bool operator<(const InOutEdgePair &other) const {return first < other.first || (first == other.first && second < other.second);}
        bool operator>(const InOutEdgePair &other) const {return first > other.first || (first == other.first && second > other.second);}
        bool operator<=(const InOutEdgePair &other) const {return first <= other.first || (first == other.first && second <= other.second);}
        bool operator>=(const InOutEdgePair &other) const {return first >= other.first || (first == other.first && second >= other.second);}
    };

    std::ostream &operator<<(std::ostream &stream, const InOutEdgePair &pair);

    class VertexResolutionPlan {
    private:
        VertexId v;
        mutable bool sorted = true;
        mutable std::vector<InOutEdgePair> edge_pairs;

        void sort() const;
    public:
        VertexResolutionPlan(Vertex &v) : v(v.getId()) {} // NOLINT(google-explicit-constructor)
        VertexResolutionPlan RC() const;

        Vertex &getCore() const {return *v;}
        void add(const InOutEdgePair &edgePair);
        void add(Edge &edge1, Edge &edge2) {add({edge1, edge2});}

        bool empty() const {return edge_pairs.empty();}
        bool incConnected(Edge &edge) const;
        bool outConnected(Edge &edge) const;
        bool allConnected() const;
        IterableStorage<std::vector<InOutEdgePair>::const_iterator> connections() const;
        IterableStorage<SkippingIterator<std::vector<InOutEdgePair>::const_iterator>> connectionsUnique() const;
    };

    std::ostream &operator<<(std::ostream &stream, const VertexResolutionPlan &vr);

    class VertexResolutionResult {
    private:
        VertexId core;
        std::unordered_map<VertexId, InOutEdgePair> new_vertices;
        std::unordered_map<EdgeId, std::unordered_map<EdgeId, VertexId>> edge_mapping;
        void innerAdd(Vertex &new_vertex, const InOutEdgePair &edgePair);
    public:
        VertexResolutionResult(Vertex &core) : core(core.getId()) {} // NOLINT(google-explicit-constructor)
        VertexResolutionResult RC() const;
        bool empty() {
            return !core->marked();
        }

        bool contains(Edge &edge1, Edge &edge2) const;
        Vertex &getCore() const {return *core;}
        Vertex &get(Edge &edge1, Edge &edge2) const;
        const InOutEdgePair &get(Vertex &new_vertex) const;
        void add(Vertex &new_vertex, const InOutEdgePair &edgePair);
        void add(Vertex &new_vertex, Edge &edge1, Edge &edge2);
        IterableStorage<TransformingIterator<typename std::unordered_map<VertexId, InOutEdgePair>::const_iterator, Vertex>> newVertices() const;
        std::unordered_map<VertexId, InOutEdgePair>::const_iterator begin() const {return new_vertices.begin();}
        std::unordered_map<VertexId, InOutEdgePair>::const_iterator end() const {return new_vertices.end();}
    };

    std::ostream &operator<<(std::ostream &stream, const VertexResolutionResult &vr);
}