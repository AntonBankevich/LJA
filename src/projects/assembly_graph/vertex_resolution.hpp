#pragma once
#include "assembly_graph_base.hpp"
#include "unordered_map"

namespace ag {
    class InOutEdgePair;
    inline std::ostream &operator<<(std::ostream &stream, const InOutEdgePair &pair);

    //TODO: move support out of EdgePair
    struct InOutEdgePair {
        friend std::ostream &operator<<(std::ostream &, const InOutEdgePair &);
    private:
        EdgeId first;
        EdgeId second;
        size_t support;
    public:
        InOutEdgePair(Edge & first, Edge &second, size_t support);
        Edge & incoming() const {return *first;};
        Edge & outgoing() const {return *second;};
//        Stable ids of the paired edges. Unlike incoming()/outgoing() these never dereference the
//        edges, so they stay valid map keys even after the edges have been freed (e.g. by
//        extendVertex inside ExtendAdjacent, before FixPlan re-maps the resolution plan).
        EdgeId incomingId() const {return first;};
        EdgeId outgoingId() const {return second;};
        InOutEdgePair RC() const {return {second->rc(), first->rc(), support};}
        Sequence getSeq() const {
            return first->rc().truncSeq().Subseq(0, 1).rc() + first->getFinish().getSeq() + second->truncSeq().Subseq(0, 1);
        }
        Vertex &middle() const {return first->getFinish();}
        size_t getSupport() const {return support;}
        void addSupport(const InOutEdgePair &other) {
            VERIFY(first == other.first && second == other.second);
            this->support += other.support;
        }

        bool operator==(const InOutEdgePair &other) const {return first == other.first && second == other.second;}
        bool operator!=(const InOutEdgePair &other) const {return first != other.first || second != other.second;}
        bool operator<(const InOutEdgePair &other) const {return first < other.first || (first == other.first && second < other.second);}
        bool operator>(const InOutEdgePair &other) const {return first > other.first || (first == other.first && second > other.second);}
        bool operator<=(const InOutEdgePair &other) const {return first <= other.first || (first == other.first && second <= other.second);}
        bool operator>=(const InOutEdgePair &other) const {return first >= other.first || (first == other.first && second >= other.second);}
    };

    std::ostream &operator<<(std::ostream &stream, const InOutEdgePair &pair);

    //Support value stored in VertexResolutionPlan is exactly the number of reads, supporting corresponding condition.
    //TODO: support this invariant through judgement structure: make initial judgement always the same and only allow to
    //add new pairs to it with 0 support. Make editable initialization private. Only constant empty initialization is
    //available from outside.
    class VertexResolutionPlan {
    private:
        VertexId v;
        mutable bool sorted = true;
        mutable std::vector<InOutEdgePair> edge_pairs;

        void sort() const;
    public:
        VertexResolutionPlan(Vertex &v) : v(v.getId()) {} // NOLINT(google-explicit-constructor)
        static VertexResolutionPlan SimplePlan(Vertex &v);

        VertexResolutionPlan RC() const;

        Vertex &getCore() const {return *v;}
        void add(const InOutEdgePair &edgePair);
        void add(Edge &edge1, Edge &edge2, size_t support = 0) {add({edge1, edge2, support});}

        bool empty() const {return edge_pairs.empty();}
        size_t incEdgeConnections(Edge &edge) const;
        size_t outEdgeConnections(Edge &edge) const;
        bool incConnected(Edge &edge) const;
        bool outConnected(Edge &edge) const;
        bool incConnected() const;
        bool outConnected() const;
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
//        Pairs resolved by a single new edge in.start -> out.finish instead of a new vertex.
        std::vector<std::pair<EdgeId, InOutEdgePair>> new_edges;
//        Number of resolved pairs (vertex or edge outcome) each incoming/outgoing edge of core takes part in.
        std::unordered_map<EdgeId, size_t> connection_counts;
        void innerAdd(Vertex &new_vertex, const InOutEdgePair &edgePair);
        void innerAdd(Edge &new_edge, const InOutEdgePair &edgePair);
    public:
        VertexResolutionResult(Vertex &core) : core(core.getId()) {} // NOLINT(google-explicit-constructor)
        VertexResolutionResult RC() const;
        bool empty() {
            return !core->marked();
        }

        bool contains(Edge &edge1, Edge &edge2) const;
        Vertex &getCore() const {return *core;}
        Vertex &get(Edge &edge1, Edge &edge2) const;
//        A new vertex is either a + core + b (core at [1, size - 1)), for a pair of branching edges, or an exact copy of
//        core (same size), for a 1-to-1 pair: then u -> copy is a suffix edge and copy -> w a prefix edge.
        const InOutEdgePair &get(Vertex &new_vertex) const;
        void add(Vertex &new_vertex, const InOutEdgePair &edgePair);
        void add(Edge &new_edge, const InOutEdgePair &edgePair);
        const std::vector<std::pair<EdgeId, InOutEdgePair>> &newEdges() const {return new_edges;}
        size_t connections(Edge &edge) const;
        // void add(Vertex &new_vertex, Edge &edge1, Edge &edge2);
        IterableStorage<TransformingIterator<typename std::unordered_map<VertexId, InOutEdgePair>::const_iterator, Vertex>> newVertices() const;
        std::unordered_map<VertexId, InOutEdgePair>::const_iterator begin() const {return new_vertices.begin();}
        std::unordered_map<VertexId, InOutEdgePair>::const_iterator end() const {return new_vertices.end();}
    };

    std::ostream &operator<<(std::ostream &stream, const VertexResolutionResult &vr);
}