#pragma once

#include "supregraph_base.hpp"
#include <vector>

namespace spg {
    struct EdgePair {
        EdgeId first;
        EdgeId second;
        EdgePair(Edge & first, Edge &second) : first(first.getId()), second(second.getId()) {
            VERIFY(this->first->getFinish() == this->second->getStart());
        }
        EdgePair RC() const {return {second->rc(), first->rc()};}
        Sequence getSeq() const {return first->getStart().getSeq() + second->truncSeq();}
        Vertex &middle() const {return first->getFinish();}

        bool operator==(const EdgePair &other) const {return first == other.first && second == other.second;}
        bool operator!=(const EdgePair &other) const {return first != other.first || second != other.second;}
        bool operator<(const EdgePair &other) const {return first < other.first || (first == other.first && second < other.second);}
        bool operator>(const EdgePair &other) const {return first > other.first || (first == other.first && second > other.second);}
        bool operator<=(const EdgePair &other) const {return first <= other.first || (first == other.first && second <= other.second);}
        bool operator>=(const EdgePair &other) const {return first >= other.first || (first == other.first && second >= other.second);}
    };

    class VertexResolutionPlan {
    private:
        VertexId v;
        mutable bool sorted = true;
        mutable std::vector<EdgePair> edge_pairs;

        void sort() const;
    public:
        VertexResolutionPlan(Vertex &v) : v(v.getId()) {} // NOLINT(google-explicit-constructor)
        VertexResolutionPlan RC() const {
            VertexResolutionPlan res(v->rc());
            for(const EdgePair &ep : edge_pairs) {
                res.add(ep.RC());
            }
            return std::move(res);
        }

        Vertex &getCore() const {return *v;}
        void add(const EdgePair &edgePair);
        void add(Edge &edge1, Edge &edge2) {add({edge1, edge2});}

        bool empty() const {return edge_pairs.empty();}
        bool incConnected(Edge &edge) const;
        bool outConnected(Edge &edge) const;
        bool allConnected() const;
        IterableStorage<std::vector<EdgePair>::const_iterator> connections() const;
        IterableStorage<SkippingIterator<std::vector<EdgePair>::const_iterator>> connectionsUnique() const;
    };

    inline std::ostream &operator<<(std::ostream &stream, const VertexResolutionPlan &vr) {
        stream << "VRResult." << vr.getCore().getId() << ":";
        for(const EdgePair & it : vr.connections()) {
            stream << "(" << it.first->getId() << "|" << it.second->getId() << ")";
        }
        return stream;
    }


    class VertexResolutionResult {
    private:
        VertexId v;
        std::unordered_map<VertexId, EdgePair> new_vertices;
        std::unordered_map<EdgeId, std::unordered_map<EdgeId, VertexId>> edge_mapping;
        void innerAdd(Vertex &new_vertex, const EdgePair &edgePair) {
            VERIFY(new_vertices.find(new_vertex.getId()) == new_vertices.end());
            new_vertices.emplace(new_vertex.getId(), edgePair);
            edge_mapping[edgePair.first][edgePair.second] = new_vertex.getId();
        }
    public:
        VertexResolutionResult(Vertex &v) : v(v.getId()) {}
        VertexResolutionResult RC() const {
            VertexResolutionResult res(v->rc());
            for(auto it : new_vertices) {
                res.innerAdd(it.first->rc(), it.second.RC());
            }
            return std::move(res);
        }

        bool contains(Edge &edge1, Edge &edge2) const {
            return edge_mapping.find(edge1.getId()) != edge_mapping.end() &&
                edge_mapping.at(edge1.getId()).find(edge2.getId()) != edge_mapping.at(edge1.getId()).end();
        }

        Vertex &getVertex() const {return *v;}

        Vertex &get(Edge &edge1, Edge &edge2) const {
            return *edge_mapping.at(edge1.getId()).at(edge2.getId());
        }

        const EdgePair &get(Vertex &new_vertex) const {
            return new_vertices.at(new_vertex.getId());
        }

        void add(Vertex &new_vertex, const EdgePair &edgePair) {
            innerAdd(new_vertex, edgePair);
            if(new_vertex == new_vertex.rc()) {
                innerAdd(new_vertex.rc(), edgePair.RC());
            }
        }

        void add(Vertex &new_vertex, Edge &edge1, Edge &edge2) {add(new_vertex, {edge1, edge2});}

        IterableStorage<TransformingIterator<std::unordered_map<VertexId, EdgePair>::const_iterator, Vertex>> newVertices() const {
            std::function<Vertex &(const std::pair<VertexId, EdgePair> &)> transform = [](const std::pair<VertexId, EdgePair> &val) ->Vertex& {
                return *val.first;
            };
            return {{new_vertices.begin(), new_vertices.end(), transform},
                    {new_vertices.end(),   new_vertices.end(), transform}};
        }
    };

    inline std::ostream &operator<<(std::ostream &stream, const VertexResolutionResult &vr) {
        stream << "VRResult." << vr.getVertex().getId() << ":";
        for(Vertex & it : vr.newVertices()) {
            stream << it.getId() << "(" << vr.get(it).first->getId() << "|" << vr.get(it).second->getId() << ")";
        }
        return stream;
    }

    class DecisionRule {
    public:
        virtual VertexResolutionPlan judge(Vertex &v) = 0;

        virtual ~DecisionRule() = default;
    };

    class RandomDecisionRule : public DecisionRule {
    public:
        VertexResolutionPlan judge(Vertex &v) override {
            VertexResolutionPlan res(v);
            auto out_it = v.begin();
            auto inc = v.incoming();
            auto in_it = inc.begin();
            while (out_it != v.end() || in_it != inc.end()) {
                if (out_it == v.end()) --out_it;
                if (in_it == inc.end()) --in_it;
                res.add(*in_it, *out_it);
                ++in_it;
                ++out_it;
            }
            return std::move(res);
        }
    };
}