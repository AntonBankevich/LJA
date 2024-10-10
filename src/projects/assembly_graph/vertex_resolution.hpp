#pragma once
#include "assembly_graph_base.hpp"
#include "unordered_map"

namespace ag {
    template<class Traits>
    struct InOutEdgePair {
    private:
        typedef typename Traits::Edge::EdgeId EdgeId;
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;

        EdgeId first;
        EdgeId second;
    public:
        InOutEdgePair(Edge & first, Edge &second) : first(first.getId()), second(second.getId()) {
            VERIFY(this->first->getFinish() == this->second->getStart());
        }
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

    template<class Traits>
    class VertexResolutionResult {
    private:
        typedef typename Traits::Edge::EdgeId EdgeId;
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex::VertexId VertexId;
        typedef typename Traits::Vertex Vertex;

        VertexId core;
        std::unordered_map<VertexId, InOutEdgePair<Traits>> new_vertices;
        std::unordered_map<EdgeId, std::unordered_map<EdgeId, VertexId>> edge_mapping;
        void innerAdd(Vertex &new_vertex, const InOutEdgePair<Traits> &edgePair) {
            VERIFY(new_vertices.find(new_vertex.getId()) == new_vertices.end());
            new_vertices.emplace(new_vertex.getId(), edgePair);
            edge_mapping[edgePair.incoming().getId()][edgePair.outgoing().getId()] = new_vertex.getId();
        }
    public:
        VertexResolutionResult(Vertex &core) : core(core.getId()) {}
        VertexResolutionResult RC() const {
            VertexResolutionResult res(core->rc());
            for(auto it : new_vertices) {
                res.innerAdd(it.first->rc(), it.second.RC());
            }
            return std::move(res);
        }

        bool contains(Edge &edge1, Edge &edge2) const {
            VERIFY(edge1.getFinish() == *core);
            VERIFY(edge2.getStart() == *core);
            return edge_mapping.find(edge1.getId()) != edge_mapping.end() &&
                   edge_mapping.at(edge1.getId()).find(edge2.getId()) != edge_mapping.at(edge1.getId()).end();
        }

        Vertex &getCore() const {return *core;}

        Vertex &get(Edge &edge1, Edge &edge2) const {
            return *edge_mapping.at(edge1.getId()).at(edge2.getId());
        }

        const InOutEdgePair<Traits> &get(Vertex &new_vertex) const {
            return new_vertices.at(new_vertex.getId());
        }

        void add(Vertex &new_vertex, const InOutEdgePair<Traits> &edgePair) {
            innerAdd(new_vertex, edgePair);
            if(*core == core->rc() && new_vertex != new_vertex.rc()) {
                innerAdd(new_vertex.rc(), edgePair.RC());
            }
        }

        void add(Vertex &new_vertex, Edge &edge1, Edge &edge2) {add(new_vertex, {edge1, edge2});}

        IterableStorage<TransformingIterator<typename std::unordered_map<VertexId, InOutEdgePair<Traits>>::const_iterator, Vertex>> newVertices() const {
            std::function<Vertex &(const std::pair<VertexId, InOutEdgePair<Traits>> &)> transform = [](const std::pair<VertexId, InOutEdgePair<Traits>> &val) ->Vertex& {
                return *val.first;
            };
            return {{new_vertices.begin(), new_vertices.end(), transform},
                    {new_vertices.end(),   new_vertices.end(), transform}};
        }
    };

    template<class Traits>
    inline std::ostream &operator<<(std::ostream &stream, const VertexResolutionResult<Traits> &vr) {
        stream << "VRResult." << vr.getCore().getId() << ":";
        for(typename Traits::Vertex & it : vr.newVertices()) {
            stream << it.getId() << "(" << vr.get(it).incoming().getId() << "|" << vr.get(it).outgoing().getId() << ")";
        }
        return stream;
    }

}