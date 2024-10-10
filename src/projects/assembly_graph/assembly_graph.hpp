#pragma once

#include "assembly_graph_base.hpp"
#include "listeners.hpp"
#include "vertex_resolution.hpp"
#include "sequences/contigs.hpp"
#include "paths.hpp"

namespace ag {
//    TODO: assume VertexData and EdgeData have default constructors. Remove them from addVertex, addEdge.
//     Fill their contents in lesteners instead.
//TODO: Make marked vertex be considered deleted immidiately (e.g. run fireDeleteVertex)
    template<class Traits>
    class AssemblyGraph : public ResolutionFire<Traits> {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Traits::VertexData VertexData;
        typedef typename Traits::EdgeData EdgeData;
        typedef std::list<Vertex> vertex_storage_type;
        typedef typename std::list<Vertex>::iterator vertex_iterator_type;
        typedef typename std::list<Vertex>::const_iterator const_vertex_iterator_type;

    private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
        vertex_storage_type vertex_list;
        int maxVId = 0;

//    Be careful since hash does not define vertex. Rc vertices share the same hash
        Vertex &innerAddVertex(typename Vertex::id_type id, bool canonical, VertexData data);

        Vertex &innerAddVertex(typename Vertex::id_type id, Sequence seq, VertexData data);

    public:

        explicit AssemblyGraph() = default;
        virtual ~AssemblyGraph();
        AssemblyGraph(AssemblyGraph &&other) = default;
        AssemblyGraph &operator=(AssemblyGraph &&other) = default;
        AssemblyGraph(const AssemblyGraph &other) noexcept = delete;
        void fillFrom(const AssemblyGraph<Traits> &other);

        size_t size() const {return vertex_list.size();}
        size_t edgeCount() const;

        void removeIsolated();
        void removeMarked();
        void resetMarkers();

        Vertex &addVertexPair(const VertexData &data, typename Vertex::id_type id = 0);
        Vertex &addSelfRCVertex(VertexData data);
        Vertex &addVertex(const Sequence &seq, const VertexData &data = {}, typename Vertex::id_type id = 0);
        Vertex &addVertex(const Sequence &seq, typename Vertex::id_type id);
        Vertex &addVertex(const Vertex &other_graph_vertex);
        //        This method should only be invoked if no graph modification is performed in parallel or if both start and
//        rc end vertices are locked by this process or otherwise prevented from modification by other processes
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {});
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {});
        void removeEdgeLockFree(Edge &edge);
        void removeEdge(Edge &edge);
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {});
        void isolateAndMark(Vertex &vertex);
//        Make sure not to perform any other graph modifications in parallel with this method since it only blocks
//        the first and the last vertices
        Edge &mergePathToEdge(const GraphPath<Traits> &path);


        IterableStorage<SkippingIterator<AssemblyGraph::vertex_iterator_type>> vertices(bool unique = false) &;
        IterableStorage<SkippingIterator<AssemblyGraph::vertex_iterator_type>> vertices(bool unique = false) && = delete;
        IterableStorage<SkippingIterator<AssemblyGraph::const_vertex_iterator_type>> vertices(bool unique = false) const &;
        IterableStorage<SkippingIterator<AssemblyGraph::const_vertex_iterator_type>> vertices(bool unique = false) const && = delete;
        IterableStorage<SkippingIterator<AssemblyGraph::vertex_iterator_type>> verticesUnique() &;
        IterableStorage<SkippingIterator<AssemblyGraph::vertex_iterator_type>> verticesUnique() && = delete;
        IterableStorage<SkippingIterator<AssemblyGraph::const_vertex_iterator_type>> verticesUnique() const &;
        IterableStorage<SkippingIterator<AssemblyGraph::const_vertex_iterator_type>> verticesUnique() const && = delete;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 4>> edges(bool unique = false) &;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 4>> edges(bool unique = false) && = delete;
        IterableStorage<ApplyingIterator<const_vertex_iterator_type, const Edge, 4>> edges(bool unique = false) const &;
        IterableStorage<ApplyingIterator<const_vertex_iterator_type, const Edge, 4>> edges(bool unique = false) const && = delete;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 4>> edgesUnique() &;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 4>> edgesUnique() && = delete;
        IterableStorage<ApplyingIterator<const_vertex_iterator_type, const Edge, 4>> edgesUnique() const &;
        IterableStorage<ApplyingIterator<const_vertex_iterator_type, const Edge, 4>> edgesUnique() const && = delete;
    };

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addSelfRCVertex(VertexData data) {
        typename Vertex::id_type id = maxVId + 1;
        Vertex &res = innerAddVertex(id, true, std::move(data));
        res.setRC(res);
        fireAddVertex(res);
        return res;
    }
    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertexPair(const VertexData &data, typename Vertex::id_type id) {
        VERIFY(id >= 0);
        if(id == 0)
            id = maxVId + 1;
        Vertex &rc = innerAddVertex(-id, false, data.RC());
        Vertex &res = innerAddVertex(id, true, std::move(data));
        res.setRC(rc);
        this->fireAddVertex(res);
        this->fireAddVertex(rc);
        return res;
    }

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertex(const Vertex &other_graph_vertex) {
        typename AssemblyGraph<Traits>::Vertex &res = addVertex(other_graph_vertex.getSeq(), other_graph_vertex, other_graph_vertex.getInnerId());
        res.updateMaxOutId(other_graph_vertex.getMaxOutId());
        res.rc().updateMaxOutId(other_graph_vertex.rc().getMaxOutId());
        return res;
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type>> AssemblyGraph<Traits>::vertices(bool unique) & {
        std::function<bool(Vertex &)> use =
                [unique](Vertex &vertex) -> bool {
                    return !unique || vertex.isCanonical();
                };
        SkippingIterator<vertex_iterator_type> begin(vertex_list.begin(), vertex_list.end(), use);
        SkippingIterator<vertex_iterator_type> end(vertex_list.end(), vertex_list.end(), use);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type>> AssemblyGraph<Traits>::vertices(bool unique) const & {
        std::function<bool(const Vertex &)> use =
                [unique](const Vertex &vertex) -> bool {
                    return !unique || vertex.isCanonical();
                };
        SkippingIterator<const_vertex_iterator_type> begin(vertex_list.begin(), vertex_list.end(), use);
        SkippingIterator<const_vertex_iterator_type> end(vertex_list.end(), vertex_list.end(), use);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type>> AssemblyGraph<Traits>::verticesUnique() &{
        return vertices(true);
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type>> AssemblyGraph<Traits>::verticesUnique() const &{
        return vertices(true);
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type, typename Traits::Edge, 4>> AssemblyGraph<Traits>::edges(bool unique) & {
        std::function<std::array<Edge*, 4>(Vertex &)> apply = [unique](Vertex &vertex) {
            if(vertex.outDeg() > 4) {
                std::cout << vertex.getSeq() << std::endl;
                for(Edge &e : vertex) {
                    std::cout << e.truncSeq() << std::endl;
                }
            }
            VERIFY(vertex.outDeg() <= 4);
            std::array<Edge*, 4> res = {};
            size_t cur = 0;
            for(Edge &edge : vertex) {
                if(!unique || edge <= edge.rc()) {
                    res[cur] = &edge;
                    cur++;
                }
            }
            return res;
        };
        ApplyingIterator<vertex_iterator_type, Edge, 4> begin(vertex_list.begin(), vertex_list.end(), apply);
        ApplyingIterator<vertex_iterator_type, Edge, 4> end(vertex_list.end(), vertex_list.end(), apply);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type, const typename Traits::Edge, 4>> AssemblyGraph<Traits>::edges(bool unique) const & {
        std::function<std::array<const Edge*, 4>(const Vertex &)> apply = [unique](const Vertex &vertex) {
            std::array<const Edge*, 4> res = {};
            size_t cur = 0;
            for(const Edge &edge : vertex) {
                if(!unique || edge <= edge.rc()) {
                    res[cur] = &edge;
                    cur++;
                }
            }
            return res;
        };
        ApplyingIterator<const_vertex_iterator_type, const Edge, 4> begin(vertex_list.begin(), vertex_list.end(), apply);
        ApplyingIterator<const_vertex_iterator_type, const Edge, 4> end(vertex_list.end(), vertex_list.end(), apply);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type, typename Traits::Edge, 4>> AssemblyGraph<Traits>::edgesUnique() &{
        return edges(true);
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type, const typename Traits::Edge, 4>> AssemblyGraph<Traits>::edgesUnique() const &{
        return edges(true);
    }

    template<class Traits>
    void AssemblyGraph<Traits>::removeIsolated() {
        vertex_storage_type newv;
        for (auto it = vertex_list.begin(); it != vertex_list.end();) {
            if (it->outDeg() == 0 && it->inDeg() == 0) {
                this->fireDeleteVertex(*it);
                it = vertex_list.erase(it);
            } else {
                ++it;
            }
        }
    }

    template<class Traits>
    void AssemblyGraph<Traits>::removeMarked() {
        for (auto it = vertex_list.begin(); it != vertex_list.end();) {
            if (it->marked()) {
                this->fireDeleteVertex(*it);
                it = vertex_list.erase(it);
            } else {
                ++it;
            }
        }
    }

    template<class Traits>
    void AssemblyGraph<Traits>::resetMarkers() {
        for(Vertex &vertex : vertices()) {
            for(Edge &edge : vertex) {
                edge.mark(EdgeMarker::common);
            }
            vertex.unmark();
        }
    }

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertex(const Sequence &seq, const VertexData &data, typename Vertex::id_type id) {
        if (id == 0) {
            if(seq <= !seq)
                id = maxVId + 1;
            else
                id = -maxVId - 1;
        }
        maxVId = std::max(std::abs(id), maxVId);
        Vertex & res = innerAddVertex(id, seq, data);
        Vertex & rc = seq == !seq ? res : innerAddVertex(-id, !seq, data.RC());
        res.setRC(rc);
        res.setSeq(seq);
        this->fireAddVertex(res);
        if(res != rc)
            this->fireAddVertex(rc);
        return res;
    }

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertex(const Sequence &seq, typename Vertex::id_type id) {
        return addVertex(seq, VertexData(), id);
    }

    template<class Traits>
    size_t AssemblyGraph<Traits>::edgeCount() const {
        size_t res = 0;
        for(auto &v : vertices())
            res += v.outDeg();
        return res;
    }

    template<class Traits>
    void AssemblyGraph<Traits>::fillFrom(const AssemblyGraph<Traits> &other) {
        for(const Vertex &v : other.verticesUnique()) {
            addVertex(v);
        }
        IdIndex<Vertex> index(vertices().begin(), vertices().end());
        for(const Edge &e : other.edgesUnique()) {
            index.getById(e.getStart().getInnerId()).addEdgeLockFree(index.getById(e.getFinish().getInnerId()), e.getSeq(), e, e.getInnerId(), e.rc().getInnerId());
        }
        for(Edge &edge : edges()) {
            edge.incCov(-edge.intCov());
        }
    }

    template<class Traits>
    AssemblyGraph<Traits>::~AssemblyGraph() {
        for(Vertex &v : verticesUnique())
            if(!v.marked())
                isolateAndMark(v);
        removeMarked();
    }

    template<class Traits>
    typename Traits::Vertex &AssemblyGraph<Traits>::innerAddVertex(typename Vertex::id_type id, bool canonical, VertexData data) {
        VERIFY(canonical == (id > 0));
        maxVId = std::max(std::abs(id), maxVId);
        vertex_list.emplace_back(id, canonical, std::move(data));
        return vertex_list.back();
    }

    template<class Traits>
    typename Traits::Vertex &AssemblyGraph<Traits>::innerAddVertex(typename Vertex::id_type id, Sequence seq, VertexData data) {
        VERIFY(seq.isCanonical() == (id > 0));
        maxVId = std::max(std::abs(id), maxVId);
        vertex_list.emplace_back(id, std::move(seq), std::move(data));
        return vertex_list.back();
    }

    template<class Traits>
    typename Traits::Edge &
    AssemblyGraph<Traits>::addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeData data,
                                           BaseEdgeId eid, BaseEdgeId rcid) {
        return addEdgeLockFree(start, end, full_sequence.Subseq(start.size()), full_sequence.rc().Subseq(end.size()), data, eid, rcid);
    }

    template<class Traits>
    typename Traits::Edge &AssemblyGraph<Traits>::addEdgeLockFree(Vertex &start, Vertex &end,
                                          const Sequence &tseq, const Sequence &rctseq,
                                          EdgeData data, BaseEdgeId eid, BaseEdgeId rcid) {
        for(Edge &edge: start) {
            if(edge.getFinish() == end && edge.truncSeq() == tseq) {
                return edge;
            }
        }
        Edge &res = start.innerAddEdge(end, tseq, data, eid);
        if(end.rc() != start || (tseq.size() > start.size() && tseq.rc().Subseq(end.size()) != rctseq.rc().Subseq(start.size()))) {
            Edge &rc_edge = end.rc().innerAddEdge(start.rc(), rctseq, data.RC(), rcid);
            res._rc = &rc_edge;
            rc_edge._rc = &res;
        } else {
            res._rc = &res;
        }
        VERIFY(res.fullSize() == res.rc().fullSize());
        this->fireAddEdge(res);
        if(res != res.rc())
            this->fireAddEdge(res.rc());
        return res;

        return addEdgeLockFree(start, end, tseq, rctseq, data, eid, rcid);
    }

    template<class Traits>
    void AssemblyGraph<Traits>::removeEdgeLockFree(Edge &edge) {
        this->fireDeleteEdge(edge);
        if(edge != edge.rc())
            this->fireDeleteEdge(edge.rc());
        if(edge != edge.rc())
            edge.getFinish().rc().innerRemoveEdge(edge.rc());
        edge.getStart().innerRemoveEdge(edge);
    }

    template<class Traits>
    void AssemblyGraph<Traits>::removeEdge(Edge &edge) {
        Locker<VertexId> locker({edge.getStart().getId(), edge.getFinish().rc().getId()});
        removeEdgeLockFree(edge);
    }

    template<class Traits>
    typename Traits::Edge &
    AssemblyGraph<Traits>::addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeData data, BaseEdgeId eid,
                                   BaseEdgeId rcid) {
        return addEdge(start, end, full_seq.Subseq(start.size()), full_seq.rc().Subseq(end.size()), data, eid, rcid);
    }

    template<class Traits>
    typename Traits::Edge &AssemblyGraph<Traits>::addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq,
                                         EdgeData data, BaseEdgeId eid, BaseEdgeId rcid) {
        Locker<VertexId> locker({start.getId(), end.rc().getId()});
        return addEdgeLockFree(start, end, tseq, rctseq, std::move(data), eid, rcid);
    }

    template<class Traits>
    void AssemblyGraph<Traits>::isolateAndMark(Vertex &vertex) {
        VERIFY(!vertex.marked());
        for(Vertex &v : ThisAndRC(vertex)) {
            while (v.outDeg() != 0) {
                removeEdgeLockFree(v.front());
            }
            v.mark();
        }
    }

    template<class Traits>
    typename Traits::Edge &AssemblyGraph<Traits>::mergePathToEdge(const GraphPath<Traits> &path) {
        VERIFY(!path.empty());
        VERIFY(path.endClosed() && path.startClosed());
        if(path.size() == 1) {
//                TODO: Make this VERIFY rather then check
            return path.frontEdge();
        }
        ag::Locker<typename Traits::Vertex::VertexId> locker({path.start().getId(), path.finish().getId()});
        for(size_t i = 1; i + 1 < path.size(); i++) {
            VERIFY(!path.getVertex(i).marked());
        }
        VERIFY(path.start() == path.finish() || path.start() == path.finish().rc() || (path.start().isJunction() && path.finish().isJunction()));
        SequenceBuilder sb;
        Sequence new_seq = path.Seq();
        EdgeData new_data = Traits::EdgeData::Merge(path.edges().begin(), path.edges().end());
        Edge &new_edge = addEdgeLockFree(path.start(), path.finish(), new_seq, new_data);
        VERIFY((path == path.RC()) == (new_edge == new_edge.rc()));
        this->fireMergePathToEdge(path, new_edge);
        if(new_edge != new_edge.rc())
            this->fireMergePathToEdge(path.RC(), new_edge.rc());

        std::vector<typename Traits::Vertex::VertexId> inner_vertices;
        for(size_t i = 1; i < path.size(); i++){
            inner_vertices.emplace_back(path.getVertex(i).getId());
        }
        for(typename Traits::Vertex::VertexId &vid : inner_vertices) {
            if(!vid->marked())
                isolateAndMark(*vid);
        }
        return new_edge;
    }
}
