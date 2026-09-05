#pragma once

#include <alignment/alignment_form.hpp>
#include "assembly_graph_base.hpp"
#include "graph_listeners.hpp"
#include "vertex_resolution.hpp"
#include "edge_code_listener.hpp"
#include "common/output_utils.hpp"
#include "sequences/contigs.hpp"


namespace ag {
//    TODO: assume VertexData and EdgeData have default constructors. Remove them from addVertex, addEdge.
//     Fill their contents in listeners instead.
//TODO: Make marked vertex be considered deleted immidiately (e.g. run fireDeleteVertex)
//    Hidden contract: every fireX listener callback triggered by an edit below fires twice — once for
//    the edit itself, once for its reverse-complement mirror (see graph_listeners.hpp's Fire dispatchers
//    and Edge/Vertex's rc() contract note in assembly_graph_base.hpp) — EXCEPT when the edge/vertex
//    involved is self-rc (equal to its own rc()), in which case it fires only once. Before concluding a
//    listener is missing some bookkeeping, check whether the mirrored rc() call supplies it instead.
    class AssemblyGraph : public ResolutionFire {
    public:
        typedef std::list<Vertex> vertex_storage_type;
        typedef typename std::list<Vertex>::iterator vertex_iterator_type;
        typedef typename std::list<Vertex>::const_iterator const_vertex_iterator_type;

    private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
        vertex_storage_type vertex_list;
        int maxVId = 0;
        EdgeCodeListener edgeCodeListener;

        Vertex &innerAddVertex(typename Vertex::id_type id, bool canonical, VertexData data);

        Vertex &innerAddVertex(typename Vertex::id_type id, Sequence seq, VertexData data);

    protected:
        Vertex &addVertex(const Sequence &seq, const VertexData &data, typename Vertex::id_type id = 0);
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        void setHash(ag::Vertex &v, hashing::htype hash) { v.hash = hash; }
    public:

        explicit AssemblyGraph() : edgeCodeListener(*this) {}
        virtual ~AssemblyGraph();
        AssemblyGraph(AssemblyGraph &&other) = default;
        AssemblyGraph &operator=(AssemblyGraph &&other) = default;
        AssemblyGraph(const AssemblyGraph &other) noexcept = delete;

        size_t size() const {return vertex_list.size();}
        size_t edgeCount() const;

//        TODO:Rework these functions
        void removeIsolated();
        void removeMarked();
        void resetMarkers();

        Vertex &addVertexPair(VertexData data, typename Vertex::id_type id = 0);
        Vertex &addSelfRCVertex(VertexData data);
        Vertex &addVertex(const Sequence &seq, typename Vertex::id_type id = 0) {return addVertex(seq, {}, id);}
        Vertex &addVertex(const Vertex &other_graph_vertex);
        Vertex &addSPGVertex(Sequence seq, bool cyclic, bool inf_left, bool inf_right, Vertex::id_type id = Vertex::id_type());
        Edge &addSPEdgeLockFree(Vertex &start, Vertex &end, ag::Edge::id_type eid = {}, ag::Edge::id_type rcid = {});
        Edge &addSPEdge(Vertex &start, Vertex &end, ag::Edge::id_type eid = {}, ag::Edge::id_type rcid = {});
        //        This method should only be invoked if no graph modification is performed in parallel or if both start and
//        rc end vertices are locked by this process or otherwise prevented from modification by other processes
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeIdType eid = {}, EdgeIdType rcid = {});
        void removeEdgeLockFree(Edge &edge);
        void removeEdge(Edge &edge);
        void isolateAndMark(Vertex &vertex);
        template<class I>
        void isolateAndMark(I begin, I end);
//        TODO: this should be run in parallel
        void resetEdgeCodes(logging::Logger &logger, size_t threads);

        //        Make sure not to perform any other graph modifications in parallel with this method since it only blocks
//        the first and the last vertices
//        Invariant relied on by listeners (e.g. AlignedReadStatisticsTracker::fireMergePathToEdge): in a
//        Supregraph, callers only ever pass paths whose edges are uniformly all-suffix or all-prefix, never
//        a mix. mergePath enforces this by only calling mergePathToEdge on such homogeneous runs and merging
//        any other (mixed) stretch into a vertex via fireMergePath instead. This is what guarantees the
//        resulting new_edge is itself isSuffix()/isPrefix() to match the run it came from, and (by RC mirroring)
//        that the mirrored call sees the opposite, homogeneous, isPrefix()/isSuffix() run.
        Edge &mergePathToEdge(const GraphPath &path);
        Edge &mergeTipsToEdge(Edge &leftEdge, Edge &rightEdge, AlignmentForm alignment);
//        TODO: make it usable in parallel when parallel vertex adding is implemented
        GraphPath splitEdge(Edge &edge, const std::vector<EdgePosition> &split_positions);

        Edge &chooseSplitColumn(Edge &leftEdge, Edge &rightEdge, AlignmentForm alignment);

        Vertex &edgeToSupreVertex(Edge &edge);
//        Unlike mergePathToEdge/splitEdge, this is not expected to be called in parallel yet: listeners'
//        fireResolveVertex implementations are not written to be thread-safe against concurrent calls.
//        TODO: make resolveVertex (and listeners' fireResolveVertex) safe to call concurrently, the way
//        MergeAllToEdges runs mergePathToEdge in parallel across disjoint unbranching paths.
        ag::VertexResolutionResult resolveVertex(Vertex &core, const VertexResolutionPlan &resolution);
        Vertex &mergePath(const GraphPath &path);
        Vertex &mergeLoop(const GraphPath &path);


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

    inline std::string GetEdgeNameForSaving(const Edge &edge) {
        VERIFY((edge.getFinish().rc().getInnerId() > 0) == edge.getFinish().rc().isCanonical());
        if(!edge.isCanonical())
            return GetEdgeNameForSaving(edge.rc());
        return edge.getInnerId().str() + "_" + edge.rc().getInnerId().str();
    }

    template<class I>
    void AssemblyGraph::isolateAndMark(I begin, I end) {
        std::vector<VertexId> to_mark;
        while(begin != end) {
            to_mark.emplace_back(begin->getId());
            ++begin;
        }
        for(VertexId &vid : to_mark) {
            if(!vid->marked()) {
                isolateAndMark(*vid);
            }
        }
    }
}
