#pragma once

#include <unordered_set>
#include <alignment/alignment_form.hpp>
#include "assembly_graph_base.hpp"
#include "graph_listeners.hpp"
#include "vertex_resolution.hpp"
#include "edge_code_listener.hpp"
#include "sample_segment_tracker.hpp"
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
//        k of the DBG the graph was built from (0 if it was not built from a DBG). It is also the minimal size of a
//        sample segment and the size of the sample segments when they are bundled together in a chunk.
        size_t k = 0;
//        The two most basic listeners, in this order: they react to every operation before any other listener, so
//        the others can rely on edge codes and sample segments being up to date.
        EdgeCodeListener edgeCodeListener;
        SampleSegmentTracker sampleSegmentTracker;

        Vertex &innerAddVertex(typename Vertex::id_type id, bool canonical, VertexData data);

        Vertex &innerAddVertex(typename Vertex::id_type id, Sequence seq, VertexData data);

    protected:
        void setHash(ag::Vertex &v, hashing::htype hash) { v.hash = hash; }
    public:

        explicit AssemblyGraph(size_t k = 0) : k(k), edgeCodeListener(*this), sampleSegmentTracker(*this) {}
        size_t getK() const {return k;}
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

//        Checks every structural invariant of a graph at rest (no operation in progress): rc symmetry of vertices,
//        edges and their sequences; every edge spells start..finish; ids and ownership; the splitting property (an
//        outgoing suffix edge is a vertex's only outgoing edge, otherwise outgoing edges have distinct first
//        nucleotides); all edges corporeal; marked vertices isolated; vertex_list sorted by id_cmp. If the structure is
//        consistent, also runs every listener's checkConsistency (e.g. EdgeCodeListener checks edge codes). Reports
//        every violation to os and returns false if any is found. Graph types with further structural invariants (e.g.
//        Supregraph) extend it.
        virtual bool checkConsistency(std::ostream &os = std::cerr);

        // This method resorves a chunk of Ids. In this case management of ids in thes segment is delegated to whomever
        // requested this chunk. Imprtant for reproducibility of results.
        // Ides reserved this way can can be used as parameters of addVertex.
        // TODO: make all calls of addVertex with specified id only use reserved ids.
        Vertex::id_type reserveVertexIds(size_t num);

        Vertex &addVertexPair(VertexData data, typename Vertex::id_type id = 0);
        Vertex &addSelfRCVertex(VertexData data);
        Vertex &addVertex(const Sequence &seq, const VertexData &data, typename Vertex::id_type id = 0);
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        void removeEdgeLockFree(Edge &edge);
        void resetEdgeCodes(logging::Logger &logger, size_t threads);//        TODO: this should be run in parallel
        //DBG specific
        Edge &mergePathToEdgeLockFree(const GraphPath &path);
        Edge &mergeTipsToEdge(Edge &leftEdge, Edge &rightEdge, AlignmentForm alignment);
        Vertex &mergePathLockFree(const GraphPath &path, Vertex::id_type id = 0);
        GraphPath splitEdge(Edge &edge, const std::vector<EdgePosition> &split_positions);//        TODO: this should run in parallel
        //SPG specific
//        Merges a path of only prefix edges into a single prefix edge start->finish. A path of only suffix edges is
//        merged by running on its RC; the returned edge always goes in the direction of the given path. Inner
//        vertices are isolated and marked.
        Edge &mergePrefixPath(const GraphPath &path);
        Vertex &edgeToSupreVertex(Edge &edge);
        ag::VertexResolutionResult resolveCoreVertex(Vertex &core, const VertexResolutionPlan &resolution); //TODO:: this should run in parallel
        Vertex &mergeLoop(const GraphPath &path); //TOOD: implement
        Vertex &flipVertexLockFree(Vertex &core);

//-------------------Derived graph operations. These call core graph operations.-------------------------
        Vertex &addVertex(const Sequence &seq, typename Vertex::id_type id = Vertex::id_type()) {return addVertex(seq, {}, id);}
        Vertex &addVertex(const Vertex &other_graph_vertex);
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addSPEdgeLockFree(Vertex &start, Vertex &end, ag::Edge::id_type eid = {}, ag::Edge::id_type rcid = {});
        Edge &addSPEdge(Vertex &start, Vertex &end, ag::Edge::id_type eid = {}, ag::Edge::id_type rcid = {}); //Deprecated
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeIdType eid = {}, EdgeIdType rcid = {});
        void removeEdge(Edge &edge);
        void isolateAndMarkLockFree(Vertex &vertex);
        template<class I>
        void isolateAndMarkLockFree(I begin, I end);
//        Removes all edges of vertex and vertex.rc() and marks them for deletion, locking both and every neighbour whose
//        edge list changes, so independent vertices can be isolated in parallel.
        void isolateAndMark(Vertex &vertex);
        ag::VertexResolutionResult resolveVertex(Vertex &core, const VertexResolutionPlan &resolution);
        Edge &mergePathToEdge(const GraphPath &path);
        Vertex &mergePath(const GraphPath &path, Vertex::id_type id = 0);
        Vertex &flipVertex(Vertex &core);
        std::unordered_map<EdgeId, EdgeId> ExtendAdjacent(Vertex &core);

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

//    True if the path is a (possibly empty) run of prefix edges followed by a (possibly empty) run of suffix edges.
//    The vertex between the two runs (the peak) then already spells the whole path.
    bool IsPrefixSuffixPath(const GraphPath &path);

    inline std::string GetEdgeNameForSaving(const Edge &edge) {
        VERIFY((edge.getFinish().rc().getInnerId() > 0) == edge.getFinish().rc().isCanonical());
        if(!edge.isCanonical())
            return GetEdgeNameForSaving(edge.rc());
        return edge.getInnerId().str() + "_" + edge.rc().getInnerId().str();
    }

    template<class I>
    void AssemblyGraph::isolateAndMarkLockFree(I begin, I end) {
        std::vector<VertexId> to_mark;
        while(begin != end) {
            to_mark.emplace_back(begin->getId());
            ++begin;
        }
        for(VertexId &vid : to_mark) {
            if(!vid->marked()) {
                isolateAndMarkLockFree(*vid);
            }
        }
    }
}
