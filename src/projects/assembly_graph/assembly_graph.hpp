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

//        The core of mergePath, flipVertex and extendVertex: replaces old_path P = p0 -> ... -> pn by new_path Q, which
//        spells the same sequence and whose vertices and edges the caller has already added. dominant_index is the index
//        in Q of the vertex that receives the information of the removed vertices. The edges around replaced ends are
//        added here: if p0 is replaced by Q's start q0 (p0 a prefix of q0), every y -> p0 becomes y -> q0; if pn is
//        replaced by Q's finish qe (pn a suffix of qe), every pn -> z becomes qe -> z; the loop pn -> p0 becomes qe -> q0
//        once. Each new edge spells its old one extended to the end of q0 (the start of qe) and inherits its virtual
//        type. Keeps the new edges non-corporeal while the listeners react (fireReplacePath) and the vertices of P that
//        are not in Q are isolated and marked.
        void replacePathLockFree(const RAGraphPath &old_path, const RAGraphPath &new_path, size_t dominant_index);

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
//        Merges a path p0 -> ... -> pn (at least one edge) into a single vertex D spelling it: the path's vertex with
//        exactly that sequence if there is one, else a new vertex. Every other vertex of the path is removed (isolated
//        and marked). A replaced p0 must have out-degree 1, and its incoming edges y -> p0 become y -> D; a replaced pn
//        must have in-degree 1, and its outgoing edges pn -> z become D -> z; both inherit their types. Inner vertices
//        must have no other edges. A cycle (p0 = pn) is merged without its last edge, which becomes the loop D -> D, as
//        does an edge pn -> p0 when both ends are replaced. Returns D.
        Vertex &mergePathLockFree(const GraphPath &path, Vertex::id_type id = 0);
        GraphPath splitEdge(Edge &edge, const std::vector<EdgePosition> &split_positions);//        TODO: this should run in parallel
        //SPG specific
        Vertex &edgeToSupreVertex(Edge &edge, Vertex::id_type id = 0);
        Vertex &edgeToSupreVertexLockFree(Edge &edge, Vertex::id_type id = 0);
//        Requires suffix edges into core, prefix edges out of it, and sample segments on every neighbour of core. A pair
//        with a single partner on either side becomes a direct edge; only pairs with several partners on both sides get
//        a new vertex (which then gets a resolution sample).
        ag::VertexResolutionResult resolveCoreVertex(Vertex &core, const VertexResolutionPlan &resolution); //TODO:: this should run in parallel
        Vertex &mergeLoop(const GraphPath &path); //TOOD: implement
//        Flips core = inc.getFinish() into u = inc.getStart() along the suffix edge inc, the only edge into core:
//        every core -> w becomes u -> w and core is isolated and marked (the merge of inc). core may have any
//        out-degree. Returns u.
        Vertex &flipVertexLockFree(Edge &inc);
//        suffix_edge = x -> w is virtually suffix (x's only outgoing edge) and neither a prefix nor a suffix edge. x is
//        replaced by a new vertex v spelling suffix_edge (x extended to the end of w): every y -> x becomes y -> v and
//        x -> w the suffix edge v -> w; x is isolated and marked. Same result as edgeToSupreVertex followed by a flip of
//        x.rc() into the supre-vertex, without the intermediate edges. Requires no edge between x and x.rc() or w.
//        A self-rc x is extended both ways, so that v is self-rc too: its only incoming edge is w.rc() -> x, the rc of
//        suffix_edge, and v spells w.rc() -> x -> w (x extended to the start of w.rc() and to the end of w); the path
//        becomes the prefix edge w.rc() -> v followed by the suffix edge v -> w, its rc. Requires w != w.rc().
//        Returns v.
        Vertex &extendVertexLockFree(Edge &suffix_edge, Vertex::id_type id = 0);

//-------------------Derived graph operations. These call core graph operations.-------------------------
        Vertex &addVertex(const Sequence &seq, typename Vertex::id_type id = Vertex::id_type()) {return addVertex(seq, {}, id);}
        Vertex &addVertex(const Vertex &other_graph_vertex);
        Edge &addEdgeLockFree(Vertex &start, Vertex &end, const Sequence &full_sequence, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &full_seq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addEdge(Vertex &start, Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data, EdgeIdType eid = {}, EdgeIdType rcid = {});
        Edge &addSPEdgeLockFree(Vertex &start, Vertex &end, ag::Edge::id_type eid = {}, ag::Edge::id_type rcid = {});
//        The same edge with the given virtual type instead of the one following from its real prefix/suffix nature.
        Edge &addSPEdgeLockFree(Vertex &start, Vertex &end, VirtualType type);
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
//        Locks the path's ends and merges it into a single edge. A path of a single edge (isSingleton) is already merged:
//        that edge is returned unchanged. The check is done under the lock, since other threads merging paths from the
//        same ends change their edge lists.
        Edge &mergePathToEdge(const GraphPath &path);
        Vertex &mergePath(const GraphPath &path, Vertex::id_type id = 0);
        Vertex &flipVertex(Edge &inc);
//        Makes every edge into core (on both strands) a suffix edge: its start is replaced by itself extended to the end
//        of core (extendVertex; a self-rc start is extended both ways). Returns old edge id -> the edge now
//        adjacent to core. Replaced vertices are only marked.
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
