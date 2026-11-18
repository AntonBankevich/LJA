#pragma once
#include "assembly_graph/assembly_graph.hpp"
#include "common/oneline_utils.hpp"
#include "common/object_id.hpp"
#include "common/logging.hpp"

namespace spg {

//    Whether the vertex inc enters can be flipped into inc's start: inc is a suffix edge, the vertex's only incoming
//    edge, and every edge leaving the vertex is a prefix edge (the contract of AssemblyGraph::fireReplacePath for a
//    replaced path finish).
    bool CanFlip(const ag::Edge &inc);

//    A Supregraph is an AssemblyGraph whose edges have all been converted via edgeToSupreVertex. A sibling of
//    dbg::SparseDBG (both extend AssemblyGraph directly) rather than a subclass of it: SparseDBG's hashing
//    machinery is DBG-specific and has no business surviving into the Supregraph. Sample segments are maintained
//    by the SampleSegmentTracker every AssemblyGraph owns, and k comes with the base graph.
    class Supregraph : public ag::AssemblyGraph {
//        The sample segment invariant of checkConsistency for one edge.
        bool segmentsFitEdge(ag::Edge &e) const;
    public:
        explicit Supregraph(ag::AssemblyGraph &&base) : ag::AssemblyGraph(std::move(base)) {}
        Supregraph(Supregraph &&other) = default;
        Supregraph &operator=(Supregraph &&other) = default;
        Supregraph(const Supregraph &other) = delete;

//        Takes ownership of an already-built (and, if it came from a SparseDBG, already disableHashing()'d)
//        AssemblyGraph, then converts every non-suffix/non-prefix edge via edgeToSupreVertex -- those are
//        the only edges edgeToSupreVertex accepts (it VERIFYs !isSuffix() && !isPrefix()); any edge that is
//        already one or the other is boundary structure the conversion must leave alone.
        static Supregraph FromDBG(logging::Logger &logger, size_t threads, ag::AssemblyGraph &&graph) {
            Supregraph spg(std::move(graph));
            std::vector<ag::EdgeId> eids = oneline::map(spg.edgesUnique().begin(), spg.edgesUnique().end(), IdTransformer<ag::Edge>());
            for (ag::EdgeId eid : eids) {
                if (!eid->isPrefix() && !eid->isSuffix())
                    spg.edgeToSupreVertex(*eid);
            }
            for (ag::Edge &e : spg.edgesUnique())
                if (e.getVirtualType() == ag::VirtualType::unset)
                    e.setVirtualType(e.isPrefix() && e.isSuffix() ? ag::VirtualType::merge :
                                     e.isSuffix() ? ag::VirtualType::suffix : ag::VirtualType::prefix);
            std::vector<ag::VertexId> vids = oneline::map(spg.vertices().begin(), spg.vertices().end(), IdTransformer<ag::Vertex>());
            for (ag::VertexId vid : vids) {
                if (vid->inDeg() == 1 && vid->outDeg() != 1 && CanFlip(vid->incFront()))
                    spg.flipVertex(vid->incFront());
            }
            spg.removeMarked();
            return spg;
        }

//        AssemblyGraph's checks plus: every edge u -> v is virtually a suffix or prefix edge (or both), i.e. MES(u) is a
//        prefix of MES(v) or MES(v) a suffix of MES(u); structurally, u has out-degree 1 or v has in-degree 1. It holds
//        once edgeToSupreVertex has run on every edge of the initial graph, and no operation may break it. Also: an edge
//        from a vertex to itself only when that loop is the vertex's whole connected component (edgeToSupreVertex turns
//        the DBG's loops into supre-vertices). Sample segments at edges: for every edge u -> v with overlap t, either the
//        last segment of u ends at the end of u and the first segment of v ends at t, or the first segment of v ends at
//        t + 1. A segment-less end (e.g. a k-mer) contributes nothing: the (k+1)-mer is the other end's segment. Loops
//        (u == v) are exempt.
        bool checkConsistency(std::ostream &os = std::cerr) override;
    };

//    What is merged of an unbranching path (from junction to junction): the path without its start if that has several
//    outgoing edges and without its finish if that has several incoming edges (no edges if nothing is left). An
//    isolated cycle is merged as is.
    ag::GraphPath MergeablePath(const ag::GraphPath &path);
//    Merges MergeablePath(path). Returns the merged vertex.
    ag::Vertex &MergePathOrLoop(const ag::GraphPath &path, Supregraph &graph, ag::Vertex::id_type id = 0);
//    Merges the unbranching path (or cycle) through the non-junction vertex with MergePathOrLoop. Returns the merged
//    vertex. Single-threaded.
    ag::Vertex &MergeUnbranching(logging::Logger &logger, Supregraph &graph, ag::Vertex &vertex);
    void MergeAllSPG(logging::Logger &logger, size_t threads, Supregraph &graph);
//    Only removes outer vertices. It is iterative and removes unsupported vertices outer layer by outer layer.
    size_t SimpleRemoveUncovered(logging::Logger &logger, size_t threads, Supregraph &spg);
}
