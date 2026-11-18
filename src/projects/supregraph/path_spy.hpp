#pragma once

#include <filesystem>

#include "unique_vertex_storage.hpp"
#include "assembly_graph/visualization.hpp"
#include "assembly_graph/data_structures/read_alignment_storage.hpp"
#include "assembly_graph/data_structures/component.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/graph_listeners.hpp"
#include "dbg/multi_graph.hpp"

namespace spg {
    class OldVertexTracker : public SPGListener {
    public:
        struct DirectEmbedding {
            VertexId inner_vertex;
            VertexId outer_vertex;
            size_t from;
            size_t to;
            DirectEmbedding(Vertex & inner_vertex, Vertex & outer_vertex, size_t from, size_t to) :
                    inner_vertex(inner_vertex.getId()), outer_vertex(outer_vertex.getId()), from(from), to(to) {}
            DirectEmbedding(VertexId inner_vertex, VertexId outer_vertex, size_t from, size_t to) :
                    inner_vertex(inner_vertex), outer_vertex(outer_vertex), from(from), to(to) {}
            bool operator<(const DirectEmbedding &other) const;
            bool operator==(const DirectEmbedding &other) const;
        };
    private:
        typedef std::vector<DirectEmbedding> EmbeddingRecord;
//        Merges and deletions fire in parallel, so the maps are concurrent. Records are unique_ptrs so that a reference
//        returned by getRecord stays valid across rehashes (see AlignedReadStorage).
        typedef libcuckoo::cuckoohash_map<ag::ConstVertexId, std::unique_ptr<EmbeddingRecord>> EmbeddingMap;
//        Every vertex ever added. Sequences are stored only in debug mode (empty otherwise) and only serve the
//        sequence checks of checkConsistency and getAllSubvertices.
        libcuckoo::cuckoohash_map<ag::ConstVertexId, Sequence> vertex_seq;
        //A vertex has a recorded self-embedding iff it was already removed from the graph
        EmbeddingMap direct_supervertices;
        EmbeddingMap direct_subvertices;
        bool debug;

        static void append(EmbeddingMap &map, ag::ConstVertexId key, const DirectEmbedding &embedding);
        static const EmbeddingRecord &getRecord(const EmbeddingMap &map, ag::ConstVertexId key);
        void addMultiEmbedding(Vertex & subvertex, Vertex & supervertex, size_t left, size_t right);
        size_t getVertexSize(VertexId vid) const;

    public:
        OldVertexTracker(ag::ResolutionFire &fire, bool debug);

//        Whether vid (live or deleted) was ever seen by the tracker. Vertices deleted before the tracker was attached
//        were not.
        bool checkExists(VertexId vid) const {return vertex_seq.contains(ag::ConstVertexId(vid));}


        void fireAddEdge(Edge &e) override;
        void fireAddVertex(Vertex &v) override;
        void fireDeleteVertex(Vertex &v) override;
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override {}
        void fireReplacePath(const ag::PathReplacement &replacement) override;
        void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) override {VERIFY(false);}
        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;

        std::vector<DirectEmbedding> getAllEmbeddings(VertexId v) const;
        std::vector<DirectEmbedding> getAllSubvertices(VertexId v) const;

//        Sequence parts (marked "debug") are checked only in debug mode.
//        1. Every vertex ever added has records; (debug) a live vertex's recorded sequence is its sequence.
//        2. Every embedding (inner, outer, from, to) has from <= to; (debug) seq(inner) == seq(outer)[from, to).
//        3. Every embedding is recorded both as a supervertex of inner and a subvertex of outer.
//        4. A vertex carries the dead marker (a self-embedding) iff it is not live.
//        5. Every live prefix (suffix) edge's start (finish) is recorded inside its finish (start).
        bool checkConsistency(std::ostream &os) override;
    };

//    "Which currently-alive vertices is this read a substring of?" A read whose MES extension fit inside
//    a single vertex is equivalent to a vertex V = the vertex its path lies within (V may already be
//    deleted, in which case the read carries a legacy path whose start is V). This walks
//    OldVertexTracker's embedding records up from V, past vertices that have since been removed, to the
//    minimal living vertices that still contain V's sequence, and returns them (there may be more than
//    one). Only meaningful for reads equivalent to a vertex (substring/legacy reads); a read spanning a
//    real path corresponds to a compact path instead and is out of scope here.
    std::vector<ag::VertexId> readSubstringVertices(const ag::AlignedRead &read, const OldVertexTracker &tracker);

//    Writes read alignments in the format polishing reads (that of PrintAlignments in polishing/perfect_alignment.hpp:
//    "read_id read_from read_to contig_id contig_from contig_to" per line, half-open, reads in storage order) straight
//    from the stored read paths instead of realigning reads. Contigs are the vertices: a canonical vertex is named by
//    its id, a non-canonical one by its (negative) id with coordinates on its own sequence, i.e. on the reverse
//    complement of the canonical contig. Read coordinates are on the read's path sequence (the corrected compressed
//    read). A read spanning several vertices gets one record per vertex of its path. A read lying inside a single
//    vertex V (possibly deleted) gets one record per living vertex containing V, found through tracker once per V.
//    Records that reach neither a read end nor a contig end on some side go to the partial file.
//    Returns {good alignments file, partial alignments file}.
    std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
    PrintStorageAlignments(logging::Logger &logger, size_t threads, const ag::AlignedReadStorage &storage,
                           const OldVertexTracker &tracker, const std::experimental::filesystem::path &dir);
}



//class PathSpy : ag::ResolutionListener {
//private:
//    ag::AlignedReadStorage *storage;
//    std::unordered_map<std::string, std::pair<size_t, std::experimental::filesystem::path>> read_dirs;
//public:
//    PathSpy(ag::AlignedReadStorage &storage, const std::experimental::filesystem::path &dir) : storage(&storage){
//        recreate_dir(dir);
//        for(ag::AlignedRead &read : *storage) {
//            std::experimental::filesystem::path new_dir = dir/read.getId();
//            recreate_dir(new_dir);
//            read_dirs[read.getId()] = {0, new_dir};
//        }
//    }
//
//    DrawRead(const ag::AlignedRead &read) {
//        auto &p = read_dirs.at(read.getId());
//        ++p.first;
//        const ag::GraphPath &path = read.getPath();
//        ag::Component component = ag::Component::neighbourhood(graph, path, 10000, path.calculateSize() * 2 + 20);
//        std::experimental::filesystem::path fpath = p.second / (stoi(p.first)+".dot");
//        ag::Printer printer;
//        printer.printDot(fpath, component);
////        Need labeler, colorer and calls to this function
//    }
//
//    virtual void fireEdgeToSupreVertex(Vertex &v, Edge &e) {
//    }
//
//    virtual void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
//
//    }
//    virtual void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {}
//    virtual void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {}
//    virtual void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
//                                     const AlignmentForm &left_al, const AlignmentForm &right_al) {}
//    virtual void fireSplitEdge(Edge &edge, const RAGraphPath &split) {}
//
//    virtual void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {}
//
//    virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {};
//
//};