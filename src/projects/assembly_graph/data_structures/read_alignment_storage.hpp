#pragma once

#include "assembly_graph/graph_listeners.hpp"
#include "assembly_graph/aligned_read.hpp"
#include "aligned_read_listeners.hpp"
#include <common/logging.hpp>
#include <common/omp_utils.hpp>
#include "libcuckoo/cuckoohash_map.hh"

namespace ag {

    class AlignedReadStorageMaintenance;

//    TODO: create a class VirtualReadStorage that combines several storages together
//  TODO: move starts to Maintenance and add AddRead method, which can only be invoked before Maintenance is enabled
    class AlignedReadStorage : public AlignedReadStorageFire {
        friend class AlignedReadStorageMaintenance;
        void Normalize(GraphPath &path) {
            size_t right_cut = path.rightCut();
            path.setCutRight(0);
            while (path.getFinish().outDeg() == 1 && !path.getFinish().frontVertex().coverage_info.empty() &&
                path.getFinish().frontVertex().size() - path.getFinish().frontVertex().coverage_info.front().finish() <=
                right_cut + path.getFinish().front().truncSize() ) {
                right_cut += path.getFinish().front().truncSize();
                path += path.getFinish().front();
            }
            path.setCutRight(right_cut);
            size_t left_cut = path.leftCut();
            path.setCutLeft(0);
            while (path.getStart().inDeg() == 1 && !path.getStart().incFrontVertex().coverage_info.empty() &&
                path.getStart().incFrontVertex().coverage_info.back().start() <= left_cut + path.getStart().incFront().rc().truncSize() ) {
                path.push_front(path.getStart().incFront());
                left_cut += path.getStart().incFront().rc().truncSize();
            }
        }
    private:
        typedef std::vector<AlignedReadDirection> DirectionRecord;
        // Both maps store unique_ptrs by value and move them around.
        // This keeps a record's address stable across cuckoohash_map rehashes, so a reference
        // returned by the accessors below stays valid even while other keys are inserted/erased,
        // and provides concurrent insert/erase/find with local locks (see SuffixTracker::edge_data).
        std::vector <AlignedRead> reads;
        libcuckoo::cuckoohash_map <ConstVertexId, std::unique_ptr<DirectionRecord>> starts;
        libcuckoo::cuckoohash_map <ConstVertexId, std::unique_ptr<DirectionRecord>> reads_inside_vertices;
        ag::AlignedReadStorageMaintenance * maintenance = nullptr;

        void updateStart(AlignedReadDirection dir);

//        Canonical getters: every read of a stored DirectionRecord (from inside this class,
//        from AlignedReadStorageMaintenance, or from outside) goes through these two pairs of
//        functions rather than calling find_fn on `starts`/`reads_inside_vertices` directly.
//        Structural map operations (insert/erase/contains/lock_table) are used directly at their
//        call sites, since AlignedReadStorageMaintenance is a friend of this class.
        DirectionRecord &getOutgoingReadsRecord(ConstVertexId vid) {
            DirectionRecord *res = nullptr;
            starts.find_fn(vid, [&res](const std::unique_ptr<DirectionRecord> &ptr) { res = ptr.get(); });
            VERIFY(res != nullptr);
            return *res;
        }
        const DirectionRecord &getOutgoingReadsRecord(ConstVertexId vid) const {
            const DirectionRecord *res = nullptr;
            starts.find_fn(vid, [&res](const std::unique_ptr<DirectionRecord> &ptr) { res = ptr.get(); });
            VERIFY(res != nullptr);
            return *res;
        }
        DirectionRecord &getSubstringReadsRecord(ConstVertexId vid) {
            DirectionRecord *res = nullptr;
            reads_inside_vertices.find_fn(vid, [&res](const std::unique_ptr<DirectionRecord> &ptr) { res = ptr.get(); });
            VERIFY(res != nullptr);
            return *res;
        }
        const DirectionRecord &getSubstringReadsRecord(ConstVertexId vid) const {
            const DirectionRecord *res = nullptr;
            reads_inside_vertices.find_fn(vid, [&res](const std::unique_ptr<DirectionRecord> &ptr) { res = ptr.get(); });
            VERIFY(res != nullptr);
            return *res;
        }
    public:
        bool checkConsistency();

        AlignedReadStorage() = default;
        AlignedReadStorage(AlignedReadStorage &&other)  noexcept;
        AlignedReadStorage(const AlignedReadStorage &other) = delete;
        AlignedReadStorage &operator=(AlignedReadStorage &&other) noexcept;
        AlignedReadStorage &operator=(const AlignedReadStorage &other) = delete;

        AlignedReadStorage(AssemblyGraph &graph, std::vector<AlignedRead> reads);
        AlignedReadStorage(logging::Logger &logger, size_t threads, AssemblyGraph &graph,
                           std::vector<AlignedRead> reads);
        virtual ~AlignedReadStorage();

//        Reads whose canonical path starts at this vertex (non-empty paths: the read leaves the vertex on
//        one of its outgoing edges). Reads starting on a specific edge are this set filtered by
//        dir.frontEdge() == edge.
        const std::vector<AlignedReadDirection> &getOutgoingReadsLockFree(const Vertex &vertex) const;
        const std::vector<AlignedReadDirection> &getOutgoingReads(const Vertex &vertex) const;
        std::vector<AlignedReadDirection> &getOutgoingReadsLockFree(const Vertex &vertex);
        std::vector<AlignedReadDirection> &getOutgoingReads(const Vertex &vertex);
        const std::vector<AlignedReadDirection> &getSubstringReadsLockFree(VertexId vertex) const;
        const std::vector<AlignedReadDirection> &getSubstringReads(VertexId vertex) const;
        std::vector<AlignedReadDirection> &getSubstringReadsLockFree(VertexId vertex);
        std::vector<AlignedReadDirection> &getSubstringReads(VertexId vertex);

        //        typename std::vector<AlignedRead>::iterator begin() { return reads.begin(); }
//        typename std::vector<AlignedRead>::iterator end() { return reads.end(); }
        typename std::vector<AlignedRead>::const_iterator begin() const { return reads.begin(); }
        typename std::vector<AlignedRead>::const_iterator end() const { return reads.end(); }
        typename std::vector<AlignedRead>::iterator begin() { return reads.begin(); }
        typename std::vector<AlignedRead>::iterator end() { return reads.end(); }
//        AlignedRead &operator[](size_t ind) { return reads[ind]; }
        const AlignedRead &operator[](size_t ind) const { return reads[ind]; }
        AlignedRead &operator[](size_t ind) { return reads[ind]; }

        size_t size() const { return reads.size(); }

        void delayedInvalidateRead(AlignedRead &read, const std::string &message);
        void rerouteRead(AlignedRead &alignedRead, GraphPath corrected, const string &message);
        bool apply(AlignedRead &alignedRead);

        void applyCorrections(logging::Logger &logger, size_t threads);

        void printReadPaths(logging::Logger &logger, const std::experimental::filesystem::path &aln_path,
                            const std::experimental::filesystem::path &gfa_path,
                            const std::experimental::filesystem::path &rp_path,
                            size_t k) const;
        void printReadFasta(logging::Logger &logger, const std::experimental::filesystem::path &path) const;
        void printFullAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const;
        void printSequences(const std::experimental::filesystem::path &path) const;


        void Save(const std::experimental::filesystem::path &path);
        void Save(std::ostream &os) const;
        static AlignedReadStorage Load(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path,
                                       AssemblyGraph &graph);
        static AlignedReadStorage Load(logging::Logger &logger, size_t threads, std::istream &is, AssemblyGraph &graph);

        static std::vector<AlignedRead> LoadReadAlignments(std::istream &is, IdIndex<Vertex> &index);
    };

//    TODO: move to dbg

    class AlignedReadStorageMaintenance : public AlignedReadStorageListener, public ResolutionListener {
        friend class AlignedReadStorage;
    private:
        AlignedReadStorage *storage;
//        Moves every read recorded as a pure substring of `from` (reads_inside_vertices[from], i.e. reads
//        that never reach an edge) onto `to`, re-expressing cut_left/cut_right for from's position-0
//        mapping onto to's position `shift` -- the exact read-side counterpart of how
//        SampleSegmentTracker moves a vertex's coverage_info (shift()/insertFront()/+=) in
//        fireResolveVertex's three unambiguous cases, where core is destroyed but has exactly one real
//        neighbor to transfer onto. Must be called with the same destination and shift
//        SampleSegmentTracker uses for `from`, so a read's recorded position never drifts from where its
//        sample segment actually lives.
        void migrateSubstringReads(Vertex &from, Vertex &to, size_t shift);
    public:
        AlignedReadStorageMaintenance(AssemblyGraph &graph, AlignedReadStorage &storage);
        AlignedReadStorageMaintenance(AlignedReadStorageMaintenance &&other) = default;

        void fireAddVertex(Vertex &vertex) override;
        void fireDeleteVertex(Vertex &vertex) override;
//        `starts` and `reads_inside_vertices` are keyed by vertex now (see fireAddVertex/fireDeleteVertex),
//        so edge addition/removal needs no bookkeeping here; the base-class no-op handlers are used.

        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
//        TODO: implement properly
        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) override;
//        Relies on merged edge codes being exact concatenations (see EdgeCodeListener): a read whose path
//        merely passes through the merged run needs no update at all, so this only patches reads indexed
//        as starting on an interior edge of the run (storage->getOutgoingReads(edge)), pushing the prefix
//        before their original start onto their path so it now starts at new_edge/path.getStart().
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
//        Hidden, unenforced caller contract on AssemblyGraph::splitEdge: split_positions must never fall
//        inside a read's span. Every read touching `edge` must, after the split, end within the first
//        piece, start within the last piece, or be fully contained in one piece — never span across a new
//        split-point vertex. This handler only repoints reads indexed as *starting* at `edge`
//        (storage->getOutgoingReads(edge)); a pass-through read that violated the contract would be left
//        with stale path codes and crash later in Vertex::getOutgoingByIterator, not here.
        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override;
        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) override;
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;

        void fireAddRead(const AlignedRead &read) override;
        void fireRerouteRead(AlignedRead &read) override;
        void fireInvalidateRead(AlignedRead &read) override;
        bool fireCheckConsistency() override;
    };

//    TODO: get rid of this function as soon as Andrey's code is out.
    void SaveReads(const std::experimental::filesystem::path &fname, const AlignedReadStorage &storage);
}
