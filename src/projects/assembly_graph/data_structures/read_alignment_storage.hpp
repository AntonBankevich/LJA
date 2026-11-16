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
//        End of the leftmost sample segment of a sample record: a kpomer chunk stands for each of its (k+1)-mers.
        size_t firstSegmentFinish(const CoverageSamples::SampleView &sample) const {
            if (sample.sample.type == CoverageSamples::SampleType::kpomer) {
                VERIFY(k > 0);
                return sample.start() + k + 1;
            }
            return sample.finish();
        }
//        Extends the path forward through unbranching steps to every vertex whose first sample segment lies inside the
//        read, passing through vertices without segments; trailing vertices without segments are dropped again. A DBG
//        path (ending in a k-mer) is already normal.
        void ExtendForward(GraphPath &path) const {
            if (path.getFinish().size() == k)
                return;
            size_t cut_right = path.rightCut();
            path.setCutRight(0);
            while (!path.empty() && path.backEdge().isSuffix())
                path.pop_back();
            while (path.getFinish().outDeg() == 1) {
                Edge &next = path.getFinish().front();
                if (next.getFinish().coverage_info.empty() || firstSegmentFinish(next.getFinish().coverage_info.front()) + next.rc().truncSize() + cut_right <= path.getFinish().size()) {
                    cut_right += next.truncSize();
                    path += next;
                } else {
                    break;
                }
            }
            while (!path.empty() && path.getFinish().coverage_info.empty()) {
                cut_right -= path.backEdge().truncSize();
                path.pop_back();
            }
            path.setCutRight(cut_right);
        }

//        Invalid paths (unaligned reads) and legacy paths (their vertex is gone) are returned as they are.
        GraphPath Normalize(GraphPath path) const {
            if (!path.valid() || path.isLegacy())
                return path;
            ExtendForward(path);
            path.swapDirection();
            ExtendForward(path);
            path.swapDirection();
            return path;
        }
    private:
        typedef std::vector<AlignedReadDirection> DirectionRecord;
        // Both maps store unique_ptrs by value and move them around.
        // This keeps a record's address stable across cuckoohash_map rehashes, so a reference
        // returned by the accessors below stays valid even while other keys are inserted/erased,
        // and provides concurrent insert/erase/find with local locks (see SuffixTracker::edge_data).
        std::vector <AlignedRead> reads;
//        k of the DBG the reads were aligned to. Kpomer sample chunks stand for all their (k+1)-mers, so read path
//        normalization and its checks need k to know where the individual sample segments are.
        size_t k = 0;
        libcuckoo::cuckoohash_map <ConstVertexId, std::unique_ptr<DirectionRecord>> starts;
        libcuckoo::cuckoohash_map <ConstVertexId, std::unique_ptr<DirectionRecord>> reads_inside_vertices;
//        Subreads of deleted vertices, keyed by the deleted vertex's id (legacyId()). Their paths are LegacyPaths.
        libcuckoo::cuckoohash_map <ConstVertexId, std::unique_ptr<DirectionRecord>> legacy_reads;
        ag::AlignedReadStorageMaintenance * maintenance = nullptr;
//        See startReadCorrection.
        bool read_correction = false;

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
//        Checks the read records without needing the graph: every valid read is filed exactly once per direction in the
//        right map (starts / reads_inside_vertices / legacy_reads) under the right vertex, with no stale entries; every
//        non-legacy path starts and ends at live vertices and can be followed edge by edge with its codes; cuts leave
//        some sequence; a k-mer path end is needed by the read, any other end has a sample segment inside the read; the
//        path equals its normal form. Reports every violation to os.
        bool checkRecords(std::ostream &os);
    public:
//        checkRecords plus every AlignedReadStorageListener's fireCheckConsistency.
        bool checkConsistency(std::ostream &os = std::cerr);

        AlignedReadStorage() = default;
        AlignedReadStorage(AlignedReadStorage &&other)  noexcept;
        AlignedReadStorage(const AlignedReadStorage &other) = delete;
        AlignedReadStorage &operator=(AlignedReadStorage &&other) noexcept;
        AlignedReadStorage &operator=(const AlignedReadStorage &other) = delete;

//        k is taken from graph.
        AlignedReadStorage(AssemblyGraph &graph, std::vector<AlignedRead> reads);
        AlignedReadStorage(logging::Logger &logger, size_t threads, AssemblyGraph &graph,
                           std::vector<AlignedRead> reads);
        virtual ~AlignedReadStorage();
        size_t getK() const {return k;}

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

//        The storage is in one of two modes. In graph mode (the default) the graph may change and reads may not. In read
//        correction mode, entered by startReadCorrection, reads may be rerouted or invalidated and the graph may not
//        change: graph listeners only maintain current paths, never pending corrections. applyCorrections applies all
//        pending corrections and returns to graph mode.
        void startReadCorrection() {VERIFY(!read_correction); read_correction = true;}
        bool inReadCorrection() const {return read_correction;}
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
//        TODO: the graph + read set save/load mechanism does not preserve sample segments (Vertex::coverage_info).
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
//        Moves the subreads of the vertex flipped along edge (its finish if edge is a suffix edge, its start if it is
//        a prefix edge) to the other end of edge, which contains it: the path is extended by edge and then loses it on
//        the other side, so its codes are only pushed and popped (their deck positions, and so the reads' suffix ends,
//        stay valid).
        void migrateSubstringReads(Edge &edge);
//        Moves the reads of path's inner vertices to target, which spells the whole path (path position 0 is target's
//        position 0). The new edge from_start leads from the path's start to target and to_finish from target to the
//        path's finish (nullptr where target is that end); together they have exactly the path's codes. A read is
//        extended edge by edge along the path to its start (and to its finish if it ends inside the path), and the new
//        edges are popped: reads reaching the finish start at target, the others, like the inner vertices' subreads,
//        become subreads of target. Codes are only pushed and popped, so their deck positions (and the reads' suffix
//        ends) stay valid.
        void absorbInnerVertices(const RAGraphPath &path, Vertex &target, Edge *from_start, Edge *to_finish);
//        Files dir, whose path has just become empty, as a subread, together with its other direction (moved out of
//        the starts of dir's start.rc()).
        void fileAsSubread(const AlignedReadDirection &dir);
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
//        Reads starting at, and subreads of, the removed vertices move to the vertex that contains them.
        void fireMergePrefixPath(const RAGraphPath &path, Edge &new_edge) override;
        void fireFlipVertex(Edge &edge, const std::unordered_map<EdgeId, EdgeId> &extensions) override;

        void fireAddRead(const AlignedRead &read) override;
        void fireRerouteRead(AlignedRead &read) override;
        void fireInvalidateRead(AlignedRead &read) override;
//        Every live vertex (and no deleted one) has a starts and a reads_inside_vertices record, plus
//        storage->checkRecords.
        bool checkConsistency(std::ostream &os) override;
    };

//    TODO: get rid of this function as soon as Andrey's code is out.
    void SaveReads(const std::experimental::filesystem::path &fname, const AlignedReadStorage &storage);
}
