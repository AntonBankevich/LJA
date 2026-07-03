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
    private:
        typedef std::vector<AlignedReadDirection> DirectionRecord;
        // Both maps store unique_ptrs by value and move them around.
        // This keeps a record's address stable across cuckoohash_map rehashes, so a reference
        // returned by the accessors below stays valid even while other keys are inserted/erased,
        // and provides concurrent insert/erase/find with local locks (see SuffixTracker::edge_data).
        std::vector <AlignedRead> reads;
        libcuckoo::cuckoohash_map <ConstEdgeId, std::unique_ptr<DirectionRecord>> starts;
        libcuckoo::cuckoohash_map <ConstVertexId, std::unique_ptr<DirectionRecord>> reads_inside_vertices;
        ag::AlignedReadStorageMaintenance * maintenance = nullptr;

        void updateStart(AlignedReadDirection dir);

//        Canonical getters: every read of a stored DirectionRecord (from inside this class,
//        from AlignedReadStorageMaintenance, or from outside) goes through these two pairs of
//        functions rather than calling find_fn on `starts`/`reads_inside_vertices` directly.
//        Structural map operations (insert/erase/contains/lock_table) are used directly at their
//        call sites, since AlignedReadStorageMaintenance is a friend of this class.
        DirectionRecord &getOutgoingReadsRecord(ConstEdgeId eid) {
            DirectionRecord *res = nullptr;
            starts.find_fn(eid, [&res](const std::unique_ptr<DirectionRecord> &ptr) { res = ptr.get(); });
            VERIFY(res != nullptr);
            return *res;
        }
        const DirectionRecord &getOutgoingReadsRecord(ConstEdgeId eid) const {
            const DirectionRecord *res = nullptr;
            starts.find_fn(eid, [&res](const std::unique_ptr<DirectionRecord> &ptr) { res = ptr.get(); });
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

        const std::vector<AlignedReadDirection> &getOutgoingReadsLockFree(const Edge &edge) const;
        const std::vector<AlignedReadDirection> &getOutgoingReads(const Edge &edge) const;
        std::vector<AlignedReadDirection> &getOutgoingReadsLockFree(const Edge &edge);
        std::vector<AlignedReadDirection> &getOutgoingReads(const Edge &edge);
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
        size_t startCnt(const Edge &edge) const {return getOutgoingReadsRecord(edge.getId()).size();}

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
    public:
        AlignedReadStorageMaintenance(AssemblyGraph &graph, AlignedReadStorage &storage);
        AlignedReadStorageMaintenance(AlignedReadStorageMaintenance &&other) = default;

        void fireAddVertex(Vertex &vertex) override;
        void fireDeleteVertex(Vertex &vertex) override;
        void fireAddEdge(Edge &edge) override;
        void fireDeleteEdge(Edge &edge) override;

        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
//        TODO: implement properly
        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) override;
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
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
