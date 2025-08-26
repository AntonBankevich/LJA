#pragma once

#include "assembly_graph/graph_listeners.hpp"
#include "assembly_graph/aligned_read.hpp"
#include "aligned_read_listeners.hpp"
#include <common/logging.hpp>
#include <common/omp_utils.hpp>

namespace ag {

    class AlignedReadStorageMaintenance;

//    TODO: create a class VirtualReadStorage that combines several storages together
//  TODO: move starts to Maintenance and add AddRead method, which can only be invoked before Maintenance is enabled
    class AlignedReadStorage : public AlignedReadStorageFire {
        friend class AlignedReadStorageMaintenance;
    private:
        std::vector <AlignedRead> reads;
        std::unordered_map <ConstEdgeId, std::vector<AlignedReadDirection>> starts;
        mutable omp_lock_t writelock = {};
        ag::AlignedReadStorageMaintenance * maintenance = nullptr;

        void updateStart(AlignedReadDirection dir);
    public:
        void lock() const { omp_set_lock(&writelock); }
        void unlock() const { omp_unset_lock(&writelock); }

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

        const std::vector<AlignedReadDirection> &getOutgoingReadsLockFree(Edge &edge) const {
            return starts.at(edge.getId());
        }
//        TODO: Make this free of global lock
        const std::vector<AlignedReadDirection> &getOutgoingReads(Edge &edge) const {
            lock();
            const std::vector<AlignedReadDirection> & res = getOutgoingReadsLockFree(edge);
            unlock();
            return res;
        }

        std::vector<AlignedReadDirection> &getOutgoingReadsLockFree(Edge &edge) {
            return starts.at(edge.getId());
        }

        std::vector<AlignedReadDirection> &getOutgoingReads(Edge &edge) {
            lock();
            std::vector<AlignedReadDirection> & res = getOutgoingReadsLockFree(edge);
            unlock();
            return res;
        }

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
        size_t startCnt(const Edge &edge) const {return starts.at(edge.getId()).size();}

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

        void fireAddEdge(Edge &edge) override {
            if(!edge.isPrefix()) {
                storage->lock();
                storage->starts[edge.getId()] = {};
                storage->unlock();
            }
        }
        void fireDeleteEdge(Edge &edge) override {
            if(!edge.isPrefix()) {
                storage->lock();
                storage->starts.erase(edge.getId());
                storage->unlock();
            }
        }
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
//        TODO: implement properly
        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) override;;
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override;
        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) override;
        void fireAddSupreVertex(Vertex &v, Edge &e) override;

        void fireAddRead(const AlignedRead &read) override;
        void fireRerouteRead(AlignedRead &read) override;
        void fireInvalidateRead(AlignedRead &read) override;
        bool fireCheckConsistency() override;
    };

//    TODO: get rid of this function as soon as Andrey's code is out.
    void SaveReads(const std::experimental::filesystem::path &fname, const AlignedReadStorage &storage);
}
