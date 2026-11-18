#pragma once
#include "assembly_graph/assembly_graph.hpp"
#include "read_alignment_storage.hpp"
#include "read_logger.hpp"
#include "fstream"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "libcuckoo/cuckoohash_map.hh"
#include <experimental/filesystem>

namespace ag {

    class SuffixTracker;
//    Contract: sequences representing the same path are considered equivalent. While graph changes,
//    merging operations change the equivalency.
//    Invariants:
//    Once equivalent->always equivalent.
//    No negative values are stored in the table.
//    Queries to equivalent sequences have non-decreasing length
//    If a sequence has non-zero multiplicity, it represents a valid path in the graph.
//    However the last edge may have incomplete code in the sequence.
//    No such guarantee is given for sequences of multiplicity 0.
//TODO: iterate only through sequences of non-zero multiplicity to encapsulate this invariant
//TODO: store pairs of NuclDeck iterators instead of seqences and make sure they stay alive. This uses a bit less memory and small object generation, less synchronization
    struct SuffixRecord {
        friend SuffixTracker;
    private:
        typedef std::vector<std::pair<Sequence, int>> Storage;
        typedef Storage::const_iterator const_iterator;
        typedef Storage::iterator iterator;
        EdgeId eid;
        Storage paths;
        size_t zero_cnt = 0;
        size_t max_suffix_len;

        void lock() const { eid->getStart().lock(); }
        void unlock() const { eid->getStart().unlock(); }
        void updateZero(size_t old_val, size_t new_val);
//        Adds diff to the count of seq; a removal (diff < 0) must find seq stored exactly, with enough count.
        void lockFreeChangePathCnt(const Sequence &seq, int diff);
        void changePathCnt(const Sequence &seq, int diff);
        void addPath(const Sequence &seq, int diff = 1);
        void directAddPath(const Sequence &seq, size_t cnt);
        void clear();

        void removeZero();
        size_t countStartsWith(const Sequence &seq) const;
        void resetCodes(Vertex &start);
    public:
        explicit SuffixRecord(Edge &edge, size_t max_suffix_len) : eid(edge.getId()), max_suffix_len(max_suffix_len) {}
        SuffixRecord(const SuffixRecord &) = delete;
        SuffixRecord(SuffixRecord &&other) noexcept = default;
        SuffixRecord &operator=(const SuffixRecord &) = delete;
        size_t getMaxSuffixLength() const {return max_suffix_len;}
        std::string str() const;
        const_iterator begin() const { return paths.begin(); }
        const_iterator end() const { return paths.end(); }
        iterator begin() { return paths.begin(); }
        iterator end() { return paths.end(); }

        size_t getNumberOfPaths() const;

        size_t countStartsWith(const GraphPath &path) const;

        bool empty() const;
        const Storage &getSuffixes() const {return paths;}
        Edge &getEdge() const {return *eid;}
//        Storage &getSuffixes() {return paths;}
    };

    inline std::ostream &operator<<(std::ostream &os, const SuffixRecord &rec) { return os << rec.str(); }

//For each vertex this structure stores subpaths of reads that start in the vertex
//For each occurence of vertex in a read only one subpath is stored
//The subpath is chosen as the shortest path such that total length of edges, starting from the second is at least max_length
//If read stops before subpath of required length is found, read suffix is stored regardless of its length
//For current implementation min_len=0 .
    class SuffixTracker : public AlignedReadStorageListener, public ResolutionListener {
    protected:
        // This map stores unique_ptrs by value and moves them around.
        // It can provide concurrent read/write with local locks
        libcuckoo::cuckoohash_map<ConstEdgeId, std::unique_ptr<SuffixRecord>> edge_data;
    public:
        typedef std::pair<const ConstEdgeId, std::unique_ptr<SuffixRecord>> EdgeDataUnit;

        AlignedReadStorage *storage;
        size_t min_len;
        size_t max_len;

    private:
//        processPath can only be called when vertex set can not be changed, so no graph modification.
//        The suffix stored for an edge of the path is the codes after it, up to the max-length limit and never past cut,
//        the logical position in the codes' deck where they ended when the read got its path (AlignedRead's suffix
//        ends). Codes appended past cut by later operations do not count, so adding and removing a read always compute
//        the same suffixes, and removal is an exact match.
        void processPath(PathPosition left, PathPosition right, NuclDeck::IteratorInfo cut, int diff);
//        Both directions of path, cut at fend in its fsplits and rend in its rsplits.
        void processPath(const GraphPath &path, NuclDeck::IteratorInfo fend, NuclDeck::IteratorInfo rend, int diff);
//        The read's path, cut at the read's suffix ends.
        void processRead(const AlignedRead &read, int diff) {
            processPath(read.getPath(), read.getSuffixFEnd(), read.getSuffixREnd(), diff);
        }
//        The edge's record, nullptr if it has none (one lookup).
        SuffixRecord *findRecord(const Edge &edge) const {
            SuffixRecord *res = nullptr;
            edge_data.find_fn(edge.getId(), [&res](const std::unique_ptr<SuffixRecord> &ptr) { res = ptr.get(); });
            return res;
        }
//        Exactly the virtually suffix edges (and the untyped edges of a DBG) hold records. Operations move a record to
//        each new edge that must hold one, from the edge it takes over from, which then holds one too; nothing is looked
//        up to decide.
        static bool holdsRecord(const Edge &edge) {
            return edge.isVirtualSuffix() || edge.getVirtualType() == VirtualType::unset;
        }
//        Gives from's record, which must exist, to `to`: whole or, with through, only the suffixes continuing through that
//        edge, with its code cut off. moveRecord takes the paths out of from's record, so it is only for a source
//        deleted by the operation (which erases its record); copyRecord leaves it intact, for a source that feeds
//        several edges.
        void moveRecord(const Edge &from, Edge &to, const Edge *through = nullptr);
        void copyRecord(const Edge &from, Edge &to, const Edge *through = nullptr);
//        The general form: with a non-empty cut, only the suffixes starting with cut, without it; cut_len is subtracted
//        from the maximal suffix length.
        void transferRecord(const Edge &from, Edge &to, bool copy, const Sequence &cut, size_t cut_len);

    public:
        SuffixTracker(AlignedReadStorage &storage, AssemblyGraph &graph, size_t _min_len, size_t _max_len);
        void fillFromStorage(logging::Logger &logger, size_t threads);

        SuffixTracker &operator=(SuffixTracker &&other) noexcept = default;
        SuffixTracker(SuffixTracker &&other)  noexcept = default;
        SuffixTracker &operator=(const SuffixTracker &other) = delete;
        SuffixTracker(const SuffixTracker &other) = delete;

//        Records exist on every edge for which holdsRecord holds.
        bool hasSuffixRecord(const Edge &edge) const {return edge_data.contains(edge.getId());}
        SuffixRecord &getSuffixRecord(const Edge &edge) {
            SuffixRecord *res = findRecord(edge);
            VERIFY(res != nullptr);
            return *res;
        }
        const SuffixRecord &getSuffixRecord(const Edge &edge) const {
            const SuffixRecord *res = findRecord(edge);
            VERIFY(res != nullptr);
            return *res;
        }
        size_t getMinLen() const { return min_len; }
        size_t getMaxLen() const { return max_len; }

        std::function<std::string(const Edge &)> labeler() const;

//        Subtracting every stored read must empty every record; the reads are added back afterwards.
        bool fireCheckConsistency(std::ostream &os) override;

        void fireAddRead(const AlignedRead &read) override;
        void fireRerouteRead(AlignedRead &read) override;
        void fireInvalidateRead(AlignedRead &read) override;

        void fireDeleteEdge(Edge &edge) override;
        void fireReplacePath(const PathReplacement &replacement) override;
//        Every edge that holdsRecord has a record.
        bool checkConsistency(std::ostream &os) override;
//        Only for supregraph edges (e.g. from ExtendAdjacent): converting a DBG (untyped edges) with a tracker attached is
//        not allowed, since its stored suffixes depend on the edges' virtual types.
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireMergeLoop(const ag::GraphPath  &path, Vertex &vertex) override {}
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override;

        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) override;
//        This methods works only for infinite extension holding. Need to rewrite to holding info in edges
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;;
    };
}
