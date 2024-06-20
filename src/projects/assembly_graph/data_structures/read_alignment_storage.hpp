#pragma once

#include "assembly_graph/graph_listeners.hpp"
#include "assembly_graph/aligned_read.hpp"
#include "aligned_read_listeners.hpp"
#include <common/logging.hpp>
#include <common/omp_utils.hpp>

namespace ag {

    template<class Traits>
    class AlignedReadStorageMaintenance;

    template<class Traits>
    class AlignedReadStorage : public AlignedReadStorageFire<Traits> {
        friend class AlignedReadStorageMaintenance<Traits>;
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Edge::ConstEdgeId ConstEdgeId;
        typedef typename Vertex::VertexId VertexId;
    private:
        std::vector <AlignedRead<Traits>> reads;
        std::unordered_map <ConstEdgeId, std::vector<AlignedReadDirection<Traits>>> starts;
        omp_lock_t writelock = {};
        ag::AlignedReadStorageMaintenance<Traits> * maintenance = nullptr;

    public:
        void lock() { omp_set_lock(&writelock); }
        void unlock() { omp_unset_lock(&writelock); }

        bool checkConsistency();

        AlignedReadStorage() = default;
        AlignedReadStorage(AlignedReadStorage &&other)  noexcept;
        AlignedReadStorage(const AlignedReadStorage &other) = delete;
        AlignedReadStorage &operator=(AlignedReadStorage &&other) noexcept;
        AlignedReadStorage &operator=(const AlignedReadStorage &other) = delete;

        AlignedReadStorage(AssemblyGraph<Traits> &graph, std::vector<AlignedRead<Traits>> reads);
        AlignedReadStorage(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph,
                           std::vector<AlignedRead<Traits>> reads);

        virtual ~AlignedReadStorage() {delete maintenance;}

        const std::vector<AlignedReadDirection<Traits>> &getOutgoingReadsLockFree(Edge &edge) const {
            return starts.at(edge.getId());
        }

//        TODO: Make this free of global lock
        const std::vector<AlignedReadDirection<Traits>> &getOutgoingReads(Edge &edge) const {
            lock();
            const std::vector<AlignedReadDirection<Traits>> & res = getOutgoingReadsLockFree(edge);
            unlock();
            return res;
        }

        std::vector<AlignedReadDirection<Traits>> &getOutgoingReadsLockFree(Edge &edge) {
            return starts.at(edge.getId());
        }

        std::vector<AlignedReadDirection<Traits>> &getOutgoingReads(Edge &edge) {
            lock();
            std::vector<AlignedReadDirection<Traits>> & res = getOutgoingReadsLockFree(edge);
            unlock();
            return res;
        }

//        typename std::vector<AlignedRead<Traits>>::iterator begin() { return reads.begin(); }
//        typename std::vector<AlignedRead<Traits>>::iterator end() { return reads.end(); }
        typename std::vector<AlignedRead<Traits>>::const_iterator begin() const { return reads.begin(); }
        typename std::vector<AlignedRead<Traits>>::const_iterator end() const { return reads.end(); }
        typename std::vector<AlignedRead<Traits>>::iterator begin() { return reads.begin(); }
        typename std::vector<AlignedRead<Traits>>::iterator end() { return reads.end(); }
//        AlignedRead<Traits> &operator[](size_t ind) { return reads[ind]; }
        const AlignedRead<Traits> &operator[](size_t ind) const { return reads[ind]; }
        AlignedRead<Traits> &operator[](size_t ind) { return reads[ind]; }
        size_t size() const { return reads.size(); }

        void updateStart(AlignedReadDirection<Traits> dir);
        size_t startCnt(const Edge &edge) const {return starts.at(edge.getId()).size();}

        void delayedInvalidateRead(AlignedRead<Traits> &read, const std::string &message);
        void rerouteRead(AlignedRead<Traits> &alignedRead, GraphPath<Traits> corrected, const string &message);
        bool apply(AlignedRead<Traits> &alignedRead);

        void applyCorrections(logging::Logger &logger, size_t threads);

        void printReadPaths(logging::Logger &logger, const std::experimental::filesystem::path &aln_path,
                            const std::experimental::filesystem::path &gfa_path,
                            const std::experimental::filesystem::path &rp_path,
                            size_t k) const;

        [[maybe_unused]] void printReadAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const;
        void printReadFasta(logging::Logger &logger, const std::experimental::filesystem::path &path) const;
        void printFullAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const;

        void Save(const std::experimental::filesystem::path &path);
        void Save(std::ostream &os) const;
        static AlignedReadStorage Load(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &path,
                                       AssemblyGraph<Traits> &graph, const IdIndex<Vertex> &index);
        static AlignedReadStorage Load(logging::Logger &logger, size_t threads, std::istream &is,
                                       AssemblyGraph<Traits> &graph, const IdIndex<Vertex> &index);

    };

//    TODO: move to dbg
    template<class Traits>
    class AlignedReadStorageMaintenance : public AlignedReadStorageListener<Traits>, public ResolutionListener<Traits> {
        friend class AlignedReadStorage<Traits>;
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Vertex::VertexId VertexId;
    private:
        AlignedReadStorage<Traits> *storage;
    public:
        AlignedReadStorageMaintenance(AssemblyGraph<Traits> &graph, AlignedReadStorage<Traits> &storage);
        AlignedReadStorageMaintenance(AlignedReadStorageMaintenance &&other) = default;

        void fireAddEdge(Edge &edge) override {
            storage->lock();
            storage->starts[edge.getId()] = {};
            storage->unlock();
        }
        void fireDeleteEdge(Edge &edge) override {
            storage->lock();
            storage->starts.erase(edge.getId());
            storage->unlock();
        }
        void fireAddSupreVertex(Vertex &v, Edge &e) {}
        void fireMergePath(const std::vector<EdgeId> &path, Vertex &new_vertex) override;
//        TODO: implement properly
        void fireMergeLoop(const GraphPath<Traits> &path, Vertex &new_vertex) override { VERIFY(false); };
        void fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) override;
        void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) override;

        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) override;

        void fireAddRead(const AlignedRead<Traits> &read) override;
        void fireRerouteRead(AlignedRead<Traits> &read) override;
        void fireInvalidateRead(AlignedRead<Traits> &read) override;
    };


    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireResolveVertex(Vertex &core, const VertexResolutionResult <Traits> &resolution) {
        for(Edge &edge : core) {
            std::vector<AlignedReadDirection<Traits>> &old_edge_rec = storage->getOutgoingReads(edge);
            for(AlignedReadDirection<Traits> &dir : old_edge_rec) {
                dir.getRead().lock();
                VERIFY(dir.frontEdge().isPrefix());
                dir.pop_front();
                if(!dir.empty()) {
                    storage->getOutgoingReads(dir.frontEdge()).emplace_back(dir);
                } else {
//                    TODO: account for reads contained within vertices.
                }
                dir.getRead().unlock();
            }
        }
        for(Edge &edge : core.incoming()) {
            std::vector<AlignedReadDirection<Traits>> &old_edge_rec = storage->getOutgoingReads(edge);
            for(AlignedReadDirection<Traits> &dir : old_edge_rec) {
                dir.getRead().lock();
                if(!dir.isSingleton() && !dir.empty()) {
                    Edge &new_start = resolution.get(dir.frontEdge(), (dir.firstPosition() + 1).nextEdge()).rc().front().rc();
                    storage->getOutgoingReads(new_start).emplace_back(dir);
                }
                dir.getRead().unlock();
            }
        }
    }

    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) {
        size_t left_skip = 0;
        GraphPath<Traits> prefix(path.front()->getStart());
        std::vector<AlignedReadDirection<Traits>> &new_edge_recs = storage->getOutgoingReads(new_edge);
        for(EdgeId eid : path) {
            for(const AlignedReadDirection<Traits> &dir : storage->getOutgoingReads(*eid)) {
                dir.getRead().lock();
                VERIFY(dir.valid())
                VERIFY(dir.getStart() == eid->getStart());
                if (!prefix.empty()) {
                    size_t cut_left = left_skip + dir.leftCut();
                    dir.setCutLeft(0);
                    dir.push_front(prefix);
                    dir.setCutLeft(cut_left);
                }
                new_edge_recs.emplace_back(dir);
                dir.getRead().unlock();
            }
            prefix += *eid;
            left_skip += eid->rc().truncSeq().size();
        }
    }

//TODO: shorten paths instead of prolonging. This will only be possible if either outgoing reads are stored in vertices
// or suffixes are not stored for vertices with outgoing prefix edges
    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireMergePath(const std::vector<EdgeId> &path, Vertex &new_vertex) {
        VERIFY(new_vertex != path.front()->getStart() && new_vertex != path.back()->getFinish());
        size_t left_skip = 0;
        Edge &first = new_vertex.rc().front().rc();
        std::vector<AlignedReadDirection<Traits>> &recs = storage->getOutgoingReads(first);
        GraphPath<Traits> prefix(path.front()->getStart());
        for(EdgeId edge : path) {
            for(AlignedReadDirection<Traits> dir : storage->getOutgoingReads(*edge)) {
                size_t cut_left = left_skip + dir.leftCut();
                dir.setCutLeft(0);
                dir.push_front(prefix);
                dir.setCutLeft(cut_left);
                VERIFY(dir.frontEdge() == *path.front());
                recs.emplace_back(dir);
            }
            prefix += *edge;
            left_skip += edge->rc().truncSeq().size();
        }
    }

    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireAddRead(const AlignedRead<Traits> &read) {
        VERIFY_MSG(false, "New reads can not be added when maintainance is already activated");
    }

    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireRerouteRead(AlignedRead<Traits> &read) {
        storage->updateStart(read.forward());
        storage->updateStart(read.backward());
    }

    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireInvalidateRead(AlignedRead<Traits> &read) {
        storage->updateStart(read.forward());
        storage->updateStart(read.backward());
    }
    template<class Traits>

    AlignedReadStorageMaintenance<Traits>::AlignedReadStorageMaintenance(AssemblyGraph<Traits> &graph,
                                                                   AlignedReadStorage<Traits> &storage) :
            AlignedReadStorageListener<Traits>(storage, "AlignedReadStorageMaintenance"), ResolutionListener<Traits>(graph, "AlignedReadStorageMaintenance"), storage(&storage) {
        for(Edge &edge : graph.edges())
            storage.starts[edge.getId()] = {};
        for(AlignedRead<Traits> &read: storage) {
            if(!read.getPath().empty()) {
                storage.starts[read.getPath().frontEdge().getId()].emplace_back(read.forward());
                storage.starts[read.getPath().backEdge().rc().getId()].emplace_back(read.backward());
            } else if(read.valid()) {
//                TODO: implement storing paths fully contained in vertices
            }
        }
    }

    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) {
        VERIFY(split.size() > 1);
        std::unordered_map<EdgeId, std::vector<AlignedReadDirection<Traits>> *> new_recs;
        for(EdgeId eid : split)
            new_recs[eid] = &storage->getOutgoingReads(*eid);
        for(AlignedReadDirection<Traits> direction : storage->getOutgoingReads(edge)) {
//            TODO: Is this lock necessary? We only ever interact with ends of reads and it looks ok to interact with
//             different ends simultaneously. The only interactions are conditions on Splits.size() (which looks ok),
//             and remapping of Splits (which may never actually be triggered in this code if we never increase the
//             size of Splits. Same question can be posed for merge and maybe others.
            direction.getRead().lock();
            EdgeId new_start;
//            This condition takes care of handling singleton paths when this listener is called for rc edge/path
            if(direction.getStart() == edge.getStart() && direction.getFinish() != split.front()->getFinish()) {
//                All directions extending beyond edge should be processed carefully to avoid spoiling NuclDeck
//                iterators stored in AlignedRead
                if (direction.getFSplits().size() > edge.getCode().size()) {
                    VERIFY(direction.leftCut() >= edge.truncSize() - split.back()->truncSize());
                    size_t left_cut = direction.leftCut();
                    direction.setCutLeft(0);
                    new_start = split.back();
                    direction.replace_front(*split.back());
                    direction.setCutLeft(left_cut + split.back()->truncSize() - edge.truncSize());
                } else {
                    Segment<Edge> seg = direction.front();
                    size_t start_pos = 0;
                    for (EdgeId eid: split) {
                        if (seg.left < start_pos + eid->truncSize()) {
                            VERIFY(seg.right <= start_pos + eid->truncSize());
                            VERIFY(new_start != split.back());
                            direction.setPath({*eid, seg.cutLeft() - start_pos, eid->truncSize() - (seg.cutLeft() - start_pos + seg.size())});
                            new_start = eid;
                            break;
                        }
                        start_pos += eid->truncSize();
                    }
                }
            } else {
                new_start = direction.getStart() == edge.getStart() ? split.front() : direction.frontEdge().getId();
            }
            direction.getRead().unlock();
            new_recs[new_start]->emplace_back(direction);
        }
    }

    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                                                 const AlignmentForm &left_al,
                                                                 const AlignmentForm &right_al) {
        VERIFY(left_al.targetLength() == right_al.targetLength());
        size_t left_skip = left.fullSize() - left_al.queryLength();
        size_t right_skip = right.fullSize() - right_al.queryLength();
        std::vector<AlignedReadDirection<Traits>> &new_rec = storage->getOutgoingReads(new_edge);
        std::vector<AlignedReadDirection<Traits>> &new_rec_rc = storage->getOutgoingReads(new_edge.rc());
//        Reads on left edge will be handled by rc
        for(AlignedReadDirection<Traits> dir : storage->getOutgoingReads(right)) {
            dir.getRead().lock();
            size_t new_left = dir.leftCut() >= right_al.queryLength() ?
                    left_skip + right_al.targetLength() + dir.leftCut() - right_al.queryLength() :
                    left_skip + right_al.lastColumnByQpos(dir.leftCut()).getTpos();
            if(dir.isSingleton()) {
                size_t old_right = right.fullSize() - dir.rightCut();
                size_t new_right = old_right >= right_al.queryLength() ?
                                   left_skip + right_al.targetLength() + old_right - right_al.queryLength() :
                                   left_skip + right_al.firstColumnByQpos(old_right).getTpos();
                dir.setCutRight(new_edge.fullSize() - new_right);
                new_rec_rc.emplace_back(dir.RC());
            }
            dir.setCutLeft(0);
            dir.forcePushFront(left);// While this operation is underway this path is disconnected.
            dir.setCutLeft(new_left);
            dir.getRead().unlock();
            new_rec.emplace_back(dir);
        }
    }

    template<class Traits>
    void AlignedReadStorageMaintenance<Traits>::fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> & graph) {
#pragma omp parallel for default(none) schedule(dynamic, 100)
        for(size_t i = 0; i < storage->size(); i++) {
            AlignedRead<Traits> &read = (*storage)[i];
            read.resetEdgeCodes();
        }
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::delayedInvalidateRead(AlignedRead<Traits> &read, const string &message) { // NOLINT(readability-convert-member-functions-to-static)
        read.delayedInvalidate();
        this->fireDelayedInvalidateRead(read, message);
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::rerouteRead(AlignedRead<Traits> &alignedRead, GraphPath <Traits> corrected,
                                          const string &message) {
        VERIFY(corrected.truncLen() >= 500);
        alignedRead.correct(std::move(corrected));
        this->fireDelayedRerouteRead(alignedRead, message);
    }

    template<class Traits>
    bool AlignedReadStorage<Traits>::apply(AlignedRead<Traits> &alignedRead) {
        if (!alignedRead.checkCorrected())
            return false;
        if(alignedRead.getCorrected().valid())
            this->fireRerouteRead(alignedRead);
        else
            this->fireInvalidateRead(alignedRead);
        alignedRead.applyCorrection();
        return true;
    }

    template<class Traits>
    [[maybe_unused]] [[maybe_unused]] void AlignedReadStorage<Traits>::printReadAlignments(logging::Logger &logger,
                                                  const std::experimental::filesystem::path &path) const {
        logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
        std::string acgt = "ACGT";
        std::ofstream os;
        os.open(path);
        for (const AlignedRead<Traits> &read: reads) {
            const GraphPath<Traits> &al = read.getPath();
            if (!al.valid())
                continue;
            os << read.id << " " << al.getStart().getInnerId()
               << " " << al.cpath().str() << "\n";
            GraphPath<Traits> rc_al = al.RC();
            os << "-" << read.id << " " << rc_al.getStart().getInnerId() << " " << rc_al.cpath().str() << "\n";
        }
        os.close();
    }

    template<class Traits>
    void
    AlignedReadStorage<Traits>::printReadFasta(logging::Logger &logger, const std::experimental::filesystem::path &path) const {
        logger.info() << "Printing reads to fasta file " << path << std::endl;
        std::string acgt = "ACGT";
        std::ofstream os;
        os.open(path);
        for (const AlignedRead<Traits> &read: reads) {
            const GraphPath<Traits> &al = read.getPath();
            if (!al.valid())
                continue;
            os << ">" << read.getId() << "\n" << read.getPath().Seq() << "\n";
        }
        os.close();
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::printReadPaths(logging::Logger &logger,
                                               const std::experimental::filesystem::path &aln_path,
                                               const std::experimental::filesystem::path &gfa_path,
                                               const std::experimental::filesystem::path &rp_path,
                                               size_t k) const {
        logger.info() << "Printing reads paths to file " << aln_path << std::endl;
        std::ofstream os;
        os.open(aln_path);
        Save(os);
        os.close();
        os.open(rp_path);
        os << gfa_path.c_str() << std::endl;
        os << aln_path.c_str() << std::endl;
        os << k << std::endl;
        os.close();
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::printFullAlignments(logging::Logger &logger,
                                                  const std::experimental::filesystem::path &path) const {
        logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
        std::ofstream os;
        os.open(path);
        for (const AlignedRead<Traits> &read: reads) {
            const GraphPath<Traits> &al = read.getPath();
            if (!al.valid())
                continue;
            os << read.getId() << " " << read.getPath().str() << "\n";
            os << "-" << read.getId() << " " << read.getPath().RC().str() << "\n";
        }
        os.close();
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::applyCorrections(logging::Logger &logger, size_t threads) {
        if (size() > 10000)
            logger.info() << "Applying corrections to reads" << std::endl;
        omp_set_num_threads(int(threads));
        ParallelCounter cnt(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(cnt)
        for (size_t i = 0; i < reads.size(); i++) { // NOLINT(modernize-loop-convert)
            if (apply(reads[i]))
                cnt += 1;
        }
        this->fireAppliedCorrections(cnt.get());
        if (size() > 10000)
            logger.info() << "Applied correction to " << cnt.get() << " reads" << std::endl;
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::Save(std::ostream &os) const {
        os << size() << "\n";
        for (const AlignedRead<Traits> &alignedRead: *this) {
            os << alignedRead << "\n";
        }
    }

    template<class Traits>
    AlignedReadStorage<Traits> AlignedReadStorage<Traits>::Load(logging::Logger &logger, size_t threads,
                                                                std::istream &is, AssemblyGraph<Traits> &graph,
                                                                const IdIndex<Vertex> &index) {
        size_t sz;
        is >> sz;
        std::vector<AlignedRead<Traits>> reads;
        for (size_t i = 0; i < sz; i++) {
            reads.emplace_back(AlignedRead<Traits>::Load(is, index));
        }
        return {logger, threads, graph, std::move(reads)};
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::updateStart(AlignedReadDirection<Traits> dir) {
        if(dir.valid() && dir.getCorrected().valid() && dir.frontEdge() == dir.getCorrected().frontEdge())
            return;
        if(dir.valid()) {
            std::vector<AlignedReadDirection<Traits>> &old = this->getOutgoingReads(dir.frontEdge());
            dir.getStart().lock();
            old.erase(std::find(old.begin(), old.end(), dir));
            dir.getStart().unlock();
        }
        if(dir.getCorrected().valid()) {
            auto &rec = this->getOutgoingReads(dir.getCorrected().frontEdge());
            Vertex &v = dir.getCorrected().getStart();
            v.lock();
            rec.emplace_back(dir);
            v.unlock();
        }
    }

    template<class Traits>
    bool AlignedReadStorage<Traits>::checkConsistency() {
        for(AlignedRead<Traits> &al : reads) {
            if(!al.valid())
                continue;
            VERIFY(starts.find(al.getPath().frontEdge().getId()) != starts.end());
            std::vector<AlignedReadDirection<Traits>> &start = starts.at(al.getPath().frontEdge().getId());
            bool f1 = false;
            for(AlignedReadDirection<Traits> &dir : start) {
                if(dir.getRead().getId() == al.getId()) {
                    f1 = true;
                    break;
                }
            }
            VERIFY(starts.find(al.getPath().backEdge().rc().getId()) != starts.end());
            std::vector<AlignedReadDirection<Traits>> &rcstart = starts.at(al.getPath().backEdge().rc().getId());
            bool f2 = false;
            for(AlignedReadDirection<Traits> &dir : rcstart) {
                if(dir.getRead().getId() == al.getId()) {
                    f2 = true;
                    break;
                }
            }
            if(!f1 || !f2) {
                VERIFY_MSG(false, al.getId() << " " << al.getPath().getStart().getId() << " " << al.getPath().getFinish().getId() << " " << al.getPath().getFSplits().str() << " " << al.getPath().getRSplits().str());
                return false;
            }
        }
        return this->fireCheckConsistency();
    }

    template<class Traits>
    AlignedReadStorage<Traits>::AlignedReadStorage(logging::Logger &logger, size_t threads,
                                                   AssemblyGraph<Traits> &graph, std::vector<AlignedRead<Traits>> reads)
            : reads(std::move(reads)) {
        omp_set_num_threads(int(threads));
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(reads)
        for(size_t i = 0; i < reads.size(); i++) {
            this->fireAddRead(reads[i]);
        }
        maintenance = new ag::AlignedReadStorageMaintenance<Traits>(graph, *this);
    }

    template<class Traits>
    AlignedReadStorage<Traits>::AlignedReadStorage(AlignedReadStorage &&other) noexcept {
        *this = std::move(other);
    }

    template<class Traits>
    AlignedReadStorage<Traits> &AlignedReadStorage<Traits>::operator=(AlignedReadStorage &&other) noexcept {
        AlignedReadStorageFire<Traits>::operator=(std::move(other));
        std::swap(reads, other.reads);
        std::swap(starts, other.starts);
        std::swap(writelock, other.writelock);
        std::swap(maintenance, other.maintenance);
        if(maintenance != nullptr)
            maintenance->storage = this;
        if(other.maintenance != nullptr)
            other.maintenance->storage = &other;
        return *this;
    }

    template<class Traits>
    AlignedReadStorage<Traits>::AlignedReadStorage(AssemblyGraph<Traits> &graph, std::vector<AlignedRead<Traits>> reads)
            : reads(std::move(reads)) {
        for(size_t i = 0; i < reads.size(); i++) {
            this->fireAddRead(reads[i]);
        }
        maintenance = new ag::AlignedReadStorageMaintenance<Traits>(graph, *this);
    }

    template<class Traits>
    void AlignedReadStorage<Traits>::Save(const std::experimental::filesystem::path &path) {
        std::ofstream os;
        os.open(path);
        Save(os);
        os.close();
    }

    template<class Traits>
    AlignedReadStorage<Traits> AlignedReadStorage<Traits>::Load(logging::Logger &logger, size_t threads,
                                                        const std::experimental::filesystem::path &path,
                                                        AssemblyGraph<Traits> &graph, const IdIndex<Vertex> &index) {
        std::ifstream is;
        is.open(path);
        AlignedReadStorage<Traits> result(Load(logger, threads, is, graph, index));
        is.close();
        return std::move(result);
    }

    template<class Traits>
    void SaveReads(const std::experimental::filesystem::path &fname, const AlignedReadStorage<Traits> &storage) {
        std::ofstream os;
        os.open(fname);
        storage.Save(os);
        os.close();
    }

    template<class Traits>
    AlignedReadStorage<Traits> LoadReads(logging::Logger &logger, size_t threads,
                                         const std::experimental::filesystem::path &fname,
                                         const IdIndex<typename Traits::Vertex> &index) {
        std::ifstream is;
        is.open(fname);
        AlignedReadStorage<Traits> res = AlignedReadStorage<Traits>::Load(logger, threads, is, index);
        is.close();
        return std::move(res);
    }
}
