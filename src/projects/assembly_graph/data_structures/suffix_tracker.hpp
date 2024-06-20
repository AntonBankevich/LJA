#pragma once
#include "assembly_graph/assembly_graph.hpp"
#include "read_alignment_storage.hpp"
#include "read_logger.hpp"
#include "fstream"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include <experimental/filesystem>

namespace ag {

    template<class Traits>
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
    template<class Traits>
    struct SuffixRecord {
        friend SuffixTracker<Traits>;
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Traits::Edge Edge;
        typedef typename Edge::EdgeId EdgeId;
    private:
        typedef std::vector<std::pair<Sequence, int>> Storage;
        typedef Storage::const_iterator const_iterator;
        typedef Storage::iterator iterator;
        EdgeId eid;
        Storage paths;
        size_t zero_cnt = 0;
        size_t max_suffix_len;
        int num_of_ends = 0;

        void lock() const { eid->getStart().lock(); }
        void unlock() const { eid->getStart().unlock(); }
        void updateZero(size_t old_val, size_t new_val);
        void lockFreeChangePathCnt(const Sequence &seq, int diff);
        void changePathCnt(const Sequence &seq, int diff);
        void addPath(const Sequence &seq, int diff = 1);
        void removePath(const Sequence &seq) {changePathCnt(seq, -1);}
        void directAddPath(const Sequence &seq, size_t cnt) {paths.emplace_back(seq, cnt);}
        void clear() { paths.clear(); }
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
        size_t countStartsWith(const GraphPath<Traits> &path) const {
            VERIFY(path.empty() || path.getStart() == eid->getFinish());
            if(path.empty()) return countStartsWith(Sequence());
            return countStartsWith(Sequence(path.getFSplits().begin(), path.getFSplits().end() - path.backEdge().getCode().size() + 1));
        }
        bool empty() const;
        const Storage &getSuffixes() const {return paths;}
        Edge &getEdge() const {return *eid;}
//        Storage &getSuffixes() {return paths;}
    };

    template<class Traits>
    inline std::ostream &operator<<(std::ostream &os, const SuffixRecord<Traits> &rec) { return os << rec.str(); }

//For each vertex this structure stores subpaths of reads that start in the vertex
//For each occurence of vertex in a read only one subpath is stored
//The subpath is chosen as the shortest path such that total length of edges, starting from the second is at least max_length
//If read stops before subpath of required length is found, read suffix is stored regardless of its length
    template<class Traits>
    class SuffixTracker : public AlignedReadStorageListener<Traits>, public ResolutionListener<Traits> {
    protected:
        std::unordered_map<typename Traits::Edge::ConstEdgeId, SuffixRecord<Traits>> edge_data;
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge::EdgeId EdgeId;
        typedef typename Traits::Vertex::VertexId VertexId;
        typedef typename Traits::Vertex::ConstVertexId ConstVertexId;
        typedef std::pair<const typename Traits::Edge::ConstEdgeId, SuffixRecord<Traits>> EdgeDataUnit;

        AlignedReadStorage<Traits> *storage;
        size_t min_len;
        size_t max_len;

    private:
//        processPath can only be called when vertex set can not be changed, so no graph modification
        void
        processPath(PathPosition<Traits> left, PathPosition<Traits> right, int diff);
        void addSubpath(PathPosition<Traits> left, PathPosition<Traits> right) {processPath(left, right, 1);}
        void removeSubpath(PathPosition<Traits> left, PathPosition<Traits> right) {processPath(left, right, -1);}

    public:
        SuffixTracker(AlignedReadStorage<Traits> &storage, AssemblyGraph<Traits> &graph, size_t _min_len, size_t _max_len);
        void fillFromStorage(logging::Logger &logger, size_t threads);

        SuffixTracker &operator=(SuffixTracker &&other) noexcept = default;
        SuffixTracker(SuffixTracker &&other)  noexcept = default;
        SuffixTracker &operator=(const SuffixTracker &other) = delete;
        SuffixTracker(const SuffixTracker &other) = delete;

        const SuffixRecord<Traits> &getSuffixRecord(const Edge &edge) const {return edge_data.at(edge.getId());}
        size_t getMinLen() const { return min_len; }
        size_t getMaxLen() const { return max_len; }

        std::function<std::string(const Edge &)> labeler() const;

        bool fireCheckConsistency() override;

        void fireAddRead(const AlignedRead<Traits> &read) override;
        void fireRerouteRead(AlignedRead<Traits> &read) override;
        void fireInvalidateRead(AlignedRead<Traits> &read) override;

        void fireAddEdge(Edge &edge) override;
        void fireDeleteEdge(Edge &edge) override;
        void fireAddSupreVertex(Vertex &v, Edge &e);
        void fireMergePath(const std::vector<EdgeId> &path, Vertex &vertex) override;
        void fireMergeLoop(const ag::GraphPath <Traits> &path, Vertex &vertex) override { VERIFY(false);}
        void fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) override;

        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) override;
//        This methods works only for infinite extension holding. Need to rewrite to holding info in edges
        void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) override;;
    };


    template<class Traits>
    bool SuffixRecord<Traits>::empty() const {
        if(paths.empty())
            return true;
        for(const auto &path: paths) {
            if(path.second > 0)
                return false;
        }
        return true;
    }

    template<class Traits>
    size_t SuffixRecord<Traits>::countStartsWith(const Sequence &seq) const {
//        lock();
        int cnt = 0;
        for (const std::pair<Sequence, int> &rec: paths) {
            if (rec.first.startsWith(seq)) {
                cnt += rec.second;
            }
        }
        if(seq.empty())
            cnt += num_of_ends;
//        unlock();
        VERIFY(cnt >= 0);
        return cnt;
    }

    template<class Traits>
    std::string SuffixRecord<Traits>::str() const {
        std::stringstream ss;
        lock();
        for (const auto &path: paths) {
            ss << path.first << " " << path.second << std::endl;
        }
        unlock();
        return ss.str();
    }

    template<class Traits>
    void SuffixRecord<Traits>::updateZero(size_t old_val, size_t new_val) {
        if(new_val == 0)
            zero_cnt++;
        if(old_val == 0)
            zero_cnt--;
    }

    template<class Traits>
    void SuffixRecord<Traits>::removeZero() {
        std::vector<std::pair<Sequence, int>> new_paths;
        for (std::pair<Sequence, int> &rec: this->paths) {
            if (rec.second != 0) {
                new_paths.emplace_back(std::move(rec.first), rec.second);
            }
        }
        std::swap(this->paths, new_paths);
        this->zero_cnt = 0;
    }

    template<class Traits>
    void SuffixRecord<Traits>::lockFreeChangePathCnt(const Sequence &seq, int diff) {
        if(seq.empty()) {
            VERIFY(num_of_ends + diff >= 0);
            num_of_ends += diff;
            return;
        }
        if (diff == 0) return;
        if(diff > 0) {
            addPath(seq, diff);
            return;
        }
        for (size_t i = 0; i < paths.size(); i++) {
            std::pair<Sequence, int> &path = paths[i];
            if(path.first == seq && path.second + diff >= 0) {
                path.second += diff;
                updateZero(path.second - diff, path.second);
                diff = 0;
                return;
            }
        }
        std::vector<size_t> subseqs;
        for (size_t i = 0; i < paths.size(); i++) {
            std::pair<Sequence, int> &path = paths[i];
            if(seq.startsWith(path.first)) {
                subseqs.push_back(i);
            }
        }
        if(diff > 0) {
            paths.emplace_back(seq, diff);
            if (zero_cnt > paths.size() / 3) {
                removeZero();
            }
            return;
        }
        std::sort(subseqs.begin(), subseqs.end(), [this](size_t a, size_t b){return paths[a].first.size() < paths[b].first.size();});
        while(!subseqs.empty() && diff != 0) {
            if(paths[subseqs.back()].second != 0) {
                if (diff + paths[subseqs.back()].second >= 0) {
                    paths[subseqs.back()].second += diff;
                    updateZero(1, paths[subseqs.back()].second);
                    diff = 0;
                } else {
                    diff += paths[subseqs.back()].second;
                    paths[subseqs.back()].second = 0;
                    updateZero(1, 0);
                }
            }
            subseqs.pop_back();
        }
        VERIFY_MSG(diff == 0, "Attempting to remove path that is not present in the record.");
        if (zero_cnt > paths.size() / 3) {
            removeZero();
        }
    }

    template<class Traits>
    void SuffixRecord<Traits>::changePathCnt(const Sequence &seq, int diff) {
        lock();
        lockFreeChangePathCnt(seq, diff);
        unlock();
    }

    template<class Traits>
    void SuffixRecord<Traits>::addPath(const Sequence &seq, int diff) {
        VERIFY(diff > 0);
        if(seq.empty()) {
            num_of_ends += diff;
            return;
        }
        for (size_t i = 0; i < paths.size(); i++) {
            std::pair<Sequence, int> &path = paths[i];
            if(seq == path.first) {
                updateZero(path.second, 1);
                path.second += diff;
                return;
            }
        }
        paths.emplace_back(seq, diff);
    }

    template<class Traits>
    void SuffixRecord<Traits>::resetCodes(Vertex &start) {
        std::vector<std::pair<Sequence, int>> old = std::move(paths);
        for(const std::pair<Sequence, int> &rec : old) {
            if(rec.second == 0)
                continue;
            GraphPath<Traits> path(start, rec.first);
            path.resetEdgeCodes();
            changePathCnt(Sequence(path.getFSplits()), rec.second);
        }
    }


    template<class Traits>
    std::function<std::string(const typename Traits::Edge & )> SuffixTracker<Traits>::labeler() const {
        return [this](const typename Traits::Edge &edge) {
            const SuffixRecord<Traits> &rec = getSuffixRecord(edge);
            std::stringstream ss;
            size_t cnt = 0;
            ss << "Ends: " << rec.num_of_ends << "\\n";
            for (const auto &ext: rec) {
                if (cnt < 30)
                    ss << ext.first << "(" << ext.second << ")\\n";
                cnt++;
            }
            if (cnt > 30) {
                ss << "and another " << (cnt - 30) << " records\\n";
            }
            return ss.str();
        };
    }

    template<class Traits>
    void SuffixTracker<Traits>::processPath(PathPosition<Traits> left, PathPosition<Traits> right, int diff) {
        PathPosition<Traits> from_pos = left;
        PathPosition<Traits> to_pos = from_pos + 1;
        size_t clen = left.nextEdge().truncSize();
        VERIFY(to_pos <= right);
        Sequence read_seq(left.getFPos(), right.getFPos());
        while (from_pos != right) {
            Edge &edge = from_pos.nextEdge();
            clen -= edge.truncSize();
            ++from_pos;
            SuffixRecord<Traits> &erec = edge_data.at(edge.getId());
            while (to_pos < right && (clen < erec.getMaxSuffixLength())) {
                Edge &new_edge = to_pos.nextEdge();
                clen += new_edge.truncSize();
                to_pos += new_edge;
            }
            VERIFY(to_pos <= right);
//            if(to_pos > right)
//                to_pos = right;
            if (clen >= min_len) {
                erec.changePathCnt(
                        read_seq.Subseq(from_pos.getFPos() - left.getFPos(), to_pos.getFPos() - left.getFPos()),
                        diff);
            }
        }
    }

    template<class Traits>
    void SuffixTracker<Traits>::fillFromStorage(logging::Logger &logger, size_t threads) {
        logger.info() << "Collecting and storing read suffixes" << std::endl;
        omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(storage) schedule(dynamic, 100)
        for (size_t i = 0; i < storage->size(); i++) {
            if(!(*storage)[i].getPath().empty())
                fireAddRead((*storage)[i]);
        }
        logger.info() << "Finished collecting and storing read suffixes" << std::endl;
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireAddRead(const AlignedRead<Traits> &read) {
        if(read.valid()) {
            addSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
            addSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
        }
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireRerouteRead(AlignedRead<Traits> &read) {
        if(read.valid()) {
            removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
            removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
        }
        if(read.getCorrected().valid()) {
            addSubpath(read.getCorrected().firstPosition(), read.getCorrected().lastPosition());
            addSubpath(read.getCorrected().lastPosition().RC(), read.getCorrected().firstPosition().RC());
        }
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireInvalidateRead(AlignedRead<Traits> &read) {
        if(read.valid()) {
            removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
            removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
        }
    }

    template<class Traits>
    SuffixTracker<Traits>::SuffixTracker(AlignedReadStorage<Traits> &storage, AssemblyGraph<Traits> &graph,
                                         size_t _min_len, size_t _max_len) :
            AlignedReadStorageListener<Traits>(storage, "SuffixTracker"), ResolutionListener<Traits>(graph, "SuffixTracker"), storage(&storage),
            min_len(_min_len), max_len(_max_len) {
        for(Edge &e : graph.edges()) {
            edge_data.emplace(std::piecewise_construct,
                         std::forward_as_tuple(e.getId()),
                         std::forward_as_tuple(e, max_len));
        }
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) {
        size_t ends = 0;
        for(EdgeId eid : path) {
            ends += edge_data.at(eid).num_of_ends;
        }
        SuffixRecord<Traits> &erec = edge_data.at(new_edge.getId());
        SuffixRecord<Traits> &lasterec = edge_data.at(path.back());
        erec.paths = std::move(lasterec.paths);
        erec.num_of_ends = ends;
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireMergePath(const std::vector<EdgeId> &path, Vertex &vertex) {
//        TODO: when we stop storing info for suffix edges, we can remove vector copy here
        size_t ends = 0;
        for(EdgeId eid : path) {
            ends += edge_data.at(eid).num_of_ends;
        }
        SuffixRecord<Traits> &lasterec = edge_data.at(path.back());
        edge_data.at(vertex.front().getId()).paths = lasterec.paths;
        SuffixRecord<Traits> &start_edge_rec = edge_data.at(vertex.rc().front().rc().getId());
        start_edge_rec.paths = std::move(lasterec.paths);
        start_edge_rec.num_of_ends = ends;
    }

    template<class Traits>
    void
    SuffixTracker<Traits>::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &,
                                               const AlignmentForm &) {
        SuffixRecord<Traits> &erec = edge_data.at(new_edge.getId());
        SuffixRecord<Traits> &rrec = edge_data.at(right.getId());
        SuffixRecord<Traits> &lrec = edge_data.at(left.getId());
        erec.paths = std::move(rrec.paths);
        erec.num_of_ends = lrec.num_of_ends + rrec.num_of_ends;
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) {
        SuffixRecord<Traits> &erec = edge_data.at(edge.getId());
        SuffixRecord<Traits> &last_rec = edge_data.at(split.back());
        last_rec.paths = std::move(erec.paths);
        for(EdgeId eid : split) {
            edge_data.at(eid->rc().getId()).num_of_ends = storage->startCnt(*eid);
        }
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
        std::function<void(size_t, Edge &)> task = [this](size_t, Edge &e) {
            edge_data.at(e.getId()).resetCodes(e.getFinish());
        };
        processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, task);
    }

    template<class Traits>
    bool SuffixTracker<Traits>::fireCheckConsistency() {
        for(auto &p : edge_data) {
            SuffixRecord<Traits> &rec = p.second;
            VERIFY(rec.num_of_ends == storage->startCnt(p.first->rc()));
        }
        for(AlignedRead<Traits> &read : *storage) {
            if(read.valid()) {
                removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
                removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
            }
        }
        for(auto &p : edge_data) {
            SuffixRecord<Traits> &rec = p.second;
            rec.removeZero();
            VERIFY(rec.begin() == rec.end());
            if(rec.begin() != rec.end())
                return false;
        }
        for(AlignedRead<Traits> &read : *storage) {
            if(read.valid()) {
                addSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
                addSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
            }
        }
        return true;
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireAddSupreVertex(Vertex &v, Edge &e) {
//        TODO: when we stop storing info for suffix edges, we can remove vector copy here
        SuffixRecord<Traits> &erec = edge_data.at(e.getId());
        edge_data.at(v.front().getId()).paths = erec.paths;
        edge_data.at(v.rc().front().rc().getId()).paths = std::move(erec.paths);
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireAddEdge(Edge &edge) {
//        TODO: get rid of these global locks
        storage->lock();
        edge_data.emplace(std::piecewise_construct,
                          std::forward_as_tuple(edge.getId()),
                          std::forward_as_tuple(edge, max_len));
        storage->unlock();
    }

    template<class Traits>
    void SuffixTracker<Traits>::fireDeleteEdge(Edge &edge) {
        storage->lock();
        edge_data.erase(edge.getId());
        storage->unlock();
    }

//    TODO: this code will break if paths can end with contracting edge. Need to make sure paths are normalized.
    template<class Traits>
    void SuffixTracker<Traits>::fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) {
        for (Edge &inc: core.incoming()) {
            const SuffixRecord<Traits> &rec = edge_data.at(inc.getId());
            std::unordered_map<EdgeId, SuffixRecord<Traits> *> rec_map;
            for(Edge &out : core) {
                if(!resolution.contains(inc, out))
                    continue;
                rec_map[out.getId()] = &edge_data.at(resolution.get(inc, out).front().getId());
                size_t len = edge_data.at(inc.getId()).max_suffix_len;
                len -= std::min(len, out.truncSize());
                rec_map[out.getId()]->max_suffix_len = len;
            }
            for (const std::pair<Sequence, size_t> &info: rec) {
                if (info.second == 0)
                    continue;
                VERIFY(!info.first.empty());
                Edge &out = core.getOutgoing(info.first[0]);
                if(out.getCode().size() < info.first.size()) {
                    rec_map.at(out.getId())->directAddPath(info.first.Subseq(out.getCode().size()), info.second);
                } else {
                    rec_map.at(out.getId())->num_of_ends += info.second;
                }
            }
        }
//                TODO: get rid of this when suffix edges do not store extensions
        for(Vertex &new_vertex : resolution.newVertices()) {
            SuffixRecord<Traits> &rec1 = edge_data.at(new_vertex.rc().front().rc().getId());
            SuffixRecord<Traits> &rec2 = edge_data.at(new_vertex.front().getId());
            rec1.paths = rec2.paths;
            rec1.max_suffix_len = rec2.max_suffix_len;
        }
    }

}
