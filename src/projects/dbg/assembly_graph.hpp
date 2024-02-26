#pragma once
#include "sequences/sequence.hpp"
#include "common/iterator_utils.hpp"
#include "common/object_id.hpp"
#include "id_index.hpp"
#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <array>
namespace ag {
    class BaseEdgeId {
    public:
        int vid;
        int eid;
        BaseEdgeId() = default;
        BaseEdgeId(int vid, int eid) : vid(vid), eid(eid) {
        }

        bool valid() const {
            VERIFY((vid != 0) == (eid != 0));
            return vid != 0;
        }

        bool operator==(const BaseEdgeId &other) const {return vid == other.vid && eid == other.eid;}
        bool operator!=(const BaseEdgeId &other) const {return vid != other.vid || eid != other.eid;}
        bool operator<(const BaseEdgeId &other) const {return vid < other.vid || (vid == other.vid && eid < other.eid);}
        bool operator>(const BaseEdgeId &other) const {return vid > other.vid || (vid == other.vid && eid > other.eid);}
        bool operator<=(const BaseEdgeId &other) const {return vid < other.vid || (vid == other.vid && eid <= other.eid);}
        bool operator>=(const BaseEdgeId &other) const {return vid > other.vid || (vid == other.vid && eid >= other.eid);}

        std::string str() const {
            return itos(vid) + "." +itos(eid);
        }
    };

    inline std::ostream &operator<<(std::ostream &os, ag::BaseEdgeId val) {
        return os << val.vid << "." << val.eid;
    }
}

template<>
inline ag::BaseEdgeId Parse<ag::BaseEdgeId>(const std::string &s) {
    int vid = ParseInt(s, 0);
    int eid = ParseInt(s, s.find('.') + 1);
    return {vid, eid};
}

namespace std {
    template <>
    class hash < ag::BaseEdgeId >{
    public:
        size_t operator()(const ag::BaseEdgeId &x ) const {
            return std::hash<int>()(x.vid * 12343251) ^ std::hash<int>()(x.eid * 1294835);
        }
    };
}


namespace ag {
    template<class T>
    class Locker {
    private:
        std::vector<T *> locked;
    public:
        explicit Locker(std::vector<T *> _vertices) : locked(std::move(_vertices)) {
            std::function<bool(T * const &v1, T * const &v2)> f = [](T * const &v1, T * const &v2) {
                return *v1 < *v2;
            };
            std::sort(locked.begin(), locked.end(), f);
            locked.erase(std::unique(locked.begin(), locked.end()), locked.end());
            for(T *v: locked) {
                v->lock();
            }
        }
        Locker(const Locker &other) = delete;
        Locker(Locker &&other) = delete;
        Locker& operator=(const Locker &other) = delete;
        Locker& operator=(Locker &&other) = delete;

        ~Locker() {
            for(T *v: locked) {
                v->unlock();
            }
        }
    };

    enum EdgeMarker {
        incorrect,
        suspicious,
        common,
        possible_break,
        correct,
        unique,
        repeat
    };

    inline bool IsMarkerCorrect(EdgeMarker marker) {
        return marker == EdgeMarker::correct || marker == EdgeMarker::unique || marker == EdgeMarker::repeat;
    }

    template<class Traits>
    class BaseVertex;
    template<class Traits>
    class BaseEdge;
    template<class Traits>
    class AssemblyGraph;



    template<class Traits>
    class BaseEdge {
    public:
        typedef BaseEdgeId id_type;
        typedef typename Traits::VertexData VertexData;
        typedef typename Traits::EdgeData EdgeData;
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef ObjectId<Edge, id_type> EdgeId;
        typedef ObjectId<const Edge, id_type> ConstEdgeId;
    private:
        id_type id;
        Vertex *start;
        Vertex *finish;
        Sequence seq;
        Edge *_rc;
        ag::EdgeMarker marker = ag::EdgeMarker::common;

        friend class BaseVertex<Traits>;

    public:
        friend class BaseVertex<Traits>;


        BaseEdge(id_type id, Vertex &_start, Vertex &_end, Sequence _seq) :
                id(std::move(id)), start(&_start), finish(&_end), seq(std::move(_seq)), _rc(nullptr) {
        }
//        BaseEdge() : id(), start(nullptr), finish(nullptr), seq() {}
        BaseEdge(const BaseEdge &other) = delete;

        virtual ~BaseEdge() = default;

        bool isCanonical() const {return *this <= rc();}

        EdgeId getId() {return {getInnerId(), dynamic_cast<Edge *>(this)};}
        ConstEdgeId getId() const {return {getInnerId(), dynamic_cast<const Edge *>(this)};}

        const Sequence &truncSeq() const { return seq; }
        const Sequence nuclLabel() const {
            if(truncSeq().empty())
                return {};
            return truncSeq().Subseq(0, 1);
        }
        Sequence kmerSeq(size_t pos) const {return fullSubseq(pos, pos);}


        const Vertex &getFinish() const {return *finish;}
        Vertex &getFinish() {return *finish;}
        const Vertex &getStart() const {return *start;}
        Vertex &getStart() {return *start;}

        void mark(ag::EdgeMarker _marker) { marker = _marker; };
        ag::EdgeMarker getMarker() const { return marker; };

//        This method should only be invoked if no graph modification is performed in parallel or if both start and
//        rc end vertices are locked by this process or otherwise prevented from modification by other processes

        Edge &rc() {
            VERIFY(!getStart().getSeq().empty());
            return *_rc;
        }

        const Edge &rc() const {
            VERIFY(!getStart().getSeq().empty());
            return *_rc;
        }

        size_t truncSize() const {
            return truncSeq().size();
        }

        bool operator==(const BaseEdge &other) const {
            return this == &other;
        }

        bool operator!=(const BaseEdge &other) const {
            return this != &other;
        }

        bool operator<(const BaseEdge &other) const {
            if(this == &other)
                return false;
            if(start != other.start)
                return *start < *other.start;
            return this->truncSeq() < other.truncSeq();
        }

        bool operator>(const BaseEdge &other) const {
            if(this == &other)
                return false;
            if(start != other.start)
                return *start > *other.start;
            return other.truncSeq() < truncSeq();
        }

        bool operator<=(const BaseEdge &other) const {
            return *this == other || *this < other;
        }

        std::string str() const {
            return getInnerId().str();
        }

        Sequence fullSubseq(size_t from, size_t to) const {
            VERIFY(start->size() > 0);
            VERIFY(from <= start->size() + truncSize());
            VERIFY(to <= start->size() + truncSize());
            if (from >= start->size()) {
                return truncSeq().Subseq(from - start->size(), to);
            } else {
                return getStart().getSeq().Subseq(from) + truncSeq().Subseq(0, to);
            }
        }

        Sequence suffix(size_t pos) const {
            VERIFY(pos <= truncSeq().size());
            size_t k = start->getSeq().size();
            if (pos >= k)
                return truncSeq().Subseq(pos - k, truncSeq().size());
            else {
                return getStart().getSeq().Subseq(pos) + truncSeq().Subseq(0, truncSeq().size());
            }
        }

        id_type getInnerId() const {
            return id;
        }

        std::string getShortId() const {
            return start->getShortId() + "." + itos(id.eid) + nuclLabel().str();
        }

        Sequence firstNucl() const {
            return truncSeq().Subseq(0, 1);
        }

        Sequence getSeq() const {
            VERIFY(!start->getSeq().empty());
            return start->getSeq() + seq;
        }

        size_t getStartSize() const {return start->size();}

        void DeleteEdge(BaseEdge &edge) {
            Locker<Vertex> locker({edge.start, edge.finish});
            DeleteEdgeLockFree(edge);
        }

        void DeleteEdgeLockFree(BaseEdge &edge) {
            Vertex &start = *edge.start;
            Vertex &rcstart = *edge.finish;
            BaseEdge &rcedge = edge.rc();
            if(edge != rcedge) {
                rcstart.innerRemoveEdge(rcedge);
            }
            start.innerRemoveEdge(edge);
        }

        size_t fullSize() const {
            return truncSeq().size() + getStartSize();
        }

        size_t innerSize() const {
            size_t full = fullSize();
            size_t vsize = getStart().size() + getFinish().size();
            return full >= vsize ? full - vsize : 0;
        }

        size_t overlapSize() const {
            size_t full = fullSize();
            size_t vsize = getStart().size() + getFinish().size();
            return full >= vsize ? 0 : vsize - full;
        }
    };

    template<class Traits>
    class BaseVertex {
    public:
        typedef typename Traits::VertexData VertexData;
        typedef typename Traits::EdgeData EdgeData;
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef ObjectId<Vertex, int> VertexId;
        typedef ObjectId<const Vertex, int> ConstVertexId;
        typedef int id_type;
    private:
        friend class AssemblyGraph<Traits>;
        friend class BaseEdge<Traits>;

        id_type id;
        Sequence seq = {};
        mutable std::list<Edge> outgoing_{};
        size_t _outDeg = 0;
        Vertex *rc_;
        omp_lock_t writelock = {};
        bool canonical;
        bool mark_ = false;
        int max_out_id = 0;

        Edge &innerAddEdge(Vertex &end, const Sequence &tseq, EdgeData data, BaseEdgeId eid = {});
        void setRC(Vertex &other) {
            rc_ = &other;
            rc_->rc_ = dynamic_cast<Vertex *>(this);
        }

    protected:
        virtual void fireAddEdge(Edge &edge) {}

    public:
        explicit BaseVertex(id_type id, Sequence seq) : id(id), seq(std::move(seq)), rc_(nullptr) {
            canonical = this->seq <= !this->seq;
            omp_init_lock(&writelock);
        }
        explicit BaseVertex(id_type id, bool canonical) : id(id), seq(), canonical(canonical), rc_(nullptr) {
            omp_init_lock(&writelock);
        }
        BaseVertex(const BaseVertex &) = delete;

        virtual ~BaseVertex()= default;

        id_type getInnerId() const {return id;};
        VertexId getId() {return {getInnerId(), dynamic_cast<Vertex *>(this)};}
        ConstVertexId getId() const {return {getInnerId(), dynamic_cast<const Vertex *>(this)};}
        std::string getShortId() const;

        void lock() { omp_set_lock(&writelock); }
        void unlock() { omp_unset_lock(&writelock); }

        bool isCanonical() const;

        void mark() { mark_ = true; }
        void unmark() { mark_ = false; }
        bool marked() const { return mark_; }

        Vertex &rc() { return *rc_; }
        const Vertex &rc() const { return *rc_; }

        void setSeq(Sequence _seq);

        size_t size() const { return seq.size(); }
        size_t getStartSize() const {return 0;};
        size_t truncSize() const {return seq.size();}
        int getMaxOutId() const {return max_out_id;}

        virtual Sequence getSeq() const { return seq; }
        Sequence truncSeq() const {return seq;}

        typename std::list<Edge>::iterator begin() const { return outgoing_.begin(); }
        typename std::list<Edge>::iterator end() const { return outgoing_.end(); }
        IterableStorage<TransformingIterator<typename std::list<Edge>::iterator, Edge>> incoming() {
            std::function<Edge &(Edge &)> transform = [](Edge &edge) -> Edge & {
                return edge.rc();
            };
            return {{rc().begin(), rc().end(), transform}, {rc().end(), rc().end(), transform}};
        }
        IterableStorage<TransformingIterator<typename std::list<Edge>::const_iterator, const Edge>> incoming() const {
            std::function<const Edge &(const Edge &)> transform = [](const Edge &edge) -> const Edge & {
                return edge.rc();
            };
            return {{rc().begin(), rc().end(), transform}, {rc().end(), rc().end(), transform}};
        }
        Edge &front() const { return outgoing_.front(); }
        Edge &back() const { return outgoing_.back(); }
        void sortOutgoing();

        Edge &getOutgoing(unsigned char c) const;
        bool hasOutgoing(unsigned char c) const;

        size_t outDeg() const { return _outDeg; }
        size_t inDeg() const { return rc_->outgoing_.size(); }
        bool isJunction() const;

        void checkConsistency() const;

//        This method should only be invoked if no graph modification is performed in parallel or if both start and
//        rc end vertices are locked by this process or otherwise prevented from modification by other processes
        Edge &addEdgeLockFree(Vertex &end, const Sequence &full_sequence, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {});
        Edge &addEdgeLockFree(Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {});
        bool removeEdgeLockFree(Edge &edge);

        Edge &addEdge(Vertex &end, const Sequence &full_seq, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {}) {
            Locker<BaseVertex<Traits>> locker({this, &end.rc()});
            return addEdgeLockFree(end, full_seq, std::move(data), eid, rcid);
        }

        Edge &addEdge(Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data = {}, BaseEdgeId eid = {}, BaseEdgeId rcid = {}) {
            Locker<BaseVertex<Traits>> locker({this, &end.rc()});
            return addEdgeLockFree(end, tseq, rctseq, std::move(data), eid, rcid);
        }

        void updateMaxOutId(int value) {max_out_id = std::max(max_out_id, value);}

        bool removeEdge(Edge &edge) {
            VERIFY(edge.getStart() == *this);
            Locker<BaseVertex<Traits>> locker({this, &edge.getFinish().rc()});
            return removeEdgeLockFree(edge);
        }

//        Use this method very carefully. It breaks basic graph assumptions but it is often important for parallelization.
//        This method removes one edge in a pair of rc edges. When a thread invokes this method no other thread should
//        try to access rc() method in the rc edge.
        bool innerRemoveEdge(Edge &edge) {
            VERIFY(edge.getStart() == *this);
            if(edge._rc != nullptr) {
                edge._rc->_rc = nullptr;
            }
            auto it = std::find(outgoing_.begin(), outgoing_.end(), edge);
            if(it == outgoing_.end())
                return false;
            outgoing_.erase(it);
            _outDeg--;
            return true;
        }

        //TODO: create method deleteEdge and do this properly
        void clear() {
            outgoing_.clear();
            _outDeg = 0;
        }

        bool operator==(const BaseVertex &other) const {return this == &other;}
        bool operator!=(const BaseVertex &other) const {return this != &other;}

        bool operator<(const BaseVertex &other) const {
            return (std::abs(id) << 1) + (id > 0) < (std::abs(other.id) << 1) + (other.id > 0);
        }

        bool operator<=(const BaseVertex &other) const {
            return *this < other || *this == other;
        }

        bool operator>(const BaseVertex &other) const {
            return (std::abs(id) << 1) + (id > 0) > (std::abs(other.id) << 1) + (other.id > 0);
        }

        bool operator>=(const BaseVertex &other) const {
            return *this > other || *this == other;
        }
    };

    template<class Traits>
    struct EdgePosition {
        typedef typename Traits::Edge Edge;
        Edge *edge;
        size_t pos;

        EdgePosition(Edge &_edge, size_t _pos) : edge(&_edge), pos(_pos) {VERIFY(pos >= 0 && pos <= edge->truncSize());}
        EdgePosition() : edge(nullptr), pos(0) {}

        Sequence kmerSeq() const {return edge->kmerSeq(pos);}
        unsigned char lastNucl() const {return edge->truncSeq()[pos - 1];}
        bool isBorder() const {return pos == 0 || pos == edge->truncSize();}

        std::vector<EdgePosition> step() const;
        EdgePosition RC() const {return {edge->rc(), edge->truncSize() - pos};}
    };


    template<class Traits>
    class AssemblyGraph {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Traits::VertexData VertexData;
        typedef typename Traits::EdgeData EdgeData;
        typedef std::list<Vertex> vertex_storage_type;
        typedef typename std::list<Vertex>::iterator vertex_iterator_type;
        typedef typename std::list<Vertex>::const_iterator const_vertex_iterator_type;
    private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
        vertex_storage_type vertex_list;
        int maxVId = 0;

//    Be careful since hash does not define vertex. Rc vertices share the same hash
        Vertex &innerAddVertex(typename Vertex::id_type id, bool canonical, VertexData data) {
            maxVId = std::max(std::abs(id), maxVId);
            vertex_list.emplace_back(id, canonical, std::move(data));
            return vertex_list.back();
        }

        Vertex &innerAddVertex(typename Vertex::id_type id, Sequence seq, VertexData data) {
            maxVId = std::max(std::abs(id), maxVId);
            vertex_list.emplace_back(id, std::move(seq), std::move(data));
            return vertex_list.back();
        }

    public:

        explicit AssemblyGraph() = default;
        AssemblyGraph(AssemblyGraph &&other) = default;
        AssemblyGraph &operator=(AssemblyGraph &&other) = default;
        AssemblyGraph(const AssemblyGraph &other) noexcept = delete;
        void fillFrom(const AssemblyGraph<Traits> &other);

        size_t size() const {return vertex_list.size();}
        size_t edgeCount() const;


        void removeIsolated();
        void removeMarked();
        void resetMarkers();

        Vertex &addVertexPair(const VertexData &data, typename Vertex::id_type id = 0);
        Vertex &addSelfRCVertex(VertexData data);
        Vertex &addVertex(const Sequence &seq, const VertexData &data = {}, typename Vertex::id_type id = 0);
        Vertex &addVertex(const Sequence &seq, typename Vertex::id_type id);
        Vertex &addVertex(const Vertex &other_graph_vertex);

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

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addSelfRCVertex(VertexData data) {
        typename Vertex::id_type id = maxVId + 1;
        Vertex &res = innerAddVertex(id, true, std::move(data));
        res.setRC(res);
        return res;
    }
    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertexPair(const VertexData &data, typename Vertex::id_type id) {
        VERIFY(id >= 0);
        if(id == 0)
            id = maxVId + 1;
        Vertex &rc = innerAddVertex(-id, false, data.RC());
        Vertex &res = innerAddVertex(id, true, std::move(data));
        res.setRC(rc);
        return res;
    }

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertex(const Vertex &other_graph_vertex) {
        typename AssemblyGraph<Traits>::Vertex &res = addVertex(other_graph_vertex.getSeq(), other_graph_vertex, other_graph_vertex.getInnerId());
        res.updateMaxOutId(other_graph_vertex.getMaxOutId());
        res.rc().updateMaxOutId(other_graph_vertex.rc().getMaxOutId());
        return res;
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type>> AssemblyGraph<Traits>::vertices(bool unique) & {
        std::function<bool(Vertex &)> use =
                [unique](Vertex &vertex) -> bool {
                    return !unique || vertex.isCanonical();
                };
        SkippingIterator<vertex_iterator_type> begin(vertex_list.begin(), vertex_list.end(), use);
        SkippingIterator<vertex_iterator_type> end(vertex_list.end(), vertex_list.end(), use);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type>> AssemblyGraph<Traits>::vertices(bool unique) const & {
        std::function<bool(const Vertex &)> use =
                [unique](const Vertex &vertex) -> bool {
                    return !unique || vertex.isCanonical();
                };
        SkippingIterator<const_vertex_iterator_type> begin(vertex_list.begin(), vertex_list.end(), use);
        SkippingIterator<const_vertex_iterator_type> end(vertex_list.end(), vertex_list.end(), use);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type>> AssemblyGraph<Traits>::verticesUnique() &{
        return vertices(true);
    }

    template<class Traits>
    IterableStorage<SkippingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type>> AssemblyGraph<Traits>::verticesUnique() const &{
        return vertices(true);
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type, typename Traits::Edge, 4>> AssemblyGraph<Traits>::edges(bool unique) & {
        std::function<std::array<Edge*, 4>(Vertex &)> apply = [unique](Vertex &vertex) {
            if(vertex.outDeg() > 4) {
                std::cout << vertex.getSeq() << std::endl;
                for(Edge &e : vertex) {
                    std::cout << e.truncSeq() << std::endl;
                }
            }
            VERIFY(vertex.outDeg() <= 4);
            std::array<Edge*, 4> res = {};
            size_t cur = 0;
            for(Edge &edge : vertex) {
                if(!unique || edge <= edge.rc()) {
                    res[cur] = &edge;
                    cur++;
                }
            }
            return res;
        };
        ApplyingIterator<vertex_iterator_type, Edge, 4> begin(vertex_list.begin(), vertex_list.end(), apply);
        ApplyingIterator<vertex_iterator_type, Edge, 4> end(vertex_list.end(), vertex_list.end(), apply);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type, const typename Traits::Edge, 4>> AssemblyGraph<Traits>::edges(bool unique) const & {
        std::function<std::array<const Edge*, 4>(const Vertex &)> apply = [unique](const Vertex &vertex) {
            std::array<const Edge*, 4> res = {};
            size_t cur = 0;
            for(const Edge &edge : vertex) {
                if(!unique || edge <= edge.rc()) {
                    res[cur] = &edge;
                    cur++;
                }
            }
            return res;
        };
        ApplyingIterator<const_vertex_iterator_type, const Edge, 4> begin(vertex_list.begin(), vertex_list.end(), apply);
        ApplyingIterator<const_vertex_iterator_type, const Edge, 4> end(vertex_list.end(), vertex_list.end(), apply);
        return {begin, end};
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::vertex_iterator_type, typename Traits::Edge, 4>> AssemblyGraph<Traits>::edgesUnique() &{
        return edges(true);
    }

    template<class Traits>
    IterableStorage<ApplyingIterator<typename AssemblyGraph<Traits>::const_vertex_iterator_type, const typename Traits::Edge, 4>> AssemblyGraph<Traits>::edgesUnique() const &{
        return edges(true);
    }

    template<class Traits>
    void AssemblyGraph<Traits>::removeIsolated() {
        vertex_storage_type newv;
        for (auto it = vertex_list.begin(); it != vertex_list.end();) {
            if (it->outDeg() == 0 && it->inDeg() == 0) {
                it = vertex_list.erase(it);
            } else {
                ++it;
            }
        }
    }

    template<class Traits>
    void AssemblyGraph<Traits>::removeMarked() {
        for (auto it = vertex_list.begin(); it != vertex_list.end();) {
            if (it->marked() || (it->inDeg() == 0 && it->outDeg() == 0)) {
                it = vertex_list.erase(it);
            } else {
                ++it;
            }
        }
    }

    template<class Traits>
    void AssemblyGraph<Traits>::resetMarkers() {
        for(Vertex &vertex : vertices()) {
            for(Edge &edge : vertex) {
                edge.mark(EdgeMarker::common);
            }
            vertex.unmark();
        }
    }

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertex(const Sequence &seq, const VertexData &data, typename Vertex::id_type id) {
        if (id == 0) {
            if(seq <= !seq)
                id = maxVId + 1;
            else
                id = -maxVId - 1;
        }
        maxVId = std::max(std::abs(id), maxVId);
        Vertex *res = nullptr;
        Vertex *rc = nullptr;
        if(seq == !seq) {
            res = &innerAddVertex(id, seq, data);
            rc = res;
        } else {
            res = &innerAddVertex(id, seq, data);;
            rc = &innerAddVertex(-id, !seq, data.RC());
        }
        res->setRC(*rc);
        res->setSeq(seq);
        return *res;
    }

    template<class Traits>
    typename AssemblyGraph<Traits>::Vertex &AssemblyGraph<Traits>::addVertex(const Sequence &seq, typename Vertex::id_type id) {
        return addVertex(seq, VertexData(), id);
    }

    template<class Traits>
    size_t AssemblyGraph<Traits>::edgeCount() const {
        size_t res = 0;
        for(auto &v : vertices())
            res += v.outDeg();
        return res;
    }

    template<class Traits>
    void AssemblyGraph<Traits>::fillFrom(const AssemblyGraph<Traits> &other) {
        for(const Vertex &v : other.verticesUnique()) {
            addVertex(v);
        }
        IdIndex<Vertex> index(vertices().begin(), vertices().end());
        for(const Edge &e : other.edgesUnique()) {
            index.getById(e.getStart().getInnerId()).addEdgeLockFree(index.getById(e.getFinish().getInnerId()), e.getSeq(), e, e.getInnerId(), e.rc().getInnerId());
        }
        for(Edge &edge : edges()) {
            edge.incCov(-edge.intCov());
        }
    }

    template<class Traits>
    std::vector<EdgePosition<Traits>> EdgePosition<Traits>::step() const {
        if (pos == edge->truncSize()) {
            std::vector<EdgePosition> res;
            auto &v = edge->getFinish();
            for (Edge &next : v) {
                res.emplace_back(next, 1);
            }
            return std::move(res);
        } else {
            return {{*edge, pos + 1}};
        }
    }


//std::ostream &operator<<(std::ostream &os, const Edge &edge) {
//    os << edge.getShortId();
//    return os;
//}

    template<class Traits>
    bool BaseVertex<Traits>::isCanonical() const {
        return canonical;
    }

    template<class Traits>
    void BaseVertex<Traits>::checkConsistency() const {
        VERIFY(isCanonical() || rc().isCanonical());
        VERIFY(*this == rc() || isCanonical() != rc().isCanonical());
        VERIFY(isCanonical() == seq <= !seq);
        for (const Edge &edge : outgoing_) {
            VERIFY(edge.isCanonical() || edge.rc().isCanonical());
//            VERIFY(edge.isCanonical() == edge.getSeq() <= !edge.getSeq());
            VERIFY(edge == edge.rc() || (edge.isCanonical() != edge.rc().isCanonical()));
            VERIFY(edge.rc().getFinish() == this->rc());
            VERIFY(std::find(edge.getFinish().rc().begin(), edge.getFinish().rc().end(), edge.rc()) != edge.getFinish().rc().end());
            VERIFY(edge.intCov() == edge.rc().intCov());
        }
    }

    template<class Traits>
    std::string BaseVertex<Traits>::getShortId() const {
        std::stringstream ss;
        ss << getInnerId() % 1000000000;
        return ss.str();
    }

    template<class Traits>
    void BaseVertex<Traits>::setSeq(Sequence _seq) {
        lock();
        if (getSeq().empty()) {
            seq = std::move(_seq);
            unlock();
            if(rc_ != nullptr) {
                rc_->lock();
                rc_->seq = !getSeq();
                rc_->unlock();
            }
        } else {
            unlock();
        }
    }

//}

//void Vertex::clearSequence() {
//    if (!seq.empty()) {
//        seq = Sequence();
//        rc_->seq = Sequence();
//    }

    template<class Traits>
    typename Traits::Edge &BaseVertex<Traits>::innerAddEdge(Vertex &end, const Sequence &tseq, EdgeData data, BaseEdgeId eid) {
        if (!eid.valid()) {
            eid = {id, max_out_id + 1};
        }
        updateMaxOutId(eid.eid);
        outgoing_.emplace_back(eid, *dynamic_cast<Vertex*>(this), end, tseq, std::move(data));
        Edge &edge = outgoing_.back();
        _outDeg++;
        fireAddEdge(edge);
        return edge;
    }


//Edge &Vertex::addEdgeLockFree(const Edge &edge) {
//    VERIFY(this == edge.getStart());
//    for (Edge &e : outgoing_) {
//        if (edge.size() <= e.size()) {
//            if (edge.truncSeq() == e.truncSeq().Subseq(0, edge.size())) {
//                return e;
//            }
//        } else if (edge.truncSeq().Subseq(0, e.size()) == e.truncSeq()) {
//            e = edge;
//            return e;
//        }
//    }
//    outgoing_.emplace_back(edge);
//    return outgoing_.back();
//}
//
//Edge &Vertex::addEdge(const Edge &e) {
//    omp_set_lock(&writelock);
//    Edge &res = addEdgeLockFree(e);
//    omp_unset_lock(&writelock);
//    return res;
//}

    template<class Traits>
    typename Traits::Edge &BaseVertex<Traits>::getOutgoing(unsigned char c) const {
        for (Edge &edge : outgoing_) {
            if (edge.truncSeq()[0] == c) {
                return edge;
            }
        }
        std::cout << getSeq() << std::endl;
        std::cout << size_t(c) << std::endl;
        for (const Edge &edge : outgoing_) {
            std::cout << edge.truncSeq() << std::endl;
        }
        VERIFY(false);
        return outgoing_.front();
    }

    template<class Traits>
    bool BaseVertex<Traits>::hasOutgoing(unsigned char c) const {
        for (const Edge &edge : outgoing_) {
            if (edge.truncSeq()[0] == c) {
                return true;
            }
        }
        return false;
    }

    template<class Traits>
    void BaseVertex<Traits>::sortOutgoing() {
        outgoing_.sort();
//    std::sort(outgoing_.begin(), outgoing_.end());
    }

    template<class Traits>
    bool BaseVertex<Traits>::isJunction() const {
        return outDeg() != 1 || inDeg() != 1;
    }

    template<class Traits>
    typename Traits::Edge &BaseVertex<Traits>::addEdgeLockFree(Vertex &end, const Sequence &full_sequence, EdgeData data, BaseEdgeId eid, BaseEdgeId rcid) {
        for(Edge &edge: outgoing_) {
            if(edge.getFinish() == end && edge.truncSeq() == full_sequence.Subseq(size())) {
                return edge;
            }
        }
        Edge &res = innerAddEdge(end, full_sequence.Subseq(size()), data, eid);
        if(end.rc() != *this || full_sequence != !full_sequence) {
            Edge &rc_edge = end.rc().innerAddEdge(rc(), full_sequence.rc().Subseq(end.size()), data.RC(), rcid);
            res._rc = &rc_edge;
            rc_edge._rc = &res;
        } else {
            res._rc = &res;
        }
        VERIFY(res.fullSize() == res.rc().fullSize());
        return res;
    }

    template<class Traits>
    typename Traits::Edge &BaseVertex<Traits>::addEdgeLockFree(Vertex &end, const Sequence &tseq, const Sequence &rctseq, EdgeData data, BaseEdgeId eid, BaseEdgeId rcid) {
        for(Edge &edge: outgoing_) {
            if(edge.getFinish() == end && edge.truncSeq() == tseq) {
                return edge;
            }
        }
        Edge &res = innerAddEdge(end, tseq, data, eid);
        if(end.rc() != *this || (tseq.size() > size() && tseq.rc().Subseq(size()) != rctseq.rc().Subseq(size()))) {
            Edge &rc_edge = end.rc().innerAddEdge(rc(),rctseq, data.RC(), rcid);
            res._rc = &rc_edge;
            rc_edge._rc = &res;
        } else {
            res._rc = &res;
        }
        VERIFY(res.fullSize() == res.rc().fullSize());
        return res;
    }

    template<class Traits>
    bool BaseVertex<Traits>::removeEdgeLockFree(Edge &edge) {
        if(edge != edge.rc()) {
            edge.rc().getStart().innerRemoveEdge(edge.rc());
        }
        return innerRemoveEdge(edge);
    }

}
