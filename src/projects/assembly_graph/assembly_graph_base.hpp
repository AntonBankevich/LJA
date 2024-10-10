#pragma once
#include "sequences/sequence.hpp"
#include "common/iterator_utils.hpp"
#include "common/object_id.hpp"
#include "common/id_index.hpp"
#include "sequences/contigs.hpp"
#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <array>
#include <stdexcept>

namespace ag {

//    TODO: Make all Traits-templated complex classes virtual descendants of Traits. This way all typedefs will be there automatically.
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

    struct EdgeSaveLabel {
        ag::BaseEdgeId fId;
        ag::BaseEdgeId rcId;
        EdgeSaveLabel(ag::BaseEdgeId fId, ag::BaseEdgeId rcId) : fId(fId), rcId(rcId) {
        }
    };

    inline std::ostream &operator<<(std::ostream &os, ag::EdgeSaveLabel val) {
        return os << val.fId << "_" << val.rcId;
    }
}

template<>
inline ag::BaseEdgeId Parse<ag::BaseEdgeId>(const std::string &s, size_t start, size_t end) {
    size_t pos = s.find('.', start);
    if(pos == size_t(-1) || pos >= end){ throw std::invalid_argument("Incorrect edge id record"); }
    int vid = Parse<int>(s, start, pos);
    int eid = Parse<int>(s, pos + 1, end);
    return {vid, eid};
}

template<>
inline ag::EdgeSaveLabel Parse<ag::EdgeSaveLabel>(const std::string &s, size_t start, size_t end) {
    size_t pos = s.find('_', start);
    if (pos == size_t(-1) || pos >= end) { throw std::invalid_argument("Incorrect edge id pair record"); }
    return {Parse<ag::BaseEdgeId>(s, start, pos), Parse<ag::BaseEdgeId>(s, pos + 1, end)};
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
        std::vector<T> locked;
    public:
        explicit Locker(std::vector<T> _vertices) : locked(std::move(_vertices)) {
            std::sort(locked.begin(), locked.end());
            locked.erase(std::unique(locked.begin(), locked.end()), locked.end());
            for(T &v: locked) {
                v->lock();
            }
        }
        Locker(const Locker &other) = delete;
        Locker(Locker &&other) = delete;
        Locker& operator=(const Locker &other) = delete;
        Locker& operator=(Locker &&other) = delete;

        ~Locker() {
            for(T &v: locked) {
                v->unlock();
            }
        }
    };

    template<class T>
    class RCIterator {
    private:
        T *val;
        bool rc;
        bool isend;
    public:
        RCIterator(T &obj, bool rc, bool isend) : val(&obj), rc(rc), isend(isend) {
        }
        static RCIterator begin(T &obj) {return {obj, false, false};}
        static RCIterator end(T &obj) {return {obj, false, true};}
        bool operator==(const RCIterator &other) const {return val == other.val && rc == other.rc && isend == other.isend;}
        bool operator!=(const RCIterator &other) const {return !(*this == other);}
        RCIterator &operator++() {
            VERIFY(!isend);
            if(rc || *val == val->rc()) {
                rc = false;
                isend = true;
            } else {
                rc = true;
            }
            return *this;
        }

        RCIterator operator++(int) const {
            RCIterator res = *this;
            ++res;
            return res;
        }
        T& operator*() const {
            VERIFY(!isend);
            if(rc)
                return val->rc();
            else
                return *val;
        }
    };

    template<class T>
    IterableStorage<RCIterator<T>> ThisAndRC(T&obj) {
        return {RCIterator<T>::begin(obj), RCIterator<T>::end(obj)};
    }

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
//    TODO: this is strange to have declaration of AssemblyGraph here. Need to make sure it does not break anything.
    template<class Traits>
    class AssemblyGraph;



    template<class Traits>
    class BaseEdge {
        friend class BaseVertex<Traits>;
        friend class AssemblyGraph<Traits>;
    public:
        typedef BaseEdgeId id_type;
        typedef typename Traits::VertexData VertexData;
        typedef typename Traits::EdgeData EdgeData;
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef ObjectId<Edge, id_type> EdgeId;
        typedef ConstObjectId<Edge, id_type> ConstEdgeId;
    private:
        id_type id;
        Vertex *start;
        Vertex *finish;
        Sequence seq;
        Edge *_rc;
        ag::EdgeMarker marker = ag::EdgeMarker::common;

    public:
        bool fire_create = false;
        bool fire_destroy = false;


        BaseEdge(id_type id, Vertex &_start, Vertex &_end, Sequence _seq) :
                id(std::move(id)), start(&_start), finish(&_end), seq(std::move(_seq)), _rc(nullptr) {
        }
        BaseEdge() : id(0, 0), start(nullptr), finish(nullptr), seq(), _rc(nullptr) {
//            TODO: Remove this!!! It exists only for Andreys code compilation but that code should be purged
        }
//        BaseEdge() : id(), start(nullptr), finish(nullptr), seq() {}
        BaseEdge(BaseEdge &&other) = delete;
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

        bool isOuter() const { return getStart().outDeg() > 1 && getFinish().inDeg() > 1; }
        bool isInner() const { return getStart().outDeg() == 1 && getFinish().inDeg() == 1; }

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

        Sequence firstNucl() const {
            return truncSeq().Subseq(0, 1);
        }

        Sequence getSeq() const {
            VERIFY(!start->getSeq().empty());
            return start->getSeq() + seq;
        }

        size_t getStartSize() const {return start->size();}

        void DeleteEdge(BaseEdge &edge) {
            Locker<typename Vertex::VertexId> locker({edge.start, edge.finish});
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

    template<class T>
    inline std::string DefaultEdgeName(BaseEdge<T> &edge) {
        return edge.getInnerId().str();
    }

    template<class T>
    inline std::string SaveEdgeName(BaseEdge<T> &edge) {
        VERIFY((edge.getFinish().rc().getInnerId() > 0) == edge.getFinish().rc().isCanonical());
        return edge.getInnerId().str() + "_" + edge.rc().getInnerId().str();
    }


    template<class Traits>
    class BaseVertex {
        friend class AssemblyGraph<Traits>;
        friend class BaseEdge<Traits>;
    public:
        typedef typename Traits::VertexData VertexData;
        typedef typename Traits::EdgeData EdgeData;
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef ObjectId<Vertex, int> VertexId;
        typedef ConstObjectId<Vertex, int> ConstVertexId;
        typedef int id_type;
    private:
        id_type id;
        Sequence seq = {};
        mutable std::list<Edge> outgoing_{};
        size_t _outDeg = 0;
        Vertex *rc_;
        omp_lock_t writelock = {};
        bool canonical;
        bool mark_ = false;
        std::array<int, 5> max_out_id = {0,0,0,0, 0};

        Edge &innerAddEdge(Vertex &end, const Sequence &tseq, EdgeData data, BaseEdgeId eid = {});
        void setRC(Vertex &other) {
            rc_ = &other;
            rc_->rc_ = dynamic_cast<Vertex *>(this);
        }
        //        This method should only be invoked if no graph modification is performed in parallel or if both start and
//        rc end vertices are locked by this process or otherwise prevented from modification by other processes
//        void removeEdgeLockFree(Edge &edge);

    protected:
        virtual void fireAddEdge(Edge &edge) {}

    public:
        bool fire_create = false;
        bool fire_destroy = false;

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

        void lock() { omp_set_lock(&writelock); }
        void unlock() { omp_unset_lock(&writelock); }

        bool isCanonical() const;
        bool isPalindrome() const;

        void mark() { mark_ = true; }
        void unmark() { mark_ = false; }
        bool marked() const { return mark_; }

        Vertex &rc() { return *rc_; }
        const Vertex &rc() const { return *rc_; }

        void setSeq(Sequence _seq);

        size_t size() const { return seq.size(); }
        size_t getStartSize() const {return 0;};
        size_t truncSize() const {return seq.size();}
        std::array<int, 5> getMaxOutId() const {return max_out_id;}

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
        bool hasOutgoingSuffix() const;

        size_t outDeg() const { return _outDeg; }
        size_t inDeg() const { return rc_->outgoing_.size(); }
        bool isJunction() const;

        void checkConsistency() const;

        void updateMaxOutId(int value) {
            max_out_id[value % 10] = std::max(max_out_id[value%10], value / 10);
        }

        void updateMaxOutId(const std::array<int, 5> other) {
            for(size_t i = 0; i < max_out_id.size(); i++)
                updateMaxOutId(i + other[i] * 10);
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
    std::ostream &operator<<(std::ostream &os, const BaseVertex<Traits> &vertex) {
        if(vertex.size() > 8)
            os << vertex.getId();
        else
            os << vertex.getSeq();
        return os;
    }

    template<class Traits>
    std::ostream &operator<<(std::ostream &os, const BaseEdge<Traits> &edge) {
        os << edge.getStart() << "->" << edge.getFinish() << "(" << edge.firstNucl() << edge.truncSize() << ")";
        return os;
    }


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
    void BaseVertex<Traits>::setSeq(Sequence _seq) {
        lock();
        VERIFY(isCanonical() == _seq.isCanonical());
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
            int code = tseq.empty() ? 4 : tseq[0];
            eid = {id, (max_out_id[code] + 1) * 10 + code};
        }
        VERIFY(eid.vid == id);
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
        size_t cnt = 0;
        for (Edge &edge : outgoing_) {
            if (edge.nuclLabel().empty() || edge.nuclLabel()[0] == c) {
                cnt++;
            }
        }
        VERIFY(cnt <= 1);
        for (Edge &edge : outgoing_) {
            if (edge.nuclLabel().empty() || edge.nuclLabel()[0] == c) {
                return edge;
            }
        }
        std::cout << c << std::endl;
        std::cout << getSeq() << std::endl;
        for (const Edge &edge : outgoing_) {
            std::cout << edge.nuclLabel() << std::endl;
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

//    template<class Traits>
//    void BaseVertex<Traits>::removeEdgeLockFree(Edge &edge) {
//        if(edge != edge.rc()) {
//            edge.rc().getStart().innerRemoveEdge(edge.rc());
//        }
//        innerRemoveEdge(edge);
//    }

    template<class Traits>
    bool BaseVertex<Traits>::hasOutgoingSuffix() const {
        for(Edge &edge : *this) {
            if(edge.truncSize() == 0)
                return true;
        }
        return false;
    }

    template<class Traits>
    bool BaseVertex<Traits>::isPalindrome() const {return *this == rc();}
}