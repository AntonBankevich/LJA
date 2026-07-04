#pragma once
#include "sequences/sequence.hpp"
#include "common/iterator_utils.hpp"
#include "common/object_id.hpp"
#include "common/omp_utils.hpp"
#include "common/id_index.hpp"
#include "sequences/contigs.hpp"
#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <common/hash_utils.hpp>

namespace ag {

//    Edge id currently is a pair: id of the start vertex (which is int) and another int number.
//    The last decimal digit of edge id is its first nucleotide.
    class EdgeIdType {
    public:
        int vid;
        int eid;
        EdgeIdType() = default;
        EdgeIdType(int vid, int eid) : vid(vid), eid(eid) {
        }

        bool valid() const {
            VERIFY((vid != 0) == (eid != 0));
            return vid != 0;
        }

        bool operator==(const EdgeIdType &other) const {return vid == other.vid && eid == other.eid;}
        bool operator!=(const EdgeIdType &other) const {return vid != other.vid || eid != other.eid;}
        bool operator<(const EdgeIdType &other) const {return vid < other.vid || (vid == other.vid && eid < other.eid);}
        bool operator>(const EdgeIdType &other) const {return vid > other.vid || (vid == other.vid && eid > other.eid);}
        bool operator<=(const EdgeIdType &other) const {return vid < other.vid || (vid == other.vid && eid <= other.eid);}
        bool operator>=(const EdgeIdType &other) const {return vid > other.vid || (vid == other.vid && eid >= other.eid);}

        std::string str() const {
            return itos(vid) + "." +itos(eid);
        }
    };

    inline std::ostream &operator<<(std::ostream &os, ag::EdgeIdType val) {
        return os << val.vid << "." << val.eid;
    }

    struct EdgeSaveLabel {
        ag::EdgeIdType fId;
        ag::EdgeIdType rcId;
        EdgeSaveLabel(ag::EdgeIdType fId, ag::EdgeIdType rcId) : fId(fId), rcId(rcId) {
        }
    };

}
namespace std {
    inline std::ostream &operator<<(std::ostream &os, ag::EdgeSaveLabel val);
}

template<>
inline ag::EdgeIdType Parse<ag::EdgeIdType>(const std::string &s, size_t start, size_t end) {
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
    return {Parse<ag::EdgeIdType>(s, start, pos), Parse<ag::EdgeIdType>(s, pos + 1, end)};
}

namespace std {
//    For some reason compiler finds this hash function only when you put it in std even though it is a bad practice.
    template <>
    class hash < ag::EdgeIdType >{
    public:
        size_t operator()(const ag::EdgeIdType &x ) const;
    };
}

namespace ag {

    template<class T>
    class RCIterator {
    private:
        T *val;
        bool rc;
        bool isend;
    public:
        typedef T value_type;
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

    class Edge;
    typedef ObjectId<Edge, EdgeIdType> EdgeId;
    typedef ConstObjectId<Edge, EdgeIdType> ConstEdgeId;

    struct EdgeData {
    protected:
        Sequence edge_code = {};
        ag::EdgeMarker marker = ag::EdgeMarker::common;
        bool corporeal = true;
    public:
//        TODO: get rid of this or at least control access
        mutable bool is_reliable = false;
        std::vector<EdgeId> label = {};
        size_t read_tail_count = 0;
        size_t read_tail_length = 0;
    protected:

//        dbg-specific fields
        size_t cov = 0;
    public:
        EdgeData() = default;
        EdgeData RC() const {return {};}
    };
    struct VertexData {
        static const hashing::htype default_hash;
    protected:
//        dbg-specific fields
        hashing::htype hash = default_hash;
        bool cyclic = false;
        bool inf_left = false;
        bool inf_right = false;
    public:
        size_t subread_length = 0;
        size_t subread_count = 0;
        size_t covering_read_count = 0;
        VertexData() = default;
        VertexData RC() const {return {*this};}
        static VertexData SPGData(bool cyclic, bool inf_left, bool inf_right);
        static VertexData DBGData(hashing::htype hash = default_hash);
    };

    class Vertex;
    class Edge;
//    TODO: this is strange to have declaration of AssemblyGraph here. Need to make sure it does not break anything.
    class AssemblyGraph;
    class EdgeCodeListener;

//    Hidden contracts:
//    - Every edge has a reverse-complement counterpart rc(): if this edge goes u->v then rc() goes
//      v.rc()->u.rc(). The whole graph is rc-symmetric, and every graph-editing operation (and every
//      listener callback it fires) is applied to both an edge and its rc() as a matched pair. See
//      graph_listeners.hpp's Fire dispatchers for the exact mirroring pattern.
//    - isPrefix()/isSuffix() are normal Supregraph structure, not degenerate/corner cases: a DBG never
//      has prefix/suffix edges, while in a Supregraph every edge is one or the other (adjacent vertices'
//      sequences are in a prefix/suffix relationship rather than a k-mer overlap).
    class Edge : public EdgeData {
        friend class Vertex;
        friend class AssemblyGraph;
        friend class EdgeCodeListener;
    public:
        typedef EdgeIdType id_type;
        typedef EdgeId pointer_type;
        typedef ConstEdgeId const_pointer_type;
        typedef std::false_type is_pointer;
    private:
        id_type id;
        Vertex *start;
        Vertex *finish;
        Sequence seq;
        Edge *_rc;

    public:
        Edge(id_type id, Vertex &_start, Vertex &_end, Sequence _seq, EdgeData data);
        Edge();
        Edge(Edge &&other) = delete;
        Edge(const Edge &other) = delete;
        virtual ~Edge() = default;

        bool isCanonical() const {return *this <= rc();}
        bool isOuter() const;
        bool isInner() const;
        bool isPrefix() const { return rc().truncSize() == 0; }
        bool isSuffix() const { return truncSize() == 0; }
        const Edge &getCanonical() const {return isCanonical() ? *this : rc();}
        Edge &getCanonical() {return isCanonical() ? *this : rc();}

        const Sequence &truncSeq() const { return seq; }
        const Sequence &getCode() const {return edge_code;}
        Sequence kmerSeq(size_t pos) const {return fullSubseq(pos, pos + getStartSize());}
        Sequence getSeq() const;// = fullSeq
        Sequence fullSeq() const;
        Sequence fullSubseq(size_t from, size_t to) const;
        Sequence suffix(size_t pos) const;
        unsigned char firstNucl() const {return seq[0];}

        size_t truncSize() const {return truncSeq().size();}
        size_t fullSize() const {return truncSeq().size() + getStartSize();}
        size_t innerSize() const;
        size_t overlapSize() const;

        EdgeId getId() {return {getInnerId(), this};}
        ConstEdgeId getId() const {return {getInnerId(), this};}
        const Vertex &getFinish() const {return *finish;}
        Vertex &getFinish() {return *finish;}
        const Vertex &getStart() const {return *start;}
        Vertex &getStart() {return *start;}
        id_type getInnerId() const {return id;}
        size_t getStartSize() const;
        ag::EdgeMarker getMarker() const { return marker; };
        Edge &rc() {return *_rc;}
        const Edge &rc() const {return *_rc;}
        std::string str() const {return getInnerId().str();}

        bool operator==(const Edge &other) const {return this == &other;}
        bool operator!=(const Edge &other) const {return this != &other;}
        bool operator<(const Edge &other) const;
        bool operator>(const Edge &other) const;
        bool operator<=(const Edge &other) const {return *this == other || *this < other;}

        void mark(ag::EdgeMarker _marker) { marker = _marker; };
//        Non-corporeal edges can not be detected by getOutgoing(c), but still will be among outgoing/incoming edges
        void setCorporeal(bool value) {corporeal = value;}

        void DeleteEdge(Edge &edge);
        void DeleteEdgeLockFree(Edge &edge);

//        dbg-specific methods
        void incCov(int64_t delta);
        size_t intCov() const {return cov;}
        void setCov(size_t val) {cov = val;}
        double getCoverage() const {return double(cov) / truncSize();}
    };

    inline std::string DefaultEdgeName(const Edge &edge) {return edge.getInnerId().str();}

    // inline std::string SaveEdgeName(const Edge &edge) {return edge.getInnerId().str() + "_" + edge.rc().getInnerId().str();}

    typedef ObjectId<Vertex, int> VertexId;
    typedef ConstObjectId<Vertex, int> ConstVertexId;

//    Like Edge, every vertex has an rc() counterpart and the graph is edited/listened-to symmetrically
//    (see the contract note above class Edge).
    class Vertex : public VertexData {
        friend class AssemblyGraph;
        friend class Edge;
    public:
        typedef int id_type;
        typedef VertexId pointer_type;
        typedef ConstVertexId const_pointer_type;

        typedef std::false_type is_pointer;
    private:
        id_type id;
        Sequence seq = {};
        mutable std::list<Edge> outgoing_{};
        size_t _outDeg = 0;
        Vertex *rc_;
        omp_lock_t writelock = {};
        bool canonical;
        bool mark_ = false;
    private:
        std::array<int, 5> max_out_id = {0,0,0,0, 0};

        Edge &innerAddEdge(Vertex &end, const Sequence &tseq, EdgeData data, EdgeIdType eid = {});
        bool innerRemoveEdge(Edge &edge);
        void setRC(Vertex &other);

    public:
        explicit Vertex(id_type id, Sequence seq, VertexData data);
        explicit Vertex(id_type id, bool canonical, VertexData data);
        Vertex(const Vertex &) = delete;
        Vertex &operator=(const Vertex &) = delete;
        virtual ~Vertex()= default;

//        Sequence methods
        virtual Sequence getSeq() const { return seq; }
        Sequence truncSeq() const {return seq;}

//        Size methods
        size_t size() const { return seq.size(); }
        size_t getStartSize() const {return 0;};
        size_t truncSize() const {return seq.size();}

//        Incident edges
        typename std::list<Edge>::iterator begin() const { return outgoing_.begin(); }
        typename std::list<Edge>::iterator end() const { return outgoing_.end(); }
        IterableStorage<TransformingIterator<typename std::list<Edge>::iterator, Edge>> incoming();
        IterableStorage<TransformingIterator<typename std::list<Edge>::const_iterator, const Edge>> incoming() const;
        Edge &front() const { return outgoing_.front(); }
        Edge &back() const { return outgoing_.back(); }
        Edge &getOutgoing(unsigned char c) const;
//        This method works even when iterator points to end when next edge is a suffix edge
        template<class I>
        Edge &getOutgoingByIterator(I iterator) const;
        bool hasOutgoing(unsigned char c) const;
        bool hasOutgoingSuffix() const;
        size_t outDeg() const { return _outDeg; }
        size_t inDeg() const { return rc_->outgoing_.size(); }
        const Vertex &getCanonical() const {return isCanonical() ? *this : rc();}
        Vertex &getCanonical() {return isCanonical() ? *this : rc();}

//        General information
        id_type getInnerId() const {return id;};
        VertexId getId() {return {getInnerId(), this};}
        ConstVertexId getId() const {return {getInnerId(), this};}
        bool marked() const { return mark_; }
        Vertex &rc() { return *rc_; }
        const Vertex &rc() const { return *rc_; }
        std::array<int, 5> getMaxOutId() const {return max_out_id;}
        bool isCyclic() const { return cyclic; }
        bool isInfLeft() const { return inf_left; }
        bool isInfRight() const { return inf_right; }
        bool isCore() const;
        bool isOuter() const;


        //        dbg-specific methods
        hashing::htype getHash() const {return hash;}

//        Property checking
        bool isCanonical() const {return canonical;}
        bool isPalindrome() const;
        bool isJunction() const;
        void checkConsistency() const;

//        Functions that change contents of the vertex
        void lock() { omp_set_lock(&writelock); }
        void unlock() { omp_unset_lock(&writelock); }
        void mark() { mark_ = true; }
        void unmark() { mark_ = false; }
        void setSeq(Sequence _seq);
        void sortOutgoing();
        void updateMaxOutId(int value) {max_out_id[value % 10] = std::max(max_out_id[value%10], value / 10);}
        void updateMaxOutId(const std::array<int, 5> other);
        //TODO: create method deleteEdge and do this properly
        void clear();

//        Comparison functions
        bool operator==(const Vertex &other) const {return this == &other;}
        bool operator!=(const Vertex &other) const {return this != &other;}
        bool operator<(const Vertex &other) const;
        bool operator<=(const Vertex &other) const;
        bool operator>(const Vertex &other) const;
        bool operator>=(const Vertex &other) const;
    };

    std::ostream &operator<<(std::ostream &os, const Vertex &vertex);
    std::ostream &operator<<(std::ostream &os, const Edge &edge);


    struct EdgePosition {
        EdgeId edge;
        size_t pos;

        EdgePosition(Edge &_edge, size_t _pos) : edge(_edge.getId()), pos(_pos) {VERIFY(pos >= 0 && pos <= edge->truncSize());}
        EdgePosition() : edge(), pos(0) {}

        Sequence kmerSeq() const {return edge->kmerSeq(pos);}
        unsigned char lastNucl() const {return edge->truncSeq()[pos - 1];}
        bool isBorder() const {return pos == 0 || pos == edge->truncSize();}

        std::vector<EdgePosition> step() const;
        EdgePosition RC() const {return {edge->rc(), edge->truncSize() - pos};}
        bool operator==(const EdgePosition &other) const {return edge == other.edge && pos == other.pos;}
        bool operator!=(const EdgePosition &other) const {return !(*this == other);}
        bool operator<(const EdgePosition &other) const {return *edge < *other.edge || (edge == other.edge && pos < other.pos); }
        bool operator>(const EdgePosition &other) const {return *edge > *other.edge || (edge == other.edge && pos > other.pos); }
    };

    std::ostream &operator<<(std::ostream &os, const EdgePosition &epos);


    template<class I>
    Edge &Vertex::getOutgoingByIterator(I iterator) const {
        size_t cnt = 0;
        for (Edge &edge : outgoing_) {
            if (edge.corporeal && (edge.isSuffix() || edge.getCode()[0] == *iterator)) {
                cnt++;
            }
        }
        VERIFY(cnt <= 1);
        for (Edge &edge : outgoing_) {
            if (edge.corporeal && (edge.isSuffix() || edge.getCode()[0] == *iterator)) {
                return edge;
            }
        }
	    std::cout << *this << std::endl;
        std::cout << *iterator << std::endl;
        std::cout << getSeq() << std::endl;
        for (const Edge &edge : outgoing_) {
            std::cout << "Outgoing code: " << edge.getCode() << std::endl;
        }
        VERIFY(false);
        return outgoing_.front();
    }

}
