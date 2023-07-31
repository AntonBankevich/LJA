//
// Created by anton on 7/22/20.
//

#pragma once
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include <common/oneline_utils.hpp>
#include <common/iterator_utils.hpp>
#include <utility>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <forward_list>
#include <list>

namespace dbg {
    enum EdgeMarker {
        incorrect,
        suspicious,
        common,
        possible_break,
        correct,
        unique,
        repeat
    };

    bool IsMarkerCorrect(EdgeMarker marker);

    class SparseDBG;
    class Vertex;
    class Edge;

//    TODO: this class should be constructed as a mixture of multiple classes each representing possible piece of informations
//that is to be stored in the edge. Corresponding information should be able to support itself during various graph
//operations.
    class EdgeData {
    private:
        Edge *edge;
        size_t cov = 0;
        EdgeMarker marker = EdgeMarker::common;
    public:
        explicit EdgeData(Edge &edge) : edge(&edge) {}
        bool is_reliable = false;
        void incCov(int delta) {cov += delta;}
        size_t intCov() const {return cov;}
        double getCoverage() const;
        void mark(EdgeMarker _marker) { marker = _marker; };
        EdgeMarker getMarker() const { return marker; };
    };

    class Edge {
    private:
        Vertex *start;
        Vertex *finish;
        Sequence seq;
        Edge *_rc;
        mutable EdgeData data;

        friend class dbg::Vertex;

    public:
        typedef std::string id_type;
        friend class Vertex;


        Edge(Vertex &_start, Vertex &_end, Sequence _seq) :
                start(&_start), finish(&_end), seq(std::move(_seq)), _rc(nullptr), data(*this) {
        }
        Edge() : start(nullptr), finish(nullptr), seq(), data(*this) {}

        bool isCanonical() const {return *this <= rc();}

        EdgeData &getData() const {return data;}

        id_type getInnerId() const;
        std::string oldId() const;
        std::string getShortId() const;
        std::string str() const;

        Sequence getSeq() const;
        const Sequence &truncSeq() const { return seq; }
        Sequence firstNucl() const;
        Sequence kmerSeq(size_t pos) const {return fullSubseq(pos, pos);}
        Sequence fullSubseq(size_t from, size_t to) const;
        Sequence suffix(size_t pos) const;
        size_t getStartSize() const;
        size_t truncSize() const;

        const Vertex &getFinish() const {return *finish;}
        Vertex &getFinish() {return *finish;}
        const Vertex &getStart() const {return *start;}
        Vertex &getStart() {return *start;}
        Edge &rc() const;

        //        This method should only be invoked if no graph modification is performed in parallel or if both start and
//        rc end vertices are locked by this process or otherwise prevented from modification by other processes
        static void DeleteEdgeLockFree(Edge &edge);
        static void DeleteEdge(Edge &edge);

        bool operator==(const Edge &other) const;
        bool operator!=(const Edge &other) const;
        bool operator<(const Edge &other) const;
        bool operator>(const Edge &other) const;
        bool operator<=(const Edge &other) const;
    };

//    std::ostream& operator<<(std::ostream& os, const Edge& edge);


    class VertexData {
    private:
        Vertex *vertex;
        std::list<Sequence> hanging{};
    public:
        VertexData(Vertex &vertex) : vertex(&vertex) {}
        const std::list<Sequence> &getHanging() const {return hanging;}
        void addOutgoingSequence(const Sequence &edge_seq);
        void fireAddEdge(Edge &edge);
        void clear() {hanging.clear();}
    };

    class Vertex {
    private:
        friend class SparseDBG;
        friend class Edge;

        mutable std::list<Edge> outgoing_{};
        size_t _outDeg = 0;
        Vertex *rc_;
        hashing::htype hash_;
        omp_lock_t writelock = {};
        bool canonical = false;
        bool mark_ = false;

        explicit Vertex(hashing::htype hash, Vertex *_rc);
        Edge &innerAddEdge(Vertex &end, const Sequence &full_sequence);

        Sequence seq;
        mutable VertexData data;

    public:
        explicit Vertex(hashing::htype hash = 0);
        Vertex(const Vertex &) = delete;

        ~Vertex();

        std::string getInnerId() const;
        std::string oldId() const;
        std::string getShortId() const;

        void lock() { omp_set_lock(&writelock); }
        void unlock() { omp_unset_lock(&writelock); }

        VertexData &getData() const {return data;}
        bool isCanonical() const;

        void mark() { mark_ = true; }
        void unmark() { mark_ = false; }
        bool marked() const { return mark_; }

        hashing::htype hash() const { return hash_; }

        Vertex &rc() { return *rc_; }
        const Vertex &rc() const { return *rc_; }

        void setSeq(Sequence _seq);

        size_t size() const { return seq.size(); }

        const Sequence &getSeq() const { return seq; }

        std::list<Edge>::iterator begin() const { return outgoing_.begin(); }
        std::list<Edge>::iterator end() const { return outgoing_.end(); }
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
        Edge &addEdgeLockFree(Vertex &end, const Sequence &full_sequence);
        Edge &addEdge(Vertex &end, const Sequence &full_sequence);
//        Use this method very careful. It breaks basic graph assumptions but it is often important for parallelization.
//        This method removes one edge in a pair of rc edges. When a thread invokes this method no other thread should
//        try to access rc() method in the rc edge.
        void innerRemoveEdge(Edge &edge);
        void clear();

        bool operator==(const Vertex &other) const;
        bool operator!=(const Vertex &other) const;
        bool operator<(const Vertex &other) const;
        bool operator>(const Vertex &other) const;
    };

    class VertexLocker {
    private:
        std::vector<Vertex *> vertices;
    public:
        explicit VertexLocker(std::vector<Vertex *> _vertices) : vertices(std::move(_vertices)) {
            std::function<bool(Vertex * const &v1, Vertex * const &v2)> f = [](Vertex * const &v1, Vertex * const &v2) {
                return *v1 < *v2;
            };
            std::sort(vertices.begin(), vertices.end(), f);
            vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
            for(Vertex *v: vertices) {
                v->lock();
            }
        }
        VertexLocker(const VertexLocker &other) = delete;
        VertexLocker(VertexLocker &&other) = delete;
        VertexLocker& operator=(const VertexLocker &other) = delete;
        VertexLocker& operator=(VertexLocker &&other) = delete;

        ~VertexLocker() {
            for(Vertex *v: vertices) {
                v->unlock();
            }
        }
    };


    struct EdgePosition {
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

    class SparseDBG {
    public:
        typedef dbg::Vertex Vertex;
        typedef dbg::Edge Edge;
        typedef std::unordered_map<hashing::htype , Vertex, hashing::alt_hasher<hashing::htype>> vertex_map_type;
        typedef std::unordered_map<hashing::htype, Vertex, hashing::alt_hasher<hashing::htype>>::iterator vertex_iterator_type;
        typedef std::unordered_map<hashing::htype, EdgePosition, hashing::alt_hasher<hashing::htype>> anchor_map_type;
    private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
        vertex_map_type v;
        anchor_map_type anchors;
        hashing::RollingHash hasher_;
        bool anchors_filled = false;

//    Be careful since hash does not define vertex. Rc vertices share the same hash
        Vertex &innerAddVertex(hashing::htype h) {
            return v.emplace(std::piecewise_construct, std::forward_as_tuple(h),
                             std::forward_as_tuple(h)).first->second;
        }

    public:

        template<class Iterator>
        SparseDBG(Iterator begin, Iterator end, hashing::RollingHash _hasher) : hasher_(_hasher) {
            while (begin != end) {
                hashing::htype hash = *begin;
                if (v.find(hash) == v.end())
                    addVertex(hash);
                ++begin;
            }
        }
        explicit SparseDBG(hashing::RollingHash _hasher) : hasher_(_hasher) {}
        SparseDBG(SparseDBG &&other) = default;
        SparseDBG &operator=(SparseDBG &&other) = default;
        SparseDBG(const SparseDBG &other) noexcept = delete;

        SparseDBG Subgraph(std::vector<Segment<Edge>> &pieces);
        SparseDBG SplitGraph(const std::vector<EdgePosition> &breaks);
        SparseDBG AddNewSequences(logging::Logger &logger, size_t threads, const std::vector<Sequence> &new_seqs);

        const hashing::RollingHash &hasher() const {return hasher_;}
        bool containsVertex(const hashing::htype &hash) const {return v.find(hash) != v.end();}
        Vertex &getVertex(const hashing::KWH &kwh);
        Vertex &getVertex(const Sequence &seq);
        Vertex &getVertex(hashing::htype hash, bool canonical = true) {return canonical ? v.find(hash)->second : v.find(hash)->second.rc();}
        Vertex &getVertex(const Vertex &other_graph_vertex);
        std::array<Vertex *, 2> getVertices(hashing::htype hash);
//        const Vertex &getVertex(const hashing::KWH &kwh) const;
        bool isAnchor(hashing::htype hash) const {return anchors.find(hash) != anchors.end();}
        EdgePosition getAnchor(const hashing::KWH &kwh);
        size_t size() const {return v.size();}
        bool alignmentReady() const {return anchors_filled;}

        void checkConsistency(size_t threads, logging::Logger &logger);
        void checkDBGConsistency(size_t threads, logging::Logger &logger);
        void checkSeqFilled(size_t threads, logging::Logger &logger);
        void fillAnchors(size_t w, logging::Logger &logger, size_t threads);
        void fillAnchors(size_t w, logging::Logger &logger, size_t threads, const std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> &to_add);
        void processRead(const Sequence &seq);
        void processFullEdgeSequence(const Sequence &old_seq);
        void removeIsolated();
        void removeMarked();

        void resetMarkers();

        void addVertex(hashing::htype h) {innerAddVertex(h);}
        Vertex &addVertex(const hashing::KWH &kwh);
        Vertex &addVertex(const Sequence &seq);
        Vertex &addVertex(const Vertex &other_graph_vertex);
//        Edge &addEdge(Vertex &from, Vertex &to, Sequence seq, int id = 0) {
//            Edge edge(from, to, seq, id);
//        }


        std::vector<hashing::KWH> extractVertexPositions(const Sequence &seq, size_t max = size_t(-1)) const;
        void printFastaOld(const std::experimental::filesystem::path &out);

        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> vertices(bool unique = false) &;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> vertices(bool unique = false) && = delete;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> verticesUnique() &;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> verticesUnique() && = delete;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edges(bool unique = false) &;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edges(bool unique = false) && = delete;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edgesUnique() &;
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edgesUnique() && = delete;
    };

}
