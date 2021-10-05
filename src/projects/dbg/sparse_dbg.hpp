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
#include <vector>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace dbg {
    class Vertex;

    class SparseDBG;

    class Edge {
    private:
        Vertex *start_;
        Vertex *end_;
        mutable size_t cov;
    public:
        mutable size_t extraInfo;
        Sequence seq;
        std::string id = "";
        friend class Vertex;
        bool is_reliable = false;
        Edge(Vertex *_start, Vertex *_end, const Sequence &_seq) :
                start_(_start), end_(_end), cov(0), extraInfo(-1), seq(_seq) {
        }
        std::string getId() const;
        std::string oldId() const;
        std::string getShortId() const;
        Vertex *end() const;
        Vertex *start() const;
        size_t getTipSize() const;
        size_t updateTipSize() const;
        void bindTip(Vertex &start, Vertex &end);
        size_t common(const Sequence &other) const;
        size_t size() const;
        double getCoverage() const;
        size_t intCov() const;
        Edge &rc() const;
        Edge &sparseRcEdge() const;
        void incCov(size_t val) const;
        Sequence firstNucl() const;
        Sequence kmerSeq(size_t pos) const;
        Sequence suffix(size_t pos) const;
        std::string str() const;
        bool operator==(const Edge &other) const;
        bool operator!=(const Edge &other) const;
        bool operator<(const Edge &other) const;
        bool operator>(const Edge &other) const;
        bool operator<=(const Edge &other) const;
    };

//    std::ostream& operator<<(std::ostream& os, const Edge& edge);



    class Vertex {
    private:
        friend class SparseDBG;
        mutable std::vector<Edge> outgoing_{};
        Vertex *rc_;
        hashing::htype hash_;
        omp_lock_t writelock = {};
        size_t coverage_ = 0;
        bool canonical = false;
        bool mark_ = false;
        explicit Vertex(hashing::htype hash, Vertex *_rc);
    public:
        Sequence seq;

        explicit Vertex(hashing::htype hash = 0);
        Vertex(const Vertex &) = delete;
        ~Vertex();

        void mark() {mark_ = true;}
        void unmark() {mark_ = false;}
        bool marked() const {return mark_;}
        hashing::htype hash() const {return hash_;}
        Vertex &rc() {return *rc_;}
        const Vertex &rc() const {return *rc_;}
        void setSequence(const Sequence &_seq);
        void lock() {omp_set_lock(&writelock);}
        void unlock() {omp_unset_lock(&writelock);}
        std::vector<Edge>::iterator begin() const {return outgoing_.begin();}
        std::vector<Edge>::iterator end() const {return outgoing_.end();}
        size_t outDeg() const {return outgoing_.size();}
        size_t inDeg() const {return rc_->outgoing_.size();}
        Edge &operator[](size_t ind) const {return outgoing_[ind];}


        size_t coverage() const;
        bool isCanonical() const;
        bool isCanonical(const Edge &edge) const;
        void clear();
        void clearOutgoing();
        void sortOutgoing();
        void checkConsistency() const;
        std::string getId() const;
        std::string oldId() const;
        std::string getShortId() const;
        void incCoverage();
        void clearSequence();
        Edge &addEdgeLockFree(const Edge &edge);
        void addEdge(const Edge &e);
        Edge &getOutgoing(unsigned char c) const;
        bool hasOutgoing(unsigned char c) const;
        bool isJunction() const;

        bool operator==(const Vertex &other) const;
        bool operator!=(const Vertex &other) const;
        bool operator<(const Vertex &other) const;
        bool operator>(const Vertex &other) const;
    };

    struct EdgePosition {
        Edge *edge;
        size_t pos;

        EdgePosition(Edge &_edge, size_t _pos) : edge(&_edge), pos(_pos) {VERIFY(pos >= 0 && pos <= edge->size());}
        EdgePosition() : edge(nullptr), pos(0) {}

        Sequence kmerSeq() const {return edge->kmerSeq(pos);}
        unsigned char lastNucl() const {return edge->seq[pos - 1];}
        bool isBorder() const {return pos == 0 || pos == edge->size();}

        std::vector<EdgePosition> step() const;
        EdgePosition RC() const {return {edge->rc(), edge->size() - pos};}
    };

    class SparseDBG {
    public:
        typedef std::unordered_map<hashing::htype , Vertex, hashing::alt_hasher<hashing::htype>> vertex_map_type;
        typedef std::unordered_map<hashing::htype, Vertex, hashing::alt_hasher<hashing::htype>>::iterator vertex_iterator_type;
        typedef std::unordered_map<hashing::htype, EdgePosition, hashing::alt_hasher<hashing::htype>> anchor_map_type;
    private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
        vertex_map_type v;
        anchor_map_type anchors;
        hashing::RollingHash hasher_;

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

        const hashing::RollingHash &hasher() const {return hasher_;}
        bool containsVertex(const hashing::htype &hash) const {return v.find(hash) != v.end();}
        Vertex &getVertex(const hashing::KWH &kwh);
        Vertex &getVertex(const Sequence &seq);
        Vertex &getVertex(hashing::htype hash) {return v.find(hash)->second;}
        Vertex &getVertex(const Vertex &other_graph_vertex);
        std::array<Vertex *, 2> getVertices(hashing::htype hash);
//        const Vertex &getVertex(const hashing::KWH &kwh) const;
        bool isAnchor(hashing::htype hash) const {return anchors.find(hash) != anchors.end();}
        EdgePosition getAnchor(const hashing::KWH &kwh);
        size_t size() const {return v.size();}

        void checkConsistency(size_t threads, logging::Logger &logger);
        void checkSeqFilled(size_t threads, logging::Logger &logger);
        void fillAnchors(size_t w, logging::Logger &logger, size_t threads);
        void fillAnchors(size_t w, logging::Logger &logger, size_t threads, const std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> &to_add);
        void processRead(const Sequence &seq);
        void processEdge(Vertex &vertex, Sequence old_seq);
        Vertex &bindTip(Vertex &start, Edge &tip);
        void removeIsolated();
        void removeMarked();

        void addVertex(hashing::htype h) {innerAddVertex(h);}
        Vertex &addVertex(const hashing::KWH &kwh);
        Vertex &addVertex(const Sequence &seq);
        Vertex &addVertex(const Vertex &other_graph_vertex);


        std::vector<hashing::KWH> extractVertexPositions(const Sequence &seq, size_t max = size_t(-1)) const;
        void printFastaOld(const std::experimental::filesystem::path &out);

        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> vertices(bool unique = false);
        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> verticesUnique();
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edges(bool unique = false);
        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edgesUnique();
        typename vertex_map_type::iterator begin() {return v.begin();}
        typename vertex_map_type::iterator end() {return v.end();}
        typename vertex_map_type::const_iterator begin() const {return v.begin();}
        typename vertex_map_type::const_iterator end() const {return v.end();}
    };
}
