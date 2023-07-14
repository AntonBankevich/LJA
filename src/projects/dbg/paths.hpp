#pragma once

#include "sparse_dbg.hpp"
namespace dbg {
    class GraphPath {
    private:
        Vertex *start_;
        std::vector<Edge *> path;
        size_t cut_left = 0;
        size_t cut_right = 0;
    public:
        typedef typename std::vector<Edge *>::iterator iterator;
        typedef typename std::vector<Edge *>::const_iterator const_iterator;
        typedef TransformingIterator<CountingIterator<size_t>, Vertex> vertex_iterator;
        typedef TransformingIterator<CountingIterator<size_t>, Edge> edge_iterator;
        typedef TransformingGenerator<CountingIterator<size_t>, Segment<Edge>> segment_iterator;

        GraphPath(Vertex &_start, std::vector<Edge *> _path, size_t left_skip, size_t rightSkip) :
                    start_(&_start), path(std::move(_path)), cut_left(left_skip), cut_right(rightSkip) {}
        GraphPath(Vertex &_start) : start_(&_start) {} // NOLINT(google-explicit-constructor)
        GraphPath(Edge &edge) : start_(&edge.getStart()), path({&edge}) {} // NOLINT(google-explicit-constructor)
        GraphPath(const Segment<Edge> &segment) : start_(&segment.contig().getStart()), // NOLINT(google-explicit-constructor)
                          path({&segment.contig()}), cut_left(segment.left), cut_right(segment.contig().size() - segment.right) {}
        GraphPath() : start_(nullptr) {}
        template<class Iterator>
        explicit GraphPath(Iterator begin, Iterator end) : cut_left(0), cut_right(0) {
            while(begin != end) {
                *this += *begin;
                ++begin;
            }
            if(size() > 0) {
                start_ = &frontEdge().getStart();
            }
        }

        static GraphPath WalkForward(Edge &start);

        Vertex &getVertex(size_t i);
        Vertex &getVertex(size_t i) const;
        Vertex &start() const {return *start_;}
        Vertex &finish() const {return path.empty() ? start() : path.back()->getFinish();}
        size_t find(Edge &edge, size_t pos = 0) const;
        size_t find(Vertex &v, size_t pos = 0) const;
        Edge &backEdge() const {return *path.back();}
        Edge &frontEdge() const {return *path.front();}
        size_t size() const {return path.size();}

//        TODO: Find a way to iterate over temporary path objects
        IterableStorage<vertex_iterator> vertices() const &;
        IterableStorage<vertex_iterator> vertices() && = delete;
        IterableStorage<edge_iterator> edges() const &;
        IterableStorage<edge_iterator> edges() && = delete;

        segment_iterator begin() const;
        segment_iterator  end() const;

        GraphPath subPath(size_t from, size_t to);
        GraphPath RC() const;
        double minCoverage() const;
        Sequence Seq() const;
        Sequence truncSeq() const;
        size_t len() const;

        Segment<Edge> back() const;
        Segment<Edge> front() const;
        Segment<Edge> operator[](size_t i) const;
        std::string str(bool show_coverage = false) const;

        bool valid() const;
        void invalidate();
        GraphPath subalignment(size_t from, size_t to) const;
        GraphPath subalignment(size_t from) const {return subalignment(from, size());}
        GraphPath reroute(size_t left, size_t right, const GraphPath &rerouting) const;
        void operator+=(const GraphPath &other);
        void operator+=(const Segment<Edge> &other);
        void operator+=(Edge &other);
        GraphPath operator+(const GraphPath &other) const;
        GraphPath operator+(const Segment<Edge> &other) const;
        GraphPath operator+(Edge &other) const;
        //TODO deprecate
        Sequence map(std::unordered_map<const Edge *, Sequence> &edge_map);

        void pop_back() {
            path.pop_back();
            cut_right = 0;
        }
        void pop_back(size_t len) {
            path.erase(path.end() - len, path.end());
            cut_right = 0;
        }
        void cutBack(size_t l);
        void cutFront(size_t l);
        GraphPath &addStep();
        GraphPath &addStep(Edge &edge);
        std::vector<GraphPath> allSteps();
        std::vector<GraphPath> allExtensions(size_t len);
        GraphPath &extend(const Sequence &seq);

        unsigned char lastNucl() const;
        size_t leftSkip() const;
        size_t rightSkip() const;
        bool endClosed() const;
        bool startClosed() const;

        Sequence truncSubseq(size_t start_position, size_t size = 10000000) const;

        bool operator==(const GraphPath &other) const {
            return start_ == other.start_ && cut_left == other.cut_left && cut_right == other.cut_right && path == other.path;
        }
        bool operator!=(const GraphPath &other) const {return !operator==(other);}
    };


    template<class U, class V>
    class PerfectAlignment {
    public:
        Segment<U> seg_from;
        Segment<V> seg_to;
        PerfectAlignment(const Segment<U> &seg_from_, const Segment<V> &seg_to_) : seg_from(seg_from_), seg_to(seg_to_) {
            VERIFY(seg_from_.size() == seg_to_.size());
        }
        size_t size() {return seg_from.size();}
        PerfectAlignment RC() const {
            return {seg_from.RC(), seg_to.RC()};
        }
        bool operator<(const PerfectAlignment<U, V> &other) const {
            if(seg_to != other.seg_to)
                return seg_to < other.seg_to;
            else
                return seg_from < other.seg_from;
        }
    };

    class GraphAligner {
    private:
        SparseDBG &dbg;
        PerfectAlignment<Contig, dbg::Edge> extendLeft(const hashing::KWH &kwh, Contig &contig) const;
        PerfectAlignment<Contig, dbg::Edge> extendRight(const hashing::KWH &kwh, Contig &contig) const;
    public:
        explicit GraphAligner(SparseDBG &dbg) : dbg(dbg) {
            VERIFY(dbg.alignmentReady());
        }

        GraphPath align(const EdgePosition &pos, const Sequence &seq) const;
        GraphPath align(const Sequence &seq, Edge *edge_to, size_t pos_to);
        GraphPath align(const Sequence &seq, const std::string &name = "") const;
        std::vector<PerfectAlignment<Contig, Edge>> carefulAlign(Contig &contig) const;
        std::vector<PerfectAlignment<Edge, Edge>> oldEdgeAlign(Edge &contig) const;
        std::vector<PerfectAlignment<Contig, Edge>> sparseAlign(Contig &contig) const;
    };
}

namespace std {
    template<class U, class V>
    std::ostream &operator<<(std::ostream &os, const dbg::PerfectAlignment<U, V> &al) {
        return os << al.seg_from << "->" << al.seg_to << std::endl;
    }
}