#pragma once

#include "sparse_dbg.hpp"
namespace dbg {
    class Path {
    private:
        Vertex *start_;
        std::vector<Edge *> path;
    public:
        typedef typename std::vector<Edge *>::iterator iterator;
        typedef typename std::vector<Edge *>::const_iterator const_iterator;

        Path(Vertex &_start, std::vector<Edge *> _path) : start_(&_start), path(std::move(_path)) {}
        explicit Path(Vertex &_start) : start_(&_start) {}
        static Path WalkForward(Edge &start);

        Edge &operator[](size_t i) {return *path[i];}
        Edge &operator[](size_t i) const {return *path[i];}
        Vertex &getVertex(size_t i);
        Vertex &getVertex(size_t i) const;
        Vertex &start() const {return *start_;}
        Vertex &finish() const {return path.empty() ? start() : *path.back()->end();}
        size_t find(Edge &edge, size_t pos = 0) const;
        size_t find(Vertex &v, size_t pos = 0) const;
        Edge &back() {return *path.back();}
        Edge &front() {return *path.front();}
        size_t size() const {return path.size();}
        iterator begin() {return path.begin();}
        iterator end() {return path.end();}
        const_iterator begin() const {return path.begin();}
        const_iterator end() const {return path.end();}

        Path subPath(size_t from, size_t to);
        Path RC();
        double minCoverage() const;
        Sequence Seq() const;
        Sequence truncSeq() const;
        size_t len() const;
        Path operator+(const Path &other) const;
        void operator+=(Edge &edge);
    };


    class GraphAlignment {
    private:
        Vertex *start_;
        std::vector<Segment<Edge>> als;
    public:
//    TODO change interface
        template<class Iterator>
        explicit GraphAlignment(Iterator begin, Iterator end) {
            while(begin != end) {
                als.emplace_back(*begin);
                ++begin;
            }
            if(!als.empty())
                start_ = als.front().contig().start();
        }
        explicit GraphAlignment(std::vector<Segment<Edge>> &&_path) : start_(_path.front().contig().start()), als(std::move(_path)) {}
        GraphAlignment(Vertex *_start, std::vector<Segment<Edge>> &&_path) : start_(_start), als(std::move(_path)) {}
        explicit GraphAlignment(Vertex &_start) : start_(&_start) {}
        GraphAlignment() : start_(nullptr) {}

        Vertex &start() const {return *start_;}
        Vertex &finish() const {return als.empty() ? *start_ : *als.back().contig().end();}
        Segment<Edge> &back() {return als.back();}
        Segment<Edge> &front() {return als.front();}
        const Segment<Edge> &back() const {return als.back();}
        const Segment<Edge> &front() const {return als.front();}
        const Segment<Edge> &operator[](size_t i) const {return als[i];}
        Segment<Edge> &operator[](size_t i) {return als[i];}
        Vertex &getVertex(size_t i) const {return i == 0 ? *start_ : *als[i - 1].contig().end();}
        Vertex &getVertex(size_t i)  {return i == 0 ? *start_ : *als[i - 1].contig().end();}
        size_t find(Edge &edge, size_t pos = 0) const;
        size_t find(Vertex &v, size_t pos = 0) const;
        typename std::vector<Segment<Edge>>::iterator begin() {return als.begin();}
        typename std::vector<Segment<Edge>>::iterator end() {return als.end();}
        typename std::vector<Segment<Edge>>::const_iterator begin() const {return als.begin();}
        typename std::vector<Segment<Edge>>::const_iterator end() const {return als.end();}
        Path path();
        size_t size() const {return als.size();}
        std::string str(bool show_coverage = false) const;
        size_t len() const;

        GraphAlignment RC() const;
        bool valid() const;
        void invalidate();
        GraphAlignment subalignment(size_t from, size_t to) const;
        GraphAlignment subalignment(size_t from) const {return subalignment(from, size());}
        GraphAlignment reroute(size_t left, size_t right, const GraphAlignment &rerouting) const;
        GraphAlignment reroute(size_t left, size_t right, const Path &rerouting) const;
        void operator+=(const Path &other);
        void operator+=(const GraphAlignment &other);
        void operator+=(const Segment<Edge> &other);
        void operator+=(Edge &other);
        GraphAlignment operator+(const GraphAlignment &other) const;
        GraphAlignment operator+(const Path &other) const;
        GraphAlignment operator+(const Segment<Edge> &other) const;
        GraphAlignment operator+(Edge &other) const;
        //TODO deprecate
        Sequence map(std::unordered_map<const Edge *, Sequence> &edge_map);

        void push_back(const Segment<Edge> &seg) {als.push_back(seg);}
        void pop_back() {als.pop_back();}
        void pop_back(size_t len) {als.erase(als.end() - len, als.end());}
        void cutBack(size_t l);
        GraphAlignment &addStep();
        GraphAlignment &addStep(Edge &edge);
        std::vector<GraphAlignment> allSteps();
        std::vector<GraphAlignment> allExtensions(size_t len);
        GraphAlignment &extend(const Sequence &seq);

        unsigned char lastNucl() const;
        size_t leftSkip() const;
        size_t rightSkip() const;
        bool endClosed() const;
        bool startClosed() const;
        double minCoverage() const;

        Sequence Seq() const;
        Sequence truncSeq() const;
        Sequence truncSeq(size_t start_position, size_t size = 10000000) const;

        bool operator==(const GraphAlignment &other) const {return start_ == other.start_ && als == other.als;}
        bool operator!=(const GraphAlignment &other) const {return !operator==(other);}
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
    };

    class GraphAligner {
    private:
        SparseDBG &dbg;

    public:
        explicit GraphAligner(SparseDBG &dbg) : dbg(dbg) {
        }

        GraphAlignment align(const EdgePosition &pos, const Sequence &seq) const;
        GraphAlignment align(const Sequence &seq, Edge *edge_to, size_t pos_to);
        GraphAlignment align(const Sequence &seq) const;
        std::vector<PerfectAlignment<Contig, Edge>> carefulAlign(Contig &contig) const;
        std::vector<PerfectAlignment<Edge, Edge>> oldEdgeAlign(Edge &contig) const;
    };
}

namespace std {
    template<class U, class V>
    std::ostream &operator<<(std::ostream &os, const dbg::PerfectAlignment<U, V> &al) {
        return os << al.seg_from << "->" << al.seg_to << std::endl;
    }
}