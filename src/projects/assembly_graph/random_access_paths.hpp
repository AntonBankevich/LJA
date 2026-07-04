#pragma once

#include "assembly_graph/assembly_graph_base.hpp"
#include "sequences/contigs.hpp"
namespace ag {
class RAPathDirection;

class RAGraphPath {
public:
    friend class RAPathIterator;

    friend class RAPathDirection;

private:
    VertexId start_;
    std::deque<EdgeId> path;
    size_t cut_left;
    size_t cut_right;

    void set(Edge &edge, size_t pos) { path[pos] = edge.getId(); }

public:
    typedef typename std::vector<EdgeId>::iterator iterator;
    typedef typename std::vector<EdgeId>::const_iterator const_iterator;
    typedef TransformingIterator<CountingIterator<size_t>, Vertex> vertex_iterator;
    typedef TransformingIterator<CountingIterator<size_t>, Edge> edge_iterator;
    typedef Generator<CountingIterator<size_t>, Segment<Edge>> segment_iterator;

    RAGraphPath(Vertex &_start, std::vector<EdgeId> _path, size_t cut_left, size_t cut_right) :
            start_(_start.getId()), path(_path.begin(), _path.end()), cut_left(cut_left),
            cut_right(cut_right) {}

    RAGraphPath(Vertex &_start, size_t cut_left = 0, size_t cut_right = 0) : start_(
            _start.getId()), // NOLINT(google-explicit-constructor)
                                                                             cut_left(cut_left), cut_right(
                    cut_right) {} // NOLINT(google-explicit-constructor)
    RAGraphPath(Edge &edge, size_t cut_left = 0, size_t cut_right = 0) : start_(
            edge.getStart().getId()), // NOLINT(google-explicit-constructor)
                                                                         path({edge.getId()}),
                                                                         cut_left(cut_left),
                                                                         cut_right(
                                                                                 cut_right) {} // NOLINT(google-explicit-constructor)
    RAGraphPath(const Segment<Edge> &segment) : start_( // NOLINT(google-explicit-constructor)
            segment.contig().getStart().getId()), // NOLINT(google-explicit-constructor)
                                                path({segment.contig().getId()}), cut_left(segment.left),
                                                cut_right(segment.contig().truncSize() - segment.right) {}

    RAGraphPath() : start_({}), cut_left(0), cut_right(0) {}

    std::vector<EdgeId> asEdgeIds() const;

    template<class Iterator>
    explicit RAGraphPath(Iterator begin, Iterator end);

    static RAGraphPath WalkForward(Edge &start);

    void normalize();

    size_t find(Edge &edge, size_t pos = 0) const;

    size_t find(Vertex &v, size_t pos = 0) const;

    Vertex &getVertex(size_t i) const;

    Edge &getEdge(size_t i) const;

    Segment<Edge> operator[](size_t i) const;

    Vertex &getStart() const { return *start_; }

    Vertex &getFinish() const { return empty() ? getStart() : backEdge().getFinish(); }

    Edge &backEdge() const { return *path.back(); }

    Edge &frontEdge() const { return *path.front(); }

    Segment<Edge> back() const;

    Segment<Edge> front() const;

    size_t size() const { return path.size(); }

    bool empty() const { return size() == 0; }

    bool isSingleton() const { return size() == 1; }

    bool valid() const;

    size_t leftCut() const;

    size_t rightCut() const;

    bool endClosed() const;

    bool startClosed() const;

    //        TODO: Find a way to iterate over temporary path objects
    IterableStorage<vertex_iterator> vertices(size_t from = 0, size_t to = size_t(-1)) const &;

    IterableStorage<vertex_iterator> innerVertices() const &;

    IterableStorage<vertex_iterator> vertices() && = delete;

    IterableStorage<edge_iterator> edges() const &;

    IterableStorage<edge_iterator> edges() && = delete;

    const std::deque<EdgeId> edgeIds() const &;

    segment_iterator begin() const;

    segment_iterator end() const;

    RAPathDirection forward();

    RAPathDirection backward();

//        Reverses piece order and RCs each piece: for path [p0, p1, ..., pk], RC() == [pk.rc(), ...,
//        p0.rc()] — so RC().frontEdge() == backEdge().rc() and RC().backEdge() == frontEdge().rc().
    RAGraphPath RC() const;

    Sequence Seq() const;

    Sequence truncSeq() const;

    RAGraphPath subPath(size_t from, size_t to) const;

    RAGraphPath subPath(size_t from) const { return subPath(from, size()); }

    size_t truncLen() const;

    size_t len() const;

    std::string str() const;


    void invalidate();

    RAGraphPath reroute(size_t left, size_t right, const RAGraphPath &rerouting) const;

    void simpleReroute(size_t left, size_t right, const std::vector<EdgeId> &edges);

    void setCutLeft(size_t value) { cut_left = value; }

    void setCutRight(size_t value) { cut_right = value; }

    void operator+=(const RAGraphPath &other);

    void operator+=(const Segment<Edge> &other);

    void operator+=(Edge &other);

    void pop_back();

    void pop_back(size_t len);

    void pop_front();

    void pop_front(size_t len);

    RAGraphPath &cutBack(size_t l);

    RAGraphPath &cutFront(size_t l);

    RAGraphPath &addStep();

    RAGraphPath &addStep(Edge &edge);

    std::vector<RAGraphPath> allSteps();

    std::vector<RAGraphPath> allExtensions(size_t len);

    RAGraphPath &extend(const Sequence &seq);

    RAGraphPath &fastExtend(const Sequence &seq);

    RAGraphPath operator+(const RAGraphPath &other) const;

    RAGraphPath operator+(const Segment<Edge> &other) const;

    RAGraphPath operator+(Edge &other) const;

    //TODO deprecate
    Sequence map(std::unordered_map<const Edge *, Sequence> &edge_map);

    bool operator==(const RAGraphPath &other) const;

    bool operator!=(const RAGraphPath &other) const;
};

class RAPathIterator {
private:
    RAGraphPath *path;
    bool rc;
    int pos;
public:
    RAPathIterator(RAGraphPath &path, bool rc, int pos) : path(&path), rc(rc), pos(pos) {
    }

public:
    typedef Edge &reference;
    typedef Edge *pointer;

    reference operator*() const { return rc ? path->path[pos]->rc() : *path->path[pos]; }

    pointer operator->() const { return rc ? &path->path[pos]->rc() : &(*path->path[pos]); }

    void set(Edge &edge) const { path->set(rc ? edge.rc() : edge, pos); }

    RAPathIterator &operator++();

    RAPathIterator &operator--();

    RAPathIterator operator+(int d) const { return {*path, rc, rc ? pos - d : pos + d}; }

    RAPathIterator operator-(int d) const { return operator+(-d); }

    RAPathIterator operator++(int) const { return *this + 1; }

    RAPathIterator operator--(int) const { return *this - 1; }

    bool operator==(const RAPathIterator &other) const;

    bool operator!=(const RAPathIterator &other) const { return !(*this == other); }

    int getPos() const { return pos; }

    std::string str() const;

//        Careful! RC does not point to the same element. Instead it makes sure that begin->end and end->begin
    RAPathIterator RC() const { return {*path, !rc, pos + (rc ? 1 : -1)}; }

    RAPathIterator SameElementRC() const { return {*path, !rc, pos}; }
};

class RAPathDirection {
private:
    RAGraphPath *path;
    bool rc;

    void setCutLeft(size_t val) const { rc ? path->cut_right = val : path->cut_left = val; }

    void setCutRight(size_t val) const { rc ? path->cut_left = val : path->cut_right = val; }

public:
    RAPathDirection(RAGraphPath &path, bool rc) : path(&path), rc(rc) {}

    Vertex &getStart() const { return rc ? path->getFinish().rc() : path->getStart(); }

    Vertex &getFinish() const { return rc ? path->getStart().rc() : path->getFinish(); }

    size_t cutRight() const { return rc ? path->leftCut() : path->rightCut(); }

    size_t cutLeft() const { return rc ? path->rightCut() : path->leftCut(); }

    void pop_front() const { return rc ? path->pop_back() : path->pop_front(); }

    void pop_back() const { return rc ? path->pop_front() : path->pop_back(); }

    bool isForward() const { return !rc; }

    size_t size() const { return path->size(); }

    RAGraphPath &getPath() const { return *path; }

    RAPathIterator end() const;

    RAPathIterator begin() const;

    RAPathDirection RC() const { return {*path, !rc}; }

    bool valid() const { return path->valid(); }

    bool empty() const { return path->empty(); }

    bool isRC() const { return rc; }

    void invalidate() { path->invalidate(); }

    void rerouteSameSize(RAPathIterator from, RAPathIterator to, const ag::RAGraphPath &alt) const;

    void cutFrontPrefix() const;

    void cutBackSuffix() const;

    void forceReroute(RAPathIterator from, RAPathIterator to, const RAGraphPath &alt) const;
};

    inline std::ostream &operator<<(std::ostream &os, const RAGraphPath &path) {return os << path.str();}

}
