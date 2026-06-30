#pragma once
#include <utility>

#include "assembly_graph_base.hpp"
#include "random_access_paths.hpp"

namespace ag {

    class PathDirection;
    class ConstPathDirection;
    class PathIterator;
    class PathVertexIterator;
    class SegmentIterator;
    class PathPosition;

//    Invariants for normal paths. They do not have to be always fulfilled but need to keep track of them.
//    1. If edge in the path is neither first, nor last, its full sequence is a substring of the path sequence.
//    2. No prefix edges on start.
//    3. No suffix edges on end.
//    4. In supregraph any path either contains core vertex or consists of a single vertex,
//    which is one of the minimal (by inclusion) vertices, containing the sequence.
//    5. In DBG any path is the shortest path containing the sequence
//    6. Legacy path may consist of a single vertex that is no longer present in the graph
    class GraphPath {
    public:
//        friend class PathIterator;
        friend class PathDirection;
        friend class PathVertexIterator;
        friend class PathPosition;

    private:
        VertexId start = {}; //The first vertex in the path
        VertexId rc_start = {}; //The rc to the last vertex in the path
        NuclDeck fsplits; //Stores concatenated codes for all non-suffix edges
        NuclDeck rsplits; //Stores concatenated codes for all non-suffix edges in rc path
//        In DBG every edge of the path is recorded in both fsplits and rsplits.
//        In SPG every edge is in exactly one of them
        size_t cut_left = 0;
        size_t cut_right = 0;

        GraphPath(VertexId start, VertexId rc_start, NuclDeck fsplits, NuclDeck rsplits, size_t cut_left = 0, size_t cut_right = 0) :
                start(start), rc_start(rc_start), fsplits(std::move(fsplits)), rsplits(std::move(rsplits)),
                cut_left(cut_left), cut_right(cut_right) {}
        GraphPath(Vertex &start, Vertex &rc_start, NuclDeck fsplits, NuclDeck rsplits, size_t cut_left = 0, size_t cut_right = 0) :
                start(start.getId()), rc_start(rc_start.getId()), fsplits(std::move(fsplits)), rsplits(std::move(rsplits)),
                cut_left(cut_left), cut_right(cut_right) {}

    public:
        GraphPath(const GraphPath &) = default;
        GraphPath(GraphPath &&other)  noexcept {
            *this = std::move(other);
        }
        GraphPath &operator=(const GraphPath &) = default;
        GraphPath &operator=(GraphPath &&other) noexcept  = default;
        explicit GraphPath(const RAGraphPath &path, size_t cut_left = 0, size_t cut_right = 0);

//        typedef typename std::vector<EdgeId>::iterator iterator;
//        typedef typename std::vector<EdgeId>::const_iterator const_iterator;
//        typedef TransformingIterator<CountingIterator<size_t>, Vertex> vertex_iterator;
//        typedef TransformingIterator<CountingIterator<size_t>, Edge> edge_iterator;
//        typedef Generator<CountingIterator<size_t>, Segment<Edge>> segment_iterator;


        static GraphPath LegacyPath(VertexId start, VertexId rc_start, size_t cut_left = 0, size_t cut_right = 0) {
            return {start.legacyId(), rc_start.legacyId(), NuclDeck(), NuclDeck(), cut_left, cut_right};
        }

        GraphPath(Vertex &start, size_t cut_left = 0, size_t cut_right = 0) : start(start.getId()), rc_start(start.rc().getId()), // NOLINT(google-explicit-constructor)
                                            cut_left(cut_left), cut_right(cut_right) {
            VERIFY(this->start == start.getId());
        }
        GraphPath(Vertex &start, const Sequence &extension);
        GraphPath(Edge &edge, size_t cut_left = 0, size_t cut_right = 0); // NOLINT(google-explicit-constructor)
        GraphPath(const Segment<Edge> &segment) : GraphPath(segment.contig(), segment.cutLeft(), segment.cutRight()) {} // NOLINT(google-explicit-constructor)

        GraphPath() : start({}), rc_start({}), fsplits(), rsplits(), cut_left(0), cut_right(0) {}

        template<class Iterator>
        explicit GraphPath(Iterator begin, Iterator end);

        std::string str() const;
        RAGraphPath asRAPath() const;
        size_t calculateSize() const;

        bool isLegacy() const {return valid() && start.isLegacy();}
        bool valid() const {return start != VertexId();}
        size_t leftCut() const {return cut_left;}
        size_t rightCut() const {return cut_right;}
        bool endClosed() const {return cut_right == 0;}
        bool startClosed() const {return cut_left == 0;}
        Edge &backEdge() const { return rc_start->outDeg() == 1 ? rc_start->front().rc() : rc_start->getOutgoingByIterator(rsplits.begin()).rc(); }
        Edge &frontEdge() const { return start->outDeg() == 1 ? start->front() : start->getOutgoingByIterator(fsplits.begin()); }
        Vertex &getStart() const { return *start; }
        Vertex &getFinish() const { return rc_start->rc(); }
        Vertex &getRCStart() const {return *rc_start;}
        Segment<Edge> getSegment(PathPosition position) const;
        bool empty() const { return fsplits.empty() && rsplits.empty(); }
//        It is important to keep track of how this method behaves when graph is modified in parallel.
        bool isSingleton() const;
        Segment<Edge> back() const;
        Segment<Edge> front() const;
        PathPosition firstPosition() const;
        PathPosition lastPosition() const;
        PathPosition endPosition() const;
        PathPosition rcEndPosition() const;
        const NuclDeck &getFSplits() const {return fsplits;}
        const NuclDeck &getRSplits() const {return rsplits;}
        bool operator==(const GraphPath &other) const;
        bool operator!=(const GraphPath &other) const { return !operator==(other); }

        PathDirection forward();
        PathDirection backward();
        ConstPathDirection forward() const;
        ConstPathDirection backward() const;

        void pop_back();
        void pop_front();
        void pop_front(Edge &edge);
        void pop_back(Edge &edge);
        void replace_back(Edge &edge);
        void replace_front(Edge &edge);
        void normalize();
        void invalidate() {*this = {};}
        void setCutLeft(size_t value);
        void setCutRight(size_t value);
        void operator+=(const GraphPath &other);
        void operator+=(const Segment<Edge> &other);
        void operator+=(Edge &other);
        void push_front(Edge &other);
        void push_front(const Segment<Edge> &other);
        void push_front(const GraphPath &other);
        void pop_back(size_t len);
        void pop_front(size_t len);
        GraphPath &shorten(PathPosition from, PathPosition to);
        GraphPath &cutBack(size_t l);
        GraphPath &cutFront(size_t l);
        ag::GraphPath &addStep() {cut_right -= 1; return *this;}
        ag::GraphPath &addStep(Edge &edge) {*this += Segment<Edge>(edge, 0, 1);return *this;}
        ag::GraphPath & extend(const Sequence &seq);
        ag::GraphPath &fastExtend(const Sequence &seq);
        void resetEdgeCodes();

        void forcePushBack(Edge &edge);
        void forcePushFront(Edge &edge);

        //        TODO: Find a way to iterate over temporary path objects
        IterableStorage<PathVertexIterator> vertices() const &;
        IterableStorage<PathVertexIterator> innerVertices() const &;
        IterableStorage<PathVertexIterator> vertices() && = delete;
        IterableStorage<PathVertexIterator> innerVertices() && = delete;
        IterableStorage<PathIterator> edges() const &;
        IterableStorage<PathIterator> edges() && = delete;
        SegmentIterator begin() const;
        SegmentIterator end() const;

//        TODO: minimaze usage of this function. Use directions instead.
        GraphPath RC() const {return {rc_start, start, rsplits, fsplits, cut_right, cut_left};}
        GraphPath subPath(PathPosition from, PathPosition to) const;
        GraphPath subPath(PathPosition from) const;
        GraphPath operator+(const GraphPath &other) const;
        GraphPath operator+(const Segment<Edge> &other) const;
        GraphPath operator+(Edge &other) const;
        GraphPath operator*(size_t mult) const;
        Sequence Seq() const;
        Sequence truncSeq() const;
        size_t truncLen() const;
        size_t len() const;
        static GraphPath Load(std::istream &os, const IdIndex<Vertex> &index);

//        These functions should return the same as if sequences of the paths were compared
        bool startsWith(const GraphPath &other) const;
        bool endsWith(const GraphPath &other) const;
        bool nonContradicts(const GraphPath &other) const;
    };


    inline std::ostream &operator<<(std::ostream &os, const GraphPath &path) {
        // Change this when vertex hierarchy is saved
        if (path.valid() && !path.isLegacy())
            return os << path.getStart().getInnerId() << " " << path.getFinish().rc().getInnerId() << " F:" <<
                      path.getFSplits() << " R:" << path.getRSplits() << " " << path.leftCut() << " " << path.rightCut();
        else
            return os << "0 0 F: R: 0 0";
    }

    class PathPosition {
        friend class GraphPath;
    protected:
        VertexId vid;
//        fpos and rpos correspond to the iteration direction.
        NuclDeck::Iterator rpos;
        NuclDeck::Iterator fpos;

        void move(int d) {
            if(d > 0)
                for(size_t i = 0; i < d; i++)
                    operator++();
            else
                for(size_t i = 0; i < -d; i++)
                    operator--();
        }
    public:
        PathPosition(VertexId cur, NuclDeck::Iterator fpos, NuclDeck::Iterator rpos) : vid(cur), fpos(fpos), rpos(rpos) {
        }
        PathPosition(Vertex &cur, NuclDeck::Iterator fpos, NuclDeck::Iterator rpos) : vid(cur.getId()), fpos(fpos), rpos(rpos) {
        }
        PathPosition(NuclDeck::Iterator fpos, NuclDeck::Iterator rpos) : vid({}), fpos(fpos), rpos(rpos) {
        }
        PathPosition(const PathPosition &) = default;
        PathPosition(PathPosition &&)  noexcept = default;
        PathPosition &operator=(const PathPosition &) = default;
        PathPosition &operator=(PathPosition &&) = default;

        Edge &nextEdge() const {return vid->getOutgoingByIterator(fpos);}
        Edge &prevEdge() const {return vid->rc().getOutgoingByIterator(rpos).rc();}
        Vertex &getVertex() const {return *vid;}
        NuclDeck::Iterator getFPos() const {return fpos;}
        NuclDeck::Iterator getRPos() const {return rpos;}


//        This method is for minor optimization. we avoid calling nextEdge() if we already know the next edge.

        PathPosition &operator+=(Edge &edge) {
//            for(size_t i = 0; i < edge.getCode().size(); i++) {
//                VERIFY(*fpos == edge.getCode()[i]);
//                ++fpos;
//            }
//            for(size_t i = 0; i < edge.rc().getCode().size(); i++) {
//                --rpos;
//                VERIFY(*rpos == edge.rc().getCode()[edge.rc().getCode().size() - 1 - i]);
//            }
//            vid = edge.getFinish().getId();
            fpos += edge.getCode().size();
            rpos -= edge.rc().getCode().size();
            vid = edge.getFinish().getId();
            return *this;
        }

        PathPosition &operator++() {
            return operator+=(nextEdge());
        }

        PathPosition &operator--() {
            Edge &edge = prevEdge();
            fpos -= edge.getCode().size();
            rpos += edge.rc().getCode().size();
            vid = edge.getStart().getId();
            return *this;
        }

        PathPosition operator+(int d) const {
            PathPosition res = *this;
            res.move(d);
            return res;
        }

        PathPosition operator-(int d) const {return *this + (-d);}
        PathPosition RC() const {return {vid->rc(), rpos, fpos};}

        bool operator==(const PathPosition &other) const {
            VERIFY(&fpos.getDeck()==&other.fpos.getDeck());
            VERIFY(&rpos.getDeck()==&other.rpos.getDeck());
            return fpos == other.fpos && rpos == other.rpos;
        }
        bool operator<=(const PathPosition &other) const {return fpos <= other.fpos && rpos >= other.rpos;}
        bool operator!=(const PathPosition &other) const { return !(*this == other); }
        bool operator<(const PathPosition &other) const {return *this <= other && *this != other;}
    };

    class PathDirection {
    private:
        GraphPath *path;
        bool rc;
    public:
        PathDirection() : path(), rc(false) {}

        PathDirection(GraphPath &path, bool rc) : path(&path), rc(rc) {}
        PathDirection &operator=(const GraphPath &other);
        PathDirection &operator=(GraphPath &&other);

        Vertex &getStart() const { return rc ? path->getFinish().rc() : path->getStart(); }
        Vertex &getFinish() const { return rc ? path->getStart().rc() : path->getFinish(); }
        Segment<Edge> front() const { return rc ? path->back().RC() : path->front(); }
        Segment<Edge> back() const { return rc ? path->front().RC() : path->back(); }
        Edge &frontEdge() const { return rc ? path->backEdge().rc() : path->frontEdge(); }
        Edge &backEdge() const { return rc ? path->frontEdge().rc() : path->backEdge(); }
        size_t rightCut() const { return rc ? path->leftCut() : path->rightCut(); }
        size_t leftCut() const { return rc ? path->rightCut() : path->leftCut(); }
        const NuclDeck &getFSplits() const {return rc ? path->getRSplits() : path->getFSplits();}
        const NuclDeck &getRSplits() const {return rc ? path->getFSplits() : path->getRSplits();}


//        const NuclDeck &getFsplits() const {return rc ? path->rsplits : path->fsplits;}
//        const NuclDeck &getRsplits() const {return rc ? path->fsplits : path->rsplits;}
        Segment<Edge> getSegment(PathPosition position) const;

        IterableStorage<PathVertexIterator> vertices();

        void pop_front() const { return rc ? path->pop_back() : path->pop_front(); }
        void pop_back() const { return rc ? path->pop_front() : path->pop_back(); }
        void pop_front(Edge &edge) const { return rc ? path->pop_back(edge.rc()) : path->pop_front(edge); }
        void pop_back(Edge &edge) const { return rc ? path->pop_front(edge.rc()) : path->pop_back(edge); }
        void push_front(Edge &edge) const {if(rc) *path += edge.rc(); else path->push_front(edge);}
        void push_back(Edge &edge) const {if(rc) path->push_front(edge.rc()); else *path += edge;}
        void replace_front(Edge &edge) const {if(rc) path->replace_back(edge.rc()); else path->replace_front(edge);}
        void replace_back(Edge &edge) const {if(rc) path->replace_front(edge.rc()); else path->replace_back(edge);}
        void operator+=(Edge &edge) const { push_back(edge);}
        void operator+=(Segment<Edge> seg) const { if(rc) path->push_front(seg.RC()); else *path += seg;}
        void operator+=(const GraphPath &other) const { if(rc) path->push_front(other.RC()); else *path += other;}
        void push_front(Segment<Edge> seg) const {if(rc) *path += seg.RC(); else path->push_front(seg);}
        void push_front(const GraphPath &other) const { if(rc) *path += other.RC(); else path->push_front(other);}
        void invalidate() const { path->invalidate(); }
        void setCutLeft(size_t val) const;
        void setCutRight(size_t val) const;
        void forcePushBack(Edge &edge) const;
        void forcePushFront(Edge &edge) const;

        PathPosition firstPosition() const {return rc ? path->lastPosition().RC() : path->firstPosition();}
        PathPosition lastPosition() const {return rc ? path->firstPosition().RC() : path->lastPosition();}
        PathPosition endPosition() const {return rc ? path->rcEndPosition() : path->endPosition();}
        PathPosition rcEndPosition() const {return rc ? path->endPosition() : path->rcEndPosition();}

        bool isForward() const { return !rc; }
        size_t calculateSize() const { return path->calculateSize(); }
        GraphPath &getPath() const { return *path; }
        bool valid() const { return path->valid(); }
        bool empty() const { return path->empty(); }
        bool isSingleton() const;
        bool isRC() const { return rc; }

        PathIterator begin() const;
        PathIterator end() const;

        PathDirection RC() const { return {*path, !rc}; }

//        void forceReroute(PathIterator from, PathIterator to, const GraphPath &alt) const;
        bool operator==(const PathDirection &other) const { return path == other.path && rc == other.rc; }
        bool operator!=(const PathDirection &other) const { return !(*this == other); }
    };

    class ConstPathDirection {
    private:
        const GraphPath *path;
        bool rc;

    public:
        ConstPathDirection() : path(), rc(false) {}
        ConstPathDirection(const PathDirection &other) : path(&other.getPath()), rc(other.isRC()) {}

        ConstPathDirection(const GraphPath &path, bool rc) : path(&path), rc(rc) {}

        Vertex &getStart() const { return rc ? path->getFinish().rc() : path->getStart(); }
        Vertex &getFinish() const { return rc ? path->getStart().rc() : path->getFinish(); }
        Segment<Edge> front() const { return rc ? path->back().RC() : path->front(); }
        Segment<Edge> back() const { return rc ? path->front().RC() : path->back(); }
        Edge &frontEdge() const { return rc ? path->backEdge().rc() : path->frontEdge(); }
        Edge &backEdge() const { return rc ? path->frontEdge().rc() : path->backEdge(); }
        size_t cutRight() const { return rc ? path->leftCut() : path->rightCut(); }
        size_t cutLeft() const { return rc ? path->rightCut() : path->leftCut(); }
        PathPosition firstPosition() const {return rc ? path->lastPosition().RC() : path->firstPosition();}
        PathPosition lastPosition() const {return rc ? path->firstPosition().RC() : path->lastPosition();}
        PathPosition endPosition() const {return rc ? path->rcEndPosition() : path->endPosition();}
        PathPosition rcEndPosition() const {return rc ? path->endPosition() : path->rcEndPosition();}
        Segment<Edge> getSegment(PathPosition position) const;
        bool isForward() const { return !rc; }
        size_t calculateSize() const { return path->calculateSize(); }
        const GraphPath &getPath() const { return *path; }
        bool valid() const { return path->valid(); }
        bool empty() const { return path->empty(); }
        bool isSingleton() const {return path->isSingleton();}
        bool isRC() const { return rc; }

        IterableStorage<PathVertexIterator> vertices();

        PathIterator begin() const;
        PathIterator end() const;

        ConstPathDirection RC() const { return {*path, !rc}; }

        bool operator==(const ConstPathDirection &other) const { return path == other.path && rc == other.rc; }
        bool operator!=(const ConstPathDirection &other) const { return !(*this == other); }
    };

    //TODO: create const versions for all iterators and other containers
//    This iterator survives extending the path to left or right. It also survives cutting assuming that the position was not cut.
    class PathIterator {
    private:
        PathPosition position;
    public:
        PathIterator(PathPosition position) : position(position) {} // NOLINT(google-explicit-constructor)

        typedef Edge &reference;
        typedef Edge *pointer;
        typedef Edge value_type;

        reference operator*() const { return position.nextEdge(); }

        pointer operator->() const { return &position.nextEdge(); }

        PathIterator &operator++() {
            ++position;
            return *this;
        }

        PathIterator &operator--() {
            --position;
            return *this;
        }

        PathIterator operator+(int d) const {
            return {position + d};
        }

        PathIterator operator-(int d) const { return operator+(-d); }

        PathIterator operator++(int) { PathIterator res = *this; ++res; return res; }

        PathIterator operator--(int) { PathIterator res = *this; --res; return res; }

//        Careful! RC does not point to the rc edge. Instead it makes sure that begin->end and end->begin
        PathIterator RC() const { return {position.RC()}; }

        PathIterator SameElementRC() const { return RC() - 1; }
        bool operator==(const PathIterator &other) const {return position == other.position;}
        bool operator!=(const PathIterator &other) const {return !(*this == other);}
    };

    class SegmentIterator {
    private:
        ConstPathDirection direction;
        PathPosition position;
    public:
        using iterator_category = std::forward_iterator_tag;
        using reference = Segment<Edge>;
        using value_type = Segment<Edge>;
        SegmentIterator(const ConstPathDirection &direction, PathPosition position) :
                direction(direction), position(position) {
        }


        reference operator*() const {
            return direction.getSegment(position);
        }

        SegmentIterator &operator++() {
            ++position;
            return *this;
        }

        SegmentIterator operator++(int) { SegmentIterator res = *this; ++res; return res; }
        bool operator==(const SegmentIterator &other) const {return direction == other.direction && position == other.position;}
        bool operator!=(const SegmentIterator &other) const {return !(*this == other);}
    };


//    Since we can not point to a vertex after the last, we have no natural end iterator.
//    Instead end iterator corresponds to fake position with invalid vid and fpos, rpos, pointing to last position.
    class PathVertexIterator {
    public:
        typedef Vertex &reference;
        typedef Vertex *pointer;
    private:
        ConstPathDirection dir;
        PathPosition position;
    public:
        PathVertexIterator(const ConstPathDirection &dir, PathPosition position) : dir(dir), position(position) {}

        reference operator*() const { return position.getVertex(); }
        pointer operator->() const { return &position.getVertex(); }
        PathVertexIterator &operator++();
        PathVertexIterator operator++(int);

        bool operator==(const PathVertexIterator &other) const;
        bool operator!=(const PathVertexIterator &other) const;
    };

    template<class Iterator>
    GraphPath::GraphPath(Iterator begin, Iterator end) : start({}), rc_start({}), cut_left(0), cut_right(0) {
        while (begin != end) {
            *this += **begin;
            ++begin;
        }
    }
}