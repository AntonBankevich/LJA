#pragma once
#include <utility>

#include "assembly_graph_base.hpp"

namespace ag {
    template<class Traits>
    class PathDirection;
    template<class Traits>
    class ConstPathDirection;
    template<class Traits>
    class PathIterator;
    template<class Traits>
    class PathVertexIterator;
    template<class Traits>
    class SegmentIterator;
    template<class Traits>
    class PathPosition;

//    Invariants for normal paths. They do not have to be always fulfilled but need to keep track of them.
//    1. If edge in the path is neither first, nor last, its full sequence is inside the sequence.
//    2. No prefix edges on start.
//    3. No suffix edges on end.
//    4. In supregraph any path either contains core vertex or consists of a single vertex,
//    which is one of the minimal (by inclusion) vertices, containing the sequence.
//    5. In DBG any path is the shortest path containing the sequence
    template<class Traits>
    class GraphPath {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Edge::EdgeId EdgeId;

//        friend class PathIterator<Traits><Traits>;
        friend class PathDirection<Traits>;
        friend class PathVertexIterator<Traits>;
        friend class PathPosition<Traits>;

    private:
        VertexId start = {}; //The first vertex in the path
        VertexId rc_start = {}; //The rc to the last vertex in the path
        NuclDeck fsplits; //Stores concatenated codes for all non-suffix edges
        NuclDeck rsplits; //Stores concatenated codes for all non-suffix edges in rc path
//        In DBG every edge of the path is recorded in both fsplits and rsplits.
//        In SPG every edge is in exactly one of them
        size_t cut_left = 0;
        size_t cut_right = 0;

        GraphPath<Traits>(Vertex &start, Vertex &rc_start, NuclDeck fsplits, NuclDeck rsplits, size_t cut_left = 0, size_t cut_right = 0) :
                start(start.getId()), rc_start(rc_start.getId()), fsplits(std::move(fsplits)), rsplits(std::move(rsplits)),
                cut_left(cut_left), cut_right(cut_right) {}

    public:
        GraphPath<Traits>(const GraphPath<Traits> &) = default;
        GraphPath<Traits>(GraphPath<Traits> &&other)  noexcept {
            *this = std::move(other);
        }
        GraphPath<Traits> &operator=(const GraphPath<Traits> &) = default;
        GraphPath<Traits> &operator=(GraphPath<Traits> &&other) noexcept  = default;
        explicit GraphPath<Traits>(std::vector<EdgeId> path, size_t cut_left = 0, size_t cut_right = 0);

//        typedef typename std::vector<EdgeId>::iterator iterator;
//        typedef typename std::vector<EdgeId>::const_iterator const_iterator;
//        typedef TransformingIterator<CountingIterator<size_t>, Vertex> vertex_iterator;
//        typedef TransformingIterator<CountingIterator<size_t>, Edge> edge_iterator;
//        typedef Generator<CountingIterator<size_t>, Segment<Edge>> segment_iterator;


        GraphPath(Vertex &start, size_t cut_left = 0, size_t cut_right = 0) : start(start.getId()), rc_start(start.rc().getId()), // NOLINT(google-explicit-constructor)
                                            cut_left(cut_left), cut_right(cut_right) {
            VERIFY(this->start == start.getId());
        }
        GraphPath(Vertex &start, const Sequence &extension);
        GraphPath(Edge &edge, size_t cut_left = 0, size_t cut_right = 0); // NOLINT(google-explicit-constructor)
        GraphPath(const Segment<Edge> &segment) : GraphPath<Traits>(segment.contig(), segment.cutLeft(), segment.cutRight()) {} // NOLINT(google-explicit-constructor)

        GraphPath() : start({}), rc_start({}), fsplits(), rsplits(), cut_left(0), cut_right(0) {}

        template<class Iterator>
        explicit GraphPath<Traits>(Iterator begin, Iterator end);

        std::string str() const;
        std::vector<EdgeId> asEdgeIds() const;
        size_t calculateSize() const;

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
        Segment<Edge> getSegment(PathPosition<Traits> position) const;
        bool empty() const { return fsplits.empty() && rsplits.empty(); }
//        It is important to keep track of how this method behaves when graph is modified in parallel.
        bool isSingleton() const;
        Segment<Edge> back() const;
        Segment<Edge> front() const;
        PathPosition<Traits> firstPosition() const;
        PathPosition<Traits> lastPosition() const;
        PathPosition<Traits> endPosition() const;
        PathPosition<Traits> rcEndPosition() const;
        const NuclDeck &getFSplits() const {return fsplits;}
        const NuclDeck &getRSplits() const {return rsplits;}
        bool operator==(const GraphPath<Traits> &other) const;
        bool operator!=(const GraphPath<Traits> &other) const { return !operator==(other); }

        PathDirection<Traits> forward() {return {*this, false};}
        PathDirection<Traits> backward() {return {*this, true};}
        ConstPathDirection<Traits> forward() const {return {*this, false};}
        ConstPathDirection<Traits> backward() const {return {*this, true};}

        void pop_back();
        void pop_front();
        void replace_back(Edge &edge);
        void replace_front(Edge &edge);
        void normalize();
        void invalidate() {*this = {};}
        void setCutLeft(size_t value);
        void setCutRight(size_t value);
        void operator+=(const GraphPath<Traits> &other);
        void operator+=(const Segment<Edge> &other);
        void operator+=(Edge &other);
        void push_front(Edge &other);
        void push_front(const Segment<Edge> &other);
        void push_front(const GraphPath<Traits> &other);
        void pop_back(size_t len);
        void pop_front(size_t len);
        GraphPath<Traits> &shorten(PathPosition<Traits> from, PathPosition<Traits> to);
        GraphPath<Traits> &cutBack(size_t l);
        GraphPath<Traits> &cutFront(size_t l);
        ag::GraphPath<Traits> &addStep() {cut_right -= 1; return *this;}
        ag::GraphPath<Traits> &addStep(Edge &edge) {*this += Segment<Edge>(edge, 0, 1);return *this;}
        ag::GraphPath<Traits> & extend(const Sequence &seq);
        ag::GraphPath<Traits> &fastExtend(const Sequence &seq);
        void resetEdgeCodes();

        void forcePushBack(Edge &edge);
        void forcePushFront(Edge &edge);

        //        TODO: Find a way to iterate over temporary path objects
        IterableStorage<PathVertexIterator<Traits>> vertices() const &;
        IterableStorage<PathVertexIterator<Traits>> innerVertices() const &;
        IterableStorage<PathVertexIterator<Traits>> vertices() && = delete;
        IterableStorage<PathVertexIterator<Traits>> innerVertices() && = delete;
        IterableStorage<PathIterator<Traits>> edges() const &;
        IterableStorage<PathIterator<Traits>> edges() && = delete;
        SegmentIterator<Traits> begin() const;
        SegmentIterator<Traits> end() const;

//        TODO: minimaze usage of this function. Use directions instead.
        GraphPath<Traits> RC() const {return {*rc_start, *start, rsplits, fsplits, cut_right, cut_left};}
        GraphPath<Traits> subPath(PathPosition<Traits> from, PathPosition<Traits> to) const;
        GraphPath<Traits> subPath(PathPosition<Traits> from) const;
        GraphPath<Traits> operator+(const GraphPath<Traits> &other) const;
        GraphPath<Traits> operator+(const Segment<Edge> &other) const;
        GraphPath<Traits> operator+(Edge &other) const;
        GraphPath<Traits> operator*(size_t mult) const;
        Sequence Seq() const;
        Sequence truncSeq() const;
        size_t truncLen() const;
        size_t len() const;
        static GraphPath<Traits> Load(std::istream &os, const IdIndex<Vertex> &index);

//        These functions should return the same as if sequences of the paths were compared
        bool startsWith(const GraphPath<Traits> &other) const;
        bool endsWith(const GraphPath<Traits> &other) const;
        bool nonContradicts(const GraphPath<Traits> &other) const;
    };

    template<class Traits>
    inline std::ostream &operator<<(std::ostream &os, const GraphPath<Traits> &path) {
        if (path.valid())
            return os << path.getStart().getInnerId() << " " << path.getFinish().rc().getInnerId() << " F:" <<
                      path.getFSplits() << " R:" << path.getRSplits() << " " << path.leftCut() << " " << path.rightCut();
        else
            return os << "0 0 F: R: 0 0";
    }

    template<class Traits>
    class PathPosition {
        friend class GraphPath<Traits>;
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Traits::Edge Edge;
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
        PathPosition<Traits>(Vertex &cur, NuclDeck::Iterator fpos, NuclDeck::Iterator rpos) : vid(cur.getId()), fpos(fpos), rpos(rpos) {
        }
        PathPosition<Traits>(NuclDeck::Iterator fpos, NuclDeck::Iterator rpos) : vid({}), fpos(fpos), rpos(rpos) {
        }
        PathPosition<Traits>(const PathPosition<Traits> &) = default;
        PathPosition<Traits>(PathPosition<Traits> &&)  noexcept = default;
        PathPosition<Traits> &operator=(const PathPosition<Traits> &) = default;
        PathPosition<Traits> &operator=(PathPosition<Traits> &&) = default;

        Edge &nextEdge() const {return vid->getOutgoingByIterator(fpos);}
        Edge &prevEdge() const {return vid->rc().getOutgoingByIterator(rpos).rc();}
        Vertex &getVertex() const {return *vid;}
        NuclDeck::Iterator getFPos() const {return fpos;}
        NuclDeck::Iterator getRPos() const {return rpos;}


//        This method is for minor optimization. we avoid calling nextEdge() if we already know the next edge.

        PathPosition<Traits> &operator+=(Edge &edge) {
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

        PathPosition<Traits> &operator++() {
            return operator+=(nextEdge());
        }

        PathPosition<Traits> &operator--() {
            Edge &edge = prevEdge();
            fpos -= edge.getCode().size();
            rpos += edge.rc().getCode().size();
            vid = edge.getStart().getId();
            return *this;
        }

        PathPosition<Traits> operator+(int d) const {
            PathPosition<Traits> res = *this;
            res.move(d);
            return res;
        }

        PathPosition<Traits> operator-(int d) const {return *this + (-d);}
        PathPosition<Traits> RC() const {return {vid->rc(), rpos, fpos};}

        bool operator==(const PathPosition<Traits> &other) const {
            VERIFY(&fpos.getDeck()==&other.fpos.getDeck());
            VERIFY(&rpos.getDeck()==&other.rpos.getDeck());
            return fpos == other.fpos && rpos == other.rpos;
        }
        bool operator<=(const PathPosition<Traits> &other) const {return fpos <= other.fpos && rpos >= other.rpos;}
        bool operator!=(const PathPosition<Traits> &other) const { return !(*this == other); }
        bool operator<(const PathPosition<Traits> &other) const {return *this <= other && *this != other;}
    };

    template<class Traits>
    class PathDirection {
    private:
        GraphPath<Traits> *path;
        bool rc;
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Edge::EdgeId EdgeId;

        PathDirection<Traits>() : path(), rc(false) {}

        PathDirection<Traits>(GraphPath<Traits> &path, bool rc) : path(&path), rc(rc) {}
        PathDirection<Traits> &operator=(const GraphPath<Traits> &other);
        PathDirection<Traits> &operator=(GraphPath<Traits> &&other);

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
        Segment<Edge> getSegment(PathPosition<Traits> position) const;

        IterableStorage<PathVertexIterator<Traits>> vertices();

        void pop_front() const { return rc ? path->pop_back() : path->pop_front(); }
        void pop_back() const { return rc ? path->pop_front() : path->pop_back(); }
        void push_front(Edge &edge) const {if(rc) *path += edge.rc(); else path->push_front(edge);}
        void push_back(Edge &edge) const {if(rc) path->push_front(edge.rc()); else *path += edge;}
        void replace_front(Edge &edge) const {if(rc) path->replace_back(edge.rc()); else path->replace_front(edge);}
        void replace_back(Edge &edge) const {if(rc) path->replace_front(edge.rc()); else path->replace_back(edge);}
        void operator+=(Edge &edge) const { push_back(edge);}
        void operator+=(Segment<Edge> seg) const { if(rc) path->push_front(seg.RC()); else *path += seg;}
        void operator+=(const GraphPath<Traits> &other) const { if(rc) path->push_front(other.RC()); else *path += other;}
        void push_front(Segment<Edge> seg) const {if(rc) *path += seg.RC(); else path->push_front(seg);}
        void push_front(const GraphPath<Traits> &other) const { if(rc) *path += other.RC(); else path->push_front(other);}
        void invalidate() const { path->invalidate(); }
        void setCutLeft(size_t val) const;
        void setCutRight(size_t val) const;
        void forcePushBack(Edge &edge) const;
        void forcePushFront(Edge &edge) const;

        PathPosition<Traits> firstPosition() const {return rc ? path->lastPosition().RC() : path->firstPosition();}
        PathPosition<Traits> lastPosition() const {return rc ? path->firstPosition().RC() : path->lastPosition();}
        PathPosition<Traits> endPosition() const {return rc ? path->rcEndPosition() : path->endPosition();}
        PathPosition<Traits> rcEndPosition() const {return rc ? path->endPosition() : path->rcEndPosition();}

        bool isForward() const { return !rc; }
        size_t calculateSize() const { return path->calculateSize(); }
        GraphPath<Traits> &getPath() const { return *path; }
        bool valid() const { return path->valid(); }
        bool empty() const { return path->empty(); }
        bool isSingleton() const;
        bool isRC() const { return rc; }

        PathIterator<Traits> begin() const;
        PathIterator<Traits> end() const;

        PathDirection<Traits> RC() const { return {*path, !rc}; }

//        void forceReroute(PathIterator<Traits><Traits> from, PathIterator<Traits><Traits> to, const GraphPath<Traits><Traits> &alt) const;
        bool operator==(const PathDirection<Traits> &other) const { return path == other.path && rc == other.rc; }
        bool operator!=(const PathDirection<Traits> &other) const { return !(*this == other); }
    };

    template<class Traits>
    class ConstPathDirection {
    private:
        const GraphPath<Traits> *path;
        bool rc;

    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Edge::EdgeId EdgeId;

        ConstPathDirection<Traits>() : path(), rc(false) {}
        ConstPathDirection<Traits>(const PathDirection<Traits> &other) : path(&other.getPath()), rc(other.isRC()) {}

        ConstPathDirection<Traits>(const GraphPath<Traits> &path, bool rc) : path(&path), rc(rc) {}

        Vertex &getStart() const { return rc ? path->getFinish().rc() : path->getStart(); }
        Vertex &getFinish() const { return rc ? path->getStart().rc() : path->getFinish(); }
        Segment<Edge> front() const { return rc ? path->back().RC() : path->front(); }
        Segment<Edge> back() const { return rc ? path->front().RC() : path->back(); }
        Edge &frontEdge() const { return rc ? path->backEdge().rc() : path->frontEdge(); }
        Edge &backEdge() const { return rc ? path->frontEdge().rc() : path->backEdge(); }
        size_t cutRight() const { return rc ? path->leftCut() : path->rightCut(); }
        size_t cutLeft() const { return rc ? path->rightCut() : path->leftCut(); }
        PathPosition<Traits> firstPosition() const {return rc ? path->lastPosition().RC() : path->firstPosition();}
        PathPosition<Traits> lastPosition() const {return rc ? path->firstPosition().RC() : path->lastPosition();}
        PathPosition<Traits> endPosition() const {return rc ? path->rcEndPosition() : path->endPosition();}
        PathPosition<Traits> rcEndPosition() const {return rc ? path->endPosition() : path->rcEndPosition();}
        Segment<Edge> getSegment(PathPosition<Traits> position) const;
        bool isForward() const { return !rc; }
        size_t calculateSize() const { return path->calculateSize(); }
        const GraphPath<Traits> &getPath() const { return *path; }
        bool valid() const { return path->valid(); }
        bool empty() const { return path->empty(); }
        bool isSingleton() const {return path->isSingleton();}
        bool isRC() const { return rc; }

        IterableStorage<PathVertexIterator<Traits>> vertices() {return {{*this, firstPosition()}, {*this, endPosition()}};}

        PathIterator<Traits> begin() const {return {firstPosition()};}
        PathIterator<Traits> end() const {return {lastPosition()};}

        ConstPathDirection<Traits> RC() const { return {*path, !rc}; }

        bool operator==(const ConstPathDirection<Traits> &other) const { return path == other.path && rc == other.rc; }
        bool operator!=(const ConstPathDirection<Traits> &other) const { return !(*this == other); }
    };

    template<class Traits>
    Segment<typename Traits::Edge> ConstPathDirection<Traits>::getSegment(PathPosition<Traits> position) const {
        Edge &next = position.nextEdge();
        size_t cut_left = firstPosition() == position ? cutLeft() : 0;
        size_t cut_right = lastPosition() == position + 1 ? cutRight() : 0;
        return {next, cut_left, next.truncSize() - cut_right};
    }

    //TODO: create const versions for all iterators and other containers
//    This iterator survives extending the path to left or right. It also survives cutting assuming that the position was not cut.
    template<class Traits>
    class PathIterator {
    private:
        PathPosition<Traits> position;
    public:
        PathIterator<Traits>(PathPosition<Traits> position) : position(position) {} // NOLINT(google-explicit-constructor)

        typedef typename Traits::Edge Edge;
        typedef Edge &reference;
        typedef Edge *pointer;
        typedef Edge value_type;

        reference operator*() const { return position.nextEdge(); }

        pointer operator->() const { return &position.nextEdge(); }

        PathIterator<Traits> &operator++() {
            ++position;
            return *this;
        }

        PathIterator<Traits> &operator--() {
            --position;
            return *this;
        }

        PathIterator<Traits> operator+(int d) const {
            return {position + d};
        }

        PathIterator<Traits> operator-(int d) const { return operator+(-d); }

        PathIterator<Traits> operator++(int) { PathIterator<Traits> res = *this; ++res; return res; }

        PathIterator<Traits> operator--(int) { PathIterator<Traits> res = *this; --res; return res; }

//        Careful! RC does not point to the rc edge. Instead it makes sure that begin->end and end->begin
        PathIterator<Traits> RC() const { return {position.RC()}; }

        PathIterator<Traits> SameElementRC() const { return RC() - 1; }
        bool operator==(const PathIterator<Traits> &other) const {return position == other.position;}
        bool operator!=(const PathIterator<Traits> &other) const {return !(*this == other);}
    };

    template<class Traits>
    class SegmentIterator {
    private:
        ConstPathDirection<Traits> direction;
        PathPosition<Traits> position;
    public:
        typedef typename Traits::Edge Edge;
        using iterator_category = std::forward_iterator_tag;
        using reference = Segment<Edge>;
        using value_type = Segment<Edge>;
        SegmentIterator<Traits>(const ConstPathDirection<Traits> &direction, PathPosition<Traits> position) :
                direction(direction), position(position) {
        }


        reference operator*() const {
            return direction.getSegment(position);
        }

        SegmentIterator<Traits> &operator++() {
            ++position;
            return *this;
        }

        SegmentIterator<Traits> operator++(int) { SegmentIterator<Traits> res = *this; ++res; return res; }
        bool operator==(const SegmentIterator<Traits> &other) const {return direction == other.direction && position == other.position;}
        bool operator!=(const SegmentIterator<Traits> &other) const {return !(*this == other);}
    };


    //    Since we can not point to a vertex after the last, we have no natural end iterator.
//    Instead end iterator corresponds to fake position with invalid vid and fpos, rpos, pointing to last position.
    template<class Traits>
    class PathVertexIterator {
    private:
        const GraphPath<Traits> * path;
        PathPosition<Traits> position;
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        PathVertexIterator<Traits>(const GraphPath<Traits> &path, PathPosition<Traits> position) :
                path(&path), position(position) {
        }

//        typedef typename Traits::Edge Edge;
        typedef Vertex &reference;
        typedef Vertex *pointer;

        reference operator*() const { return position.getVertex(); }
        pointer operator->() const { return &position.getVertex(); }
        PathVertexIterator<Traits> &operator++() {
            if(position == path->lastPosition())
                position = path->endPosition();
            else
                ++position;
            return *this;
        }
        PathVertexIterator<Traits> operator++(int) { PathVertexIterator<Traits> res = *this; ++res; return std::move(res); }
        bool operator==(const PathVertexIterator<Traits> &other) const {return path == path && position == other.position;}
        bool operator!=(const PathVertexIterator<Traits> &other) const {return !(*this == other);}
    };

    template<class Traits>
    IterableStorage<PathVertexIterator<Traits>> GraphPath<Traits>::vertices() const &{
        if(valid())
            return {{*this, firstPosition()}, {*this, {endPosition()}}};
        else
            return {{*this, endPosition()}, {*this, {endPosition()}}};
    }

    template<class Traits>
    IterableStorage<PathVertexIterator<Traits>> GraphPath<Traits>::innerVertices() const & {
        if(empty())
            return {{*this, endPosition()},{*this, endPosition()}};
        return {{*this, firstPosition() + 1},{*this, lastPosition()}};
    }

    template<class Traits>
    IterableStorage<PathIterator<Traits>> GraphPath<Traits>::edges() const & {
        return {{firstPosition()}, {lastPosition()}};
    }

    template<class Traits>
    SegmentIterator<Traits> GraphPath<Traits>::begin() const {return {forward(), forward().firstPosition()};}

    template<class Traits>
    SegmentIterator<Traits> GraphPath<Traits>::end() const {return {forward(), forward().lastPosition()};}

    template<class Traits>
    Sequence GraphPath<Traits>::truncSeq() const {
        SequenceBuilder sb;
        PathPosition<Traits> first = firstPosition();
        PathPosition<Traits> last = lastPosition();
        for (PathPosition<Traits> pp = first; pp != last;) {
            Edge &edge = pp.nextEdge();
            size_t left = pp == first ? cut_left : 0;
            ++pp;
            size_t right = pp == last ? cut_right : 0;
            sb.append(edge.truncSeq().Subseq(left, edge.truncSize() - right));
        }
        return sb.BuildSequence();
    }

    template<class Traits>
    Sequence GraphPath<Traits>::Seq() const {
        if(empty()) {
            return getStart().getSeq().Subseq(cut_left, getStart().size() - cut_right);
        }
        SequenceBuilder sb;
        size_t first_right = isSingleton() ? cut_right : 0;
        if(cut_left >= getStart().size())
            sb.append(frontEdge().truncSeq().Subseq(cut_left - getStart().size(), frontEdge().truncSize() - first_right));
        else
            sb.append(getStart().getSeq().Subseq(cut_left) + frontEdge().truncSeq().Subseq(0, frontEdge().truncSize() - first_right));
        PathPosition<Traits> first = firstPosition();
        PathPosition<Traits> last = lastPosition();
        for (PathPosition<Traits> pp = first + 1; pp != last;) {
            Edge &edge = pp.nextEdge();
            ++pp;
            size_t right = pp == last ? cut_right : 0;
            sb.append(edge.truncSeq().Subseq(0, edge.truncSize() - right));
        }
        return sb.BuildSequence();
    }

    template<class Traits>
    size_t GraphPath<Traits>::truncLen() const {
        return valid() ? len() - start->size() : 0;
    }

    template<class Traits>
    size_t GraphPath<Traits>::len() const {
        size_t res = start->size();
        for (Edge &edge : edges())
            res += edge.truncSize();
        return res - cut_left - cut_right;
    }

    template<class Traits>
    Segment<typename Traits::Edge> GraphPath<Traits>::back() const {
        return {backEdge(), (isSingleton() ? leftCut() : 0), backEdge().truncSize() - rightCut()};
    }

    template<class Traits>
    Segment<typename Traits::Edge> GraphPath<Traits>::front() const {
        return {frontEdge(), leftCut(), isSingleton() ? frontEdge().truncSize() - rightCut() : frontEdge().truncSize()};
    }

//    template<class Traits>
//    Segment<typename Traits::Edge> GraphPath<Traits>::operator[](size_t i) const {
//        return {getEdge(i), i == 0 ? cutLeft() : 0, i == size() - 1 ? backEdge().truncSize() - cutRight() : getEdge(i).truncSize()};
//    }

    template<class Traits>
    PathPosition<Traits> GraphPath<Traits>::firstPosition() const {
        return valid() ? PathPosition<Traits>(getStart(), fsplits.begin(), rsplits.end()) : endPosition();
    }

    template<class Traits>
    PathPosition<Traits> GraphPath<Traits>::lastPosition() const {
        return valid() ? PathPosition<Traits>(getFinish(), fsplits.end(), rsplits.begin()) : endPosition();
    }

    template<class Traits>
    PathPosition<Traits> GraphPath<Traits>::endPosition() const {
        return {fsplits.end(), rsplits.begin()};
    }

    template<class Traits>
    PathPosition<Traits> GraphPath<Traits>::rcEndPosition() const {
        return {rsplits.end(), fsplits.begin()};
    }

    template<class Traits>
    GraphPath<Traits> GraphPath<Traits>::subPath(PathPosition<Traits> from) const { return subPath(from, lastPosition()); }

    template<class Traits>
    GraphPath<Traits> GraphPath<Traits>::subPath(PathPosition<Traits> from, PathPosition<Traits> to) const {
//        Handles the case when the path is a vertex segment
        if(from == firstPosition() && to == lastPosition())
            return {*this};
        if(!valid()) {
            VERIFY(from == endPosition());
            VERIFY(to == endPosition());
            return {};
        }
        GraphPath<Traits> res(from.getVertex());
        PathPosition<Traits> cur = from;
        while(cur != to) {
            res += cur.nextEdge();
            ++cur;
        }
        if(from == firstPosition() && to != firstPosition())
            res.setCutLeft(cut_left);
        if(to == lastPosition() && from != lastPosition())
            res.setCutRight(cut_right);
        return std::move(res);
    }

    template<class Traits>
    void GraphPath<Traits>::operator+=(const GraphPath<Traits> &other) {
        VERIFY(this != &other);
        if(!other.valid()) {
            return;
        } else if(!valid()) {
            *this = other;
        } else if(cut_right > 0 && !empty()) {
            Edge &old = backEdge();
            fsplits.replaceBack(old.getCode().size(), other.fsplits);
            rsplits.replaceFront(old.rc().getCode().size(), other.rsplits);
        } else {
            fsplits.push_back(other.fsplits);
            rsplits.push_front(other.rsplits);
        }
        rc_start = other.rc_start;
        cut_right = other.cut_right;
    }

    template<class Traits>
    void GraphPath<Traits>::operator+=(const Segment<Edge> &other) {
        if(!valid()) {
            *this = {other};
            return;
        }
        VERIFY((cut_right == 0 && (other.cutLeft() == 0 || empty())) || (!empty() && backEdge() == other.contig() && cut_right + other.left == backEdge().truncSize()));
        if(cut_right == 0) {
            VERIFY(cut_right == 0 && (other.cutLeft() == 0 || empty()));
            VERIFY(getFinish() == other.contig().getStart());
            *this += other.contig();
        }
        cut_right = other.cutRight();
    }

    template<class Traits>
    void GraphPath<Traits>::operator+=(Edge &other) {
        VERIFY(cut_right == 0);
        VERIFY(!valid() || getFinish() == other.getStart());
        if(!valid())
            start = other.getStart().getId();
        fsplits.push_back(other.getCode());
        rsplits.push_front(other.rc().getCode());
        rc_start = other.getFinish().rc().getId();
    }

    template<class Traits>
    GraphPath<Traits> GraphPath<Traits>::operator+(const GraphPath<Traits> &other) const {
        GraphPath<Traits> res = *this;
        res += other;
        return std::move(res);
    }

    template<class Traits>
    GraphPath<Traits> GraphPath<Traits>::operator+(const Segment<Edge> &other) const {
        GraphPath<Traits> res = *this;
        res += other;
        return std::move(res);
    }

    template<class Traits>
    GraphPath<Traits> GraphPath<Traits>::operator+(Edge &other) const {
        GraphPath<Traits> res = *this;
        res += other;
        return std::move(res);
    }

    template<class Traits>
    void GraphPath<Traits>::pop_back(size_t len) {
        for (size_t i = 0; i < len; i++)
            pop_back();
    }

    template<class Traits>
    void GraphPath<Traits>::pop_front(size_t len) {
        for (size_t i = 0; i < len; i++)
            pop_front();
    }

    template<class Traits>
    GraphPath<Traits> &GraphPath<Traits>::cutBack(size_t l) {
        VERIFY(l <= len());
        size_t expected = len() - l;
        l += cut_right;
        while(!empty() && l >= backEdge().truncSize()) {
            l -= backEdge().truncSize();
            pop_back();
        }
        cut_right = l;
        VERIFY(len() == expected);
        return *this;
    }

    template<class Traits>
    GraphPath<Traits> &GraphPath<Traits>::cutFront(size_t l) {
        VERIFY(l <= len());
        size_t expected = len() - l;
        l += cut_left;
        while(!empty() && l >= frontEdge().rc().truncSize()) {
            l -= frontEdge().rc().truncSize();
            pop_front();
        }
        cut_left = l;
        VERIFY(len() == expected);
        return *this;
    }

    template<class Traits>
    PathIterator<Traits> PathDirection<Traits>::begin() const {
        return {firstPosition()};
    }

    template<class Traits>
    PathIterator<Traits> PathDirection<Traits>::end() const {
        return {lastPosition()};
    }

    template<class Traits>
    IterableStorage<PathVertexIterator<Traits>> PathDirection<Traits>::vertices() {
        return {{*this, firstPosition()}, {*this, endPosition()}};
    }

    template<class Traits>
    void PathDirection<Traits>::setCutLeft(size_t val) const {
        if(rc)
            path->setCutRight(val);
        else
            path->setCutLeft(val);
    }

    template<class Traits>
    void PathDirection<Traits>::setCutRight(size_t val) const {
        if(rc)
            path->setCutLeft(val);
        else
            path->setCutRight(val);
    }

    template<class Traits>
    Segment<typename Traits::Edge> PathDirection<Traits>::getSegment(PathPosition<Traits> position) const {
        Edge &next = position.nextEdge();
        size_t cut_left = firstPosition() == position ? leftCut() : 0;
        size_t cut_right = lastPosition() == position + 1 ? rightCut() : 0;
        return {next, cut_left, next.truncSize() - cut_right};
    }

    template<class Traits>
    PathDirection<Traits> &PathDirection<Traits>::operator=(GraphPath<Traits> &&other) {
        if(rc)
            *path = other.RC();
        else
            *path = std::move(other);
    }

    template<class Traits>
    PathDirection<Traits> &PathDirection<Traits>::operator=(const GraphPath<Traits> &other) {
        if(rc)
            *path = other.RC();
        else
            *path = other;
    }

    template<class Traits>
    void PathDirection<Traits>::forcePushBack(Edge &edge) const {
        if(rc)
            path->forcePushFront(edge.rc());
        else
            path->forcePushBack(edge);
    }

    template<class Traits>
    void PathDirection<Traits>::forcePushFront(Edge &edge) const {
        if(rc)
            path->forcePushBack(edge.rc());
        else
            path->forcePushFront(edge);
    }

    template<class Traits>
    bool PathDirection<Traits>::isSingleton() const {
        return !empty() && getFSplits().size() == frontEdge().getCode().size() && getRSplits().size() == frontEdge().rc().getCode().size();
    }

    template<class Traits>
    std::string GraphPath<Traits>::str() const {
        if (!valid())
            return "";
        std::stringstream ss;
        ss << leftCut() << "[" << getStart().getInnerId() << "(" << getStart().size() << ")";
        for (const Edge &edge: edges()) {
            ss << "->" << edge.rc().getCode() << edge.rc().truncSize() << "(" << edge.getInnerId().eid << "|" << edge.getCoverage()<< "|" <<
                    edge.rc().getInnerId().eid << ")" << edge.getCode() << edge.truncSize() << "->" << edge.getFinish().getInnerId() << "(" <<
                    edge.getFinish().size() << ")";
        }
        ss << "]" << rightCut();
        return ss.str();
    }

    template<class Traits>
    Segment<typename Traits::Edge> GraphPath<Traits>::getSegment(PathPosition<Traits> position) const {
        Edge &next = position.nextEdge();
        size_t seg_cut_left = firstPosition() == position ? leftCut() : 0;
        size_t seg_cut_right = lastPosition() == position + 1 ? rightCut() : 0;
        return {next, seg_cut_left, next.truncSize() - seg_cut_right};
    }

    template<class Traits>
    void GraphPath<Traits>::push_front(Edge &other) {
        VERIFY(cut_left == 0);
        VERIFY(!valid() || getStart() == other.getFinish());
        if(!valid())
            rc_start = other.getFinish().rc().getId();
        fsplits.push_front(other.getCode());
        rsplits.push_back(other.rc().getCode());
        start = other.getStart().getId();
    }

    template<class Traits>
    void GraphPath<Traits>::push_front(const Segment<Edge> &other) {
        if(!valid()) {
            *this = {other};
            return;
        }
        VERIFY((cut_left == 0 && (other.cutRight() == 0 || empty())) || (!empty() && frontEdge() == other.contig() && cut_left == other.right));
        if(cut_left == 0) {
            push_front(other.contig());
        }
        cut_left = other.cutLeft();
    }

    template<class Traits>
    void GraphPath<Traits>::push_front(const GraphPath<Traits> &other) {
        if(!valid()) {
            *this = other;
        } else if(cut_left > 0 && !empty()) {
            Edge &old = frontEdge();
            fsplits.replaceFront(old.getCode().size(), other.fsplits);
            rsplits.replaceBack(old.rc().getCode().size(), other.rsplits);
        } else {
            fsplits.push_front(other.fsplits);
            rsplits.push_back(other.rsplits);
        }
        start = other.start;
        cut_left = other.cut_left;
    }

    template<class Traits>
    void GraphPath<Traits>::normalize() {
        if(empty())
            return;
        while(frontEdge().isPrefix())
            pop_front();
        while(backEdge().isSuffix())
            pop_back();
    }

    template<class Traits>
    void GraphPath<Traits>::pop_back() {
        Edge &back = backEdge();
        fsplits.pop_back(back.getCode().size());
        rsplits.pop_front(back.rc().getCode().size());
        size_t lost_len = back.truncSize();
        cut_right -= std::min(cut_right, lost_len);
        rc_start = back.getStart().rc().getId();
    }

    template<class Traits>
    void GraphPath<Traits>::pop_front() {
        Edge &front = frontEdge();
        fsplits.pop_front(front.getCode().size());
        rsplits.pop_back(front.rc().getCode().size());
        size_t lost_len = front.rc().truncSize();
        cut_left -= std::min(cut_left, lost_len);
        start = front.getFinish().getId();
    }

    template<class Traits>
    template<class Iterator>
    GraphPath<Traits>::GraphPath(Iterator begin, Iterator end) : start({}), rc_start({}), cut_left(0), cut_right(0) {
        while (begin != end) {
            *this += **begin;
            ++begin;
        }
    }

    template<class Traits>
    void GraphPath<Traits>::setCutLeft(size_t value) {
        cut_left = value;
    }

    template<class Traits>
    void GraphPath<Traits>::setCutRight(size_t value) {
        cut_right = value;
    }

    template<class Traits>
    GraphPath<Traits> &GraphPath<Traits>::extend(const Sequence &seq) {
        VERIFY(valid());
        for (size_t cpos = 0; cpos < seq.size(); cpos++) {
            unsigned char c = seq[cpos];
            if (endClosed()) {
                Vertex &v = getFinish();
                if (v.hasOutgoing(c)) {
                    Edge &edge = v.getOutgoing(c);
                    addStep(edge);
                } else {
                    invalidate();
                    return *this;
                }
            } else {
                if (backEdge().truncSeq()[backEdge().truncSize() - rightCut()] == c) {
                    addStep();
                } else {
                    invalidate();
                    return *this;
                }
            }
        }

        return *this;
    }

    template<class Traits>
    GraphPath<Traits> &GraphPath<Traits>::fastExtend(const Sequence &seq) {
        size_t pos = rightCut();
        cut_right = 0;
        while (pos < seq.size()) {
            Edge &e = getFinish().getOutgoing(seq[pos]);
            *this += e;
            pos += e.truncSize();
        }
        cut_right += pos - seq.size();
        return *this;
    }

    template<class Traits>
    bool GraphPath<Traits>::operator==(const GraphPath<Traits> &other) const {
        return start == other.start && rc_start == other.rc_start &&
               cut_left == other.cut_left && cut_right == other.cut_right &&
               fsplits == other.fsplits && rsplits == other.rsplits;
    }

    template<class Traits>
    GraphPath<Traits>::GraphPath(std::vector<EdgeId> path, size_t cut_left, size_t cut_right) : GraphPath() { // NOLINT(google-explicit-constructor)
        for(EdgeId eid : path) {
            *this += *eid;
        }
        setCutLeft(cut_left);
        setCutRight(cut_right);
    }

    template<class Traits>
    GraphPath<Traits>::GraphPath(Vertex &start, const Sequence &extension) : GraphPath(start) {
        size_t pos = 0;
        while(pos < extension.size()) {
            Edge &next = getFinish().getOutgoing(extension[pos]);
            *this += next;
            pos += next.getCode().size();
        }
    }

    template<class Traits>
    GraphPath<Traits> GraphPath<Traits>::Load(std::istream &os, const IdIndex<Vertex> &index) {
        typename Vertex::id_type startId;
        typename Vertex::id_type rcStartId;
        size_t left = 0;
        size_t right = 0;
        std::string fpath;
        std::string rpath;
        os >> startId >> rcStartId >> fpath >> rpath >> left >> right;
        fpath = fpath.substr(2);
        rpath = rpath.substr(2);
        if (startId == 0 && fpath.empty()) {
            return {};
        }
        return {index.getById(startId), index.getById(rcStartId), NuclDeck(fpath), NuclDeck(rpath), left, right};
    }

    template<class Traits>
    bool GraphPath<Traits>::startsWith(const GraphPath<Traits> &other) const {
        if(other.getStart() != getStart())
            return false;
        if(other.leftCut() != leftCut())
            return false;
        if(!fsplits.startsWith(other.fsplits) || !rsplits.endsWith(other.rsplits))
            return false;
        if(fsplits.size() == other.fsplits.size() && rsplits.size() == other.rsplits.size() && rightCut() > other.rightCut())
            return false;
        return true;
    }

    template<class Traits>
    bool GraphPath<Traits>::endsWith(const GraphPath<Traits> &other) const {
        if(other.getFinish() != getFinish())
            return false;
        if(other.rightCut() != rightCut())
            return false;
        if(!fsplits.endsWith(other.fsplits) || !rsplits.startsWith(other.rsplits))
            return false;
        if(fsplits.size() == other.fsplits.size() && rsplits.size() == other.rsplits.size() && leftCut() > other.leftCut())
            return false;
        return true;
    }

    template<class Traits>
    bool GraphPath<Traits>::nonContradicts(const GraphPath<Traits> &other) const {
        return startsWith(other) || other.startsWith(*this);
    }

    template<class Traits>
    void GraphPath<Traits>::forcePushBack(Edge &edge) {
        if(empty()) {
            start= edge.getStart().getId();
        }
        fsplits.push_back(edge.getCode());
        rsplits.push_front(edge.rc().getCode());
        rc_start = edge.getFinish().rc().getId();
    }

    template<class Traits>
    void GraphPath<Traits>::forcePushFront(Edge &edge) {
        if(empty()) {
            rc_start = edge.getFinish().rc().getId();
        }
        fsplits.push_front(edge.getCode());
        rsplits.push_back(edge.rc().getCode());
        start = edge.getStart().getId();
    }

    template<class Traits>
    size_t GraphPath<Traits>::calculateSize() const {
        size_t res = 0;
        for(Edge &edge : edges()) {
            res++;
        }
        return res;
    }

    template<class Traits>
    void GraphPath<Traits>::resetEdgeCodes() {
        if(!valid())
            return;
        NuclDeck newf;
        NuclDeck newr;
        for(Edge &edge : edges()) {
            newf.push_back(edge.firstNucl());
            newr.push_front(edge.rc().firstNucl());
        }
        fsplits = newf;
        rsplits = newr;
    }

    template<class Traits>
    GraphPath<Traits> GraphPath<Traits>::operator*(size_t mult) const {
        VERIFY(getStart() == getFinish());
        GraphPath<Traits> res(getStart());
        for(size_t i = 0; i < mult; i++)
            res += *this;
        return std::move(res);
    }

    template<class Traits>
    std::vector<typename Traits::Edge::EdgeId> GraphPath<Traits>::asEdgeIds() const {
        std::vector<EdgeId> res;
        for(Edge &edge : edges())
            res.template emplace_back(edge.getId());
        return std::move(res);
    }

    template<class Traits>
    void GraphPath<Traits>::replace_back(Edge &edge) {
        Edge &old = backEdge();
        fsplits.replaceBack(old.getCode().size(), edge.getCode());
        rsplits.replaceFront(old.rc().getCode().size(), edge.rc().getCode());
        rc_start = edge.getFinish().rc().getId();
    }

    template<class Traits>
    void GraphPath<Traits>::replace_front(Edge &edge) {
        Edge &old = frontEdge();
        fsplits.replaceFront(old.getCode().size(), edge.getCode());
        rsplits.replaceBack(old.rc().getCode().size(), edge.rc().getCode());
        start = edge.getStart().getId();
    }

    template<class Traits>
    bool GraphPath<Traits>::isSingleton() const {
        return !empty() && fsplits.size() == frontEdge().getCode().size() && rsplits.size() == frontEdge().rc().getCode().size();
    }

//    TODO: optimize to O(1)
    template<class Traits>
    GraphPath<Traits> &GraphPath<Traits>::shorten(PathPosition<Traits> from, PathPosition<Traits> to) {
        while(from != firstPosition())
            pop_front();
        while(to != lastPosition())
            pop_back();
        return *this;
    }

    template<class Traits>
    GraphPath<Traits>::GraphPath(Edge &edge, size_t cut_left, size_t cut_right) : GraphPath<Traits>() { // NOLINT(google-explicit-constructor)
        *this += edge;
        this->cut_left = cut_left;
        this->cut_right = cut_right;
    }

}