#pragma once

#include "assembly_graph/assembly_graph.hpp"
#include "sequences/contigs.hpp"

namespace ag {
    template<class Traits>
    class PathDirection;
    template<class Traits>
    class PathIterator;

    template<class Traits>
    class GraphPath {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Vertex::VertexId VertexId;
        typedef typename Edge::EdgeId EdgeId;
        friend class PathIterator<Traits>;
        friend class PathDirection<Traits>;
    private:
        VertexId start_;
        std::vector<EdgeId> path;
        size_t skip_left = 0;
        size_t skip_right = 0;
        size_t cut_left;
        size_t cut_right;

        void set(Edge &edge, size_t pos) {path[pos + skip_left] = edge.getId();}
    public:
        typedef typename std::vector<EdgeId>::iterator iterator;
        typedef typename std::vector<EdgeId>::const_iterator const_iterator;
        typedef TransformingIterator<CountingIterator<size_t>, Vertex> vertex_iterator;
        typedef TransformingIterator<CountingIterator<size_t>, Edge> edge_iterator;
        typedef Generator<CountingIterator<size_t>, Segment<Edge>> segment_iterator;

        GraphPath(Vertex &_start, std::vector<EdgeId> _path, size_t cut_left, size_t cut_right, size_t skip_left = 0, size_t skip_right = 0) :
                    start_(_start.getId()), path(std::move(_path)), cut_left(cut_left), cut_right(cut_right), skip_left(skip_left), skip_right(skip_right) {}
        GraphPath(Vertex &_start, size_t cut_left = 0, size_t cut_right = 0) : start_(_start.getId()), cut_left(cut_left), cut_right(cut_right) {} // NOLINT(google-explicit-constructor)
        GraphPath(Edge &edge) : start_(edge.getStart().getId()), path({edge.getId()}), cut_left(0), cut_right(0) {} // NOLINT(google-explicit-constructor)
        GraphPath(const Segment<Edge> &segment) : start_(segment.contig().getStart().getId()), // NOLINT(google-explicit-constructor)
                          path({segment.contig().getId()}), cut_left(segment.left), cut_right(segment.contig().truncSize() - segment.right) {}
        GraphPath() : start_({}), cut_left(0), cut_right(0) {}
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

        Vertex &getVertex(size_t i) const;
        Edge &getEdge(size_t i) const {
            if(i >= size()) {
                std::cout << i << " " << path.size() << " " << skip_left << " " << skip_right << std::endl;
            }
            VERIFY(i < size()); return *path[skip_left + i];
        }
        Vertex &start() const {return *start_;}
        Vertex &finish() const {return empty() ? start() : backEdge().getFinish();}
        size_t find(Edge &edge, size_t pos = 0) const;
        size_t find(Vertex &v, size_t pos = 0) const;
        Edge &backEdge() const {return *path[path.size() - skip_right - 1];}
        Edge &frontEdge() const {return *path[skip_left];}
        size_t size() const {return path.size() - skip_left - skip_right;}
        bool empty() const {return size() == 0;}

    //        TODO: Find a way to iterate over temporary path objects
        IterableStorage<vertex_iterator> vertices() const &;
        IterableStorage<vertex_iterator> vertices() && = delete;
        IterableStorage<edge_iterator> edges() const &;
        IterableStorage<edge_iterator> edges() && = delete;

        segment_iterator begin() const;
        segment_iterator  end() const;

        GraphPath RC() const;
        double minCoverage() const;
        Sequence Seq() const;
        Sequence truncSeq() const;
        size_t truncLen() const;
        size_t len() const;

        Segment<Edge> back() const;
        Segment<Edge> front() const;
        Segment<Edge> operator[](size_t i) const;
        std::string covStr(bool show_coverage = false) const;
        std::string str() const;
        std::string lenStr() const;

        bool valid() const;
        void invalidate();
        GraphPath subPath(size_t from, size_t to) const;
        GraphPath subPath(size_t from) const {return subPath(from, size());}
        GraphPath reroute(size_t left, size_t right, const GraphPath &rerouting) const;
        void simpleReroute(size_t left, size_t right, const std::vector<EdgeId> &edges);
        void operator+=(const GraphPath &other);
        void operator+=(const Segment<Edge> &other);
        void operator+=(Edge &other);
        GraphPath operator+(const GraphPath &other) const;
        GraphPath operator+(const Segment<Edge> &other) const;
        GraphPath operator+(Edge &other) const;
        //TODO deprecate
        Sequence map(std::unordered_map<const Edge *, Sequence> &edge_map);

        void pop_back() {
            skip_right++;
            VERIFY(skip_left + skip_right <= path.size());
            size_t new_cut_size = path[path.size() - skip_right]->truncSize();
            path[path.size() - skip_right] = {};
            cut_right -= std::min(cut_right, new_cut_size);
        }
        void pop_back(size_t len) {
            for(size_t i = 0; i < len; i++)
                pop_back();
        }
        void pop_front() {
            skip_left++;
            VERIFY(skip_left + skip_right <= path.size());
            size_t new_cut_size = path[skip_left - 1]->rc().truncSize();
            start_ = path[skip_left - 1]->getFinish().getId();
            path[skip_left - 1] = {};
            cut_left -= std::min(new_cut_size, cut_left);
        }
        void pop_front(size_t len) {
            for(size_t i = 0; i < len; i++)
                pop_front();
        }
        GraphPath &cutBack(size_t l);
        GraphPath &cutFront(size_t l);
        GraphPath &addStep();
        GraphPath &addStep(Edge &edge);
        std::vector<GraphPath> allSteps();
        std::vector<GraphPath> allExtensions(size_t len);
        GraphPath &extend(const Sequence &seq);

        GraphPath &fastExtend(const Sequence &seq);

        size_t cutLeft() const;
        size_t cutRight() const;
        bool endClosed() const;
        bool startClosed() const;

        Sequence truncSubseq(size_t start_position, size_t size = 10000000) const;

        bool operator==(const GraphPath &other) const {
            return start_ == other.start_ && cut_left == other.cut_left && cut_right == other.cut_right && size() == other.size() &&
                        std::equal(path.begin() + skip_left, path.end() - skip_right, other.path.begin() + other.skip_left);
        }
        bool operator!=(const GraphPath &other) const {return !operator==(other);}
    };

    template<class Traits>
    class PathIterator {
    private:
        GraphPath<Traits> *path;
        bool rc;
        int pos;
    public:
        PathIterator(GraphPath<Traits> &path, bool rc, int pos) : path(&path), rc(rc), pos(pos) {
        }
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Edge &reference;
        typedef typename Traits::Edge *pointer;

        reference operator*() const {return rc ? path->path[pos]->rc() : *path->path[pos];}
        pointer operator->() const {return rc ? &path->path[pos]->rc() : &(*path->path[pos]);}
        void set(Edge &edge) const {path->set(rc ? edge.rc() : edge, pos - path->skip_left);}
        PathIterator &operator++() {VERIFY(pos >= -1); rc ? pos-- : pos++; return *this;}
        PathIterator &operator--() {VERIFY(pos >= -1); rc ? pos++ : pos--; return *this;}
        PathIterator operator+(int d) const {return {*path, rc, rc ? pos - d : pos + d};}
        PathIterator operator-(int d) const { return operator+(-d); }
        PathIterator operator++(int) {return *this + 1;}
        PathIterator operator--(int) {return *this - 1;}
        bool operator==(const PathIterator &other) const { return rc == other.rc && path == other.path && pos == other.pos; }
        bool operator!=(const PathIterator &other) const { return !(*this == other); }
        std::string str() const {
            std::stringstream ss;
            ss << path->lenStr() << ":" << pos << "(" << (rc ? "B" : "F") << ")";
            return ss.str();
        }
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

        PathDirection(GraphPath<Traits> &path, bool rc) : path(&path), rc(rc) {}
        Vertex &getgetStart() const { return rc ? path->finish().rc() : path->start(); }
        Vertex &getgetFinish() const { return rc ? path->start().rc() : path->finish(); }
        size_t cutRight() const {return rc ? path->cutLeft() : path->cutRight();}
        size_t cutLeft() const {return rc ? path->cutRight() : path->cutLeft();}
        void pop_front() const {return rc ? path->pop_back() : path->pop_front();}
        void pop_back() const {return rc ? path->pop_front() : path->pop_back();}
        bool isForward() const {return !rc;}
        size_t size() const {return path->size();}
        GraphPath<Traits> &getPath() {return *path;}
        PathIterator<Traits> end() const {
            return {*path, rc, int(rc ? path->skip_left - 1: path->path.size() - path->skip_right)};
        }
        PathIterator<Traits> begin() const {
            return {*path, rc, int(rc ? path->path.size() - path->skip_right - 1: path->skip_left)};
        }
        PathDirection RC() const { return {*path, !rc}; }
        bool valid() const { return path->valid(); }
        bool empty() const {return path->empty();}

        void invalidate() {
            path->invalidate();
        }
        void rerouteSameSize(PathIterator<Traits> from, PathIterator<Traits> to, const std::vector<EdgeId> &alt) {
            size_t cur = 0;
            while(from != to && cur < alt.size()) {
                from.set(*alt[cur]);
                cur++;
                ++from;
            }
            VERIFY(from == to);
            VERIFY(cur == alt.size());
            if(!empty())
                path->start_ = path->frontEdge().getStart().getId();
        }
        void cutFrontPrefix() {
            VERIFY(!path->empty());
            if (!rc) {
                path->start = path->frontEdge().getFinish();
                path->path.pop_front();
            } else {
                path->path.pop_back();
            }
        }

        void cutBackSuffix() {
            RC().cutFrontPrefix();
        }
    };
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::WalkForward(Edge &start) {
    GraphPath<Graph> res(start);
    Vertex *next = &start.getFinish();
    VERIFY(next != nullptr);
    while(*next != start.getStart() && *next != start.getStart().rc() && !next->isJunction()) {
        VERIFY(next != nullptr);
        res += next->front();
        next = &res.finish();
    }
    return std::move(res);
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::subPath(size_t from, size_t to) const {
    if(!valid()) {
        VERIFY(from == 0 && to == 0);
        return {};
    }
    if (from == to) {
        if ((from == 0 && cutLeft() > 0) || (to == size() && cutRight() > 0)) {
            return {};
        } else {
            return GraphPath<Graph>(getVertex(from));
        }
    } else {
        return {getVertex(from), std::vector<EdgeId>(path.begin() + skip_left + from, path.begin() + skip_left + to),
                from == 0 ? cutLeft() : 0, to == size() ? cutRight() : 0};
    }
}

template<class Graph>
typename Graph::Vertex &ag::GraphPath<Graph>::getVertex(size_t i) const {
    VERIFY(valid());
    VERIFY(i <= size());
    if (i == 0)
        return *start_;
    else
        return path[skip_left + i - 1]->getFinish();
}

template<class Graph>
size_t ag::GraphPath<Graph>::find(Edge &edge, size_t pos) const {
    while(pos < size() && edge != getEdge(pos))
        pos++;
    if(pos == size())
        return -1;
    return pos;
}

template<class Graph>
size_t ag::GraphPath<Graph>::find(Vertex &v, size_t pos) const {
    while(pos <= size() && v != getVertex(pos))
        pos++;
    if(pos > size())
        return -1;
    return pos;
}

template<class Graph>
double ag::GraphPath<Graph>::minCoverage() const {
    double res = 100000000;
    for (const Edge &edge : edges()) {
        res = std::min(edge->getCoverage(), res);
    }
    return res;
}

template<class Graph>
size_t ag::GraphPath<Graph>::truncLen() const {
    size_t res = 0;
    for (Edge &edge : edges())
        res += edge.truncSize();
    if(res == 0)
        return 0;
    return res - cut_left - cut_right;
}

template<class Graph>
IterableStorage<typename ag::GraphPath<Graph>::vertex_iterator> ag::GraphPath<Graph>::vertices() const & {
    std::function<Vertex &(size_t)> transformer = [this](size_t ind)->Vertex &{return getVertex(ind);};
    CountingIterator<size_t> end_it = CountingIterator<size_t>(valid()? path.size() - skip_right + 1 : 0);
    vertex_iterator vbegin (CountingIterator<size_t>(skip_left), end_it, transformer);
    vertex_iterator vend(end_it, end_it, transformer);
    return {vbegin, vend};
}

template<class Graph>
IterableStorage<typename ag::GraphPath<Graph>::edge_iterator> ag::GraphPath<Graph>::edges() const &{
    std::function<Edge &(size_t)> transformer = [this](size_t ind)->Edge &{return *path[ind];};
    CountingIterator<size_t> end_it = CountingIterator<size_t>(path.size() - skip_right);
    edge_iterator ebegin (CountingIterator<size_t>(skip_left), end_it, transformer);
    edge_iterator eend(end_it, end_it, transformer);
    return {ebegin, eend};
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::RC() const {
    if(!valid())
        return {};
    std::vector<EdgeId> res;
    for(auto it  = path.rbegin(); it != path.rend(); ++it) {
        res.emplace_back((*it)->rc().getId());
    }
    return {finish().rc(), res, cutRight(), cutLeft(), skip_right, skip_left};
}

template<class Graph>
void ag::GraphPath<Graph>::invalidate() {
    start_ = {};
    path.clear();
    cut_left = 0;
    cut_right = 0;
    skip_left = 0;
    skip_right = 0;
}

template<class Graph>
bool ag::GraphPath<Graph>::valid() const {
    VERIFY(start_.valid() || size() == 0);
    return start_.valid();
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::cutBack(size_t l) {
    VERIFY(l <= len());
    size_t expected = len() - l;
    size_t cur_cut = 0;
    size_t cut = 0;
    l += cut_right;
    cut_right = 0;
    while(cur_cut < size() && l >= getEdge(size() - 1 - cur_cut).truncSize()) {
        if(getEdge(size() - 1 - cur_cut).truncSize() == 0) {
            cur_cut++;
        } else {
            l -= getEdge(size() - 1 - cur_cut).truncSize();
            cur_cut++;
            cut = cur_cut;
        }
    }
    pop_back(cut);
    cut_right = l;
    VERIFY(len() == expected);
    return *this;
}

//template<class Graph>
//ag::GraphPath<Graph> &ag::GraphPath<Graph>::uniqueExtendBack(size_t l) {
//    if(cut_right != 0) {
//        size_t tmp = std::min(l, cut_right);
//        l -= tmp;
//        cut_right -= tmp;
//    }
//    while(l > 0) {
//        VERIFY(finish().outDeg() == 1);
//        Edge &e = finish().front();
//        size_t tmp = std::min(e.truncSize(), l);
//        *this += e;
//        cutBack(e.truncSize() - tmp);
//    }
//    return *this;
//}


////TODO: Optimize
//template<class Graph>
//ag::GraphPath<Graph> &ag::GraphPath<Graph>::uniqueExtendFront(size_t l) {
//    *this = this->RC().uniqueExtendBack(l);
//    return *this;
//}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::cutFront(size_t l) {
    VERIFY(l <= len());
    size_t expected = len() - l;
    size_t cur_cut = 0;
    size_t cut = 0;
    l += cut_left;
    cut_left = 0;
    while(cur_cut < size() && l >= getEdge(cur_cut).rc().truncSize()) {
        if(getEdge(cur_cut).rc().truncSize() == 0) {
            cur_cut++;
        } else {
            l -= getEdge(cur_cut).rc().truncSize();
            cur_cut++;
            cut = cur_cut;
        }
    }
    pop_front(cut);
    cut_left = l;
    VERIFY(len() == expected);
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::addStep() {
    cut_right -= 1;
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::addStep(Edge &edge) {
    *this += Segment<Edge>(edge, 0, 1);
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::extend(const Sequence &seq) {
    VERIFY(valid());
    for (size_t cpos = 0; cpos < seq.size(); cpos++) {
        unsigned char c = seq[cpos];
        if (endClosed()) {
            Vertex &v = finish();
            if (v.hasOutgoing(c)) {
                Edge &edge = v.getOutgoing(c);
                addStep(edge);
            } else {
                invalidate();
                return *this;
            }
        } else {
            if (back().contig().truncSeq()[back().right] == c) {
                addStep();
            } else {
                invalidate();
                return *this;
            }
        }
    }
    return *this;
}

template<class Graph>
bool ag::GraphPath<Graph>::endClosed() const {
    return valid() && cutRight() == 0;
}

template<class Graph>
bool ag::GraphPath<Graph>::startClosed() const {
    return valid() && cutLeft() == 0;
}

//template<class Graph>
//unsigned char ag::GraphPath<Graph>::lastNucl() const {
//    Segment<Edge> seg = back();
//    return seg.truncSeq()[seg.right - 1];
//}

template<class Graph>
size_t ag::GraphPath<Graph>::cutLeft() const {
    return cut_left;
}

template<class Graph>
size_t ag::GraphPath<Graph>::cutRight() const {
    return cut_right;
}

template<class Graph>
std::vector<ag::GraphPath<Graph>> ag::GraphPath<Graph>::allSteps() {
    if (size() != 0 && cut_right > 0) {
        ag::GraphPath<Graph> copy = *this;
        return {std::move(copy.addStep())};
    }
    std::vector<ag::GraphPath<Graph>> res;
    Vertex &end = finish();
    for (Edge &edge : end) {
        GraphPath<Graph> copy = *this;
        res.emplace_back(std::move(copy.addStep(edge)));
    }
    return res;
}

template<class Graph>
std::vector<ag::GraphPath<Graph>> ag::GraphPath<Graph>::allExtensions(size_t len) {
    std::vector<GraphPath<Graph>> res = {*this};
    size_t left = 0;
    size_t right = 1;
    for (size_t l = 0; l < len; l++) {
        for (size_t i = left; i < right; i++) {
            std::vector<GraphPath<Graph>> tmp = res[i].allSteps();
            res.insert(res.end(), tmp.begin(), tmp.end());
        }
        left = right;
        right = res.size();
    }
    return std::move(res);
}

template<class Graph>
Sequence ag::GraphPath<Graph>::map(std::unordered_map<const Edge *, Sequence> &edge_map) {
    SequenceBuilder sb;
    bool start = true;
    for (Segment<Edge> seg : *this) {
        auto it = edge_map.find(&seg.contig());
        if (it == edge_map.end()) {
            if (start) {
                sb.append((start_->getSeq() + seg.contig().truncSeq()).Subseq(seg.left, seg.right + start_->getSeq().size()));
                start = false;
            } else {
                sb.append(seg.truncSeq());
            }
        } else {
            size_t left = start_->getSeq().size();
            if (start) {
                left = 0;
            }
            size_t right = start_->getSeq().size();
            size_t sz = it->second.size() - start_->getSeq().size();
            if (seg.left == 0 && seg.right == seg.contig().size()) {
                right += sz;
            } else if (seg.left == 0) {
                right += std::min(sz, seg.right);
            } else if (seg.right == seg.contig().size()) {
                left += sz - std::min(sz, seg.size());
                right += sz;
            } else {
                size_t l = seg.left * sz / seg.contig().size();
                left += l;
                right += std::min(l + seg.size(), sz);
            }
            sb.append(it->second.Subseq(left, right));
            start = false;
        }
    }
    return sb.BuildSequence();
}

template<class Graph>
Sequence ag::GraphPath<Graph>::Seq() const {
    if (!valid())
        return {};
    if(size() == 0) {
        return start_->getSeq().Subseq(cutLeft(), start_->size() - cutRight());
    }
    Edge & edge = frontEdge();
    Sequence seq = edge.getSeq();
    if(size() == 1) {
        return seq.Subseq(cutLeft(), seq .size() - cutRight());
    }

    SequenceBuilder sb;
    sb.append(seq.Subseq(cutLeft()));
    for(size_t i = skip_left + 1; i + 1 < path.size() - skip_right; i++) {
        sb.append(operator[](i).truncSeq());
    }
    sb.append(backEdge().truncSeq().Subseq(0, backEdge().truncSize() - cutRight()));
    return sb.BuildSequence();
}

template<class Graph>
Sequence ag::GraphPath<Graph>::truncSeq() const {
    SequenceBuilder sb;
    for (Segment<Edge> seg : *this) {
        sb.append(seg.truncSeq());
    }
    return sb.BuildSequence();
}

template<class Graph>
Sequence ag::GraphPath<Graph>::truncSubseq(size_t start_position, size_t sz) const {
    SequenceBuilder sb;
    for (size_t i = start_position; i < this->size(); i++) {
        Segment<Edge> seg = (*this)[i];
        if (seg.size() >= sz) {
            sb.append(seg.shrinkRightToLen(sz).truncSeq());
            sz = 0;
            break;
        } else {
            sb.append(seg.truncSeq());
            sz -= seg.size();
        }
    }
    return sb.BuildSequence();
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::reroute(size_t left, size_t right, const GraphPath<Graph> &rerouting) const {
    VERIFY(left == 0 || getVertex(left) == rerouting.start());
    VERIFY(right == size() || getVertex(right) == rerouting.finish());
    GraphPath<Graph> res;
    res += subPath(0, left);
    res += rerouting;
    res += subPath(right, size());
    return std::move(res);
}

template<class Graph>
void ag::GraphPath<Graph>::operator+=(const ag::GraphPath<Graph> &other) {
    if(other.size() == 0)
        return;
    if(!valid()) {
        *this = other;
        return;
    }
    for (Segment<Edge> al : other) {
        operator+=(al);
    }
}

template<class Graph>
void ag::GraphPath<Graph>::operator+=(const Segment<Edge> &other) {
    if(!valid()) {
        *this = {other};
        return;
    }
    if(cut_right == 0) {
        VERIFY(other.left == 0 && finish() == other.contig().getStart());
        if(skip_right == 0)
            path.emplace_back(other.contig().getId());
        else {
            path[path.size() - skip_right] = other.contig().getId();
            skip_right -= 1;
        }
        cut_right = other.contig().truncSize() - other.right;
    } else {
        VERIFY(cut_right == other.contig().truncSize() - other.left && other.contig() == backEdge());
    }
    cut_right = other.contig().truncSize() - other.right;
}

template<class Graph>
void ag::GraphPath<Graph>::operator+=(Edge &other) {
    GraphPath<Graph>::operator+=(Segment<Edge>(other, 0, other.truncSize()));
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::operator+(const GraphPath<Graph> &other) const {
    GraphPath<Graph> res = *this;
    res += other;
    return std::move(res);
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::operator+(const Segment<Edge> &other) const {
    GraphPath<Graph> res = *this;
    res += other;
    return std::move(res);
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::operator+(Edge &other) const {
    GraphPath<Graph> res = *this;
    res += other;
    return std::move(res);
}

template<class Graph>
std::string ag::GraphPath<Graph>::covStr(bool show_coverage) const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << cutLeft() << " " << start().getInnerId();
    for(const Segment<Edge> &seg : *this) {
        ss << " " << seg.size() << "/" <<seg.contig().truncSize() << seg.contig().nuclLabel() ;
        if(show_coverage) {
            ss << "(" << seg.contig().getCoverage() << ")";
        }
        ss << " " << seg.contig().getFinish().getInnerId();
    }
    ss << " " << cutRight();
    return ss.str();
}

template<class Graph>
std::string ag::GraphPath<Graph>::str() const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << cutLeft() << " " << start().getInnerId();
    for(const Segment<Edge> &seg : *this) {
        ss << " " << seg.size() << "/" <<seg.contig().truncSize() << seg.contig().nuclLabel() ;
        ss << " " << seg.contig().getFinish().getInnerId();
    }
    ss << " " << cutRight();
    return ss.str();
}

template<class Graph>
std::string ag::GraphPath<Graph>::lenStr() const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << skip_left << "_" << cutLeft() << "[" << start().getInnerId() << "(" << start().size() << ")";
    for(const typename Graph::Edge &edge : edges()) {
        ss << " -> " << edge.getInnerId().eid << edge.nuclLabel() << "(" << edge.rc().truncSize() << "|" << edge.truncSize() << ")->"  << edge.getFinish().getInnerId() << "(" << edge.getFinish().size() << ")" ;
    }
    ss << "]" << cutRight() << "_" << skip_right;
    return ss.str();
}

template<class Graph>
Segment<typename Graph::Edge> ag::GraphPath<Graph>::back() const {
    return {backEdge(), (size() == 1 ? cutLeft() : 0), backEdge().truncSize() - cutRight()};
}

template<class Graph>
Segment<typename Graph::Edge> ag::GraphPath<Graph>::front() const {
    return {frontEdge(), cutLeft(), size() == 1 ? frontEdge().truncSize() - cutRight() : frontEdge().truncSize()};
}

template<class Graph>
Segment<typename Graph::Edge> ag::GraphPath<Graph>::operator[](size_t i) const {
    return {getEdge(i), i == 0 ? cutLeft() : 0, i == size() - 1 ? backEdge().truncSize() - cutRight() : getEdge(i).truncSize()};
}

template<class Graph>
typename ag::GraphPath<Graph>::segment_iterator ag::GraphPath<Graph>::begin() const {
    std::function<Segment<Edge> (size_t)> transformer = [this](size_t ind)->Segment<Edge> {
        return operator[](ind);
    };
    return {CountingIterator<size_t>(0), {size()}, transformer};
}

template<class Graph>
typename ag::GraphPath<Graph>::segment_iterator ag::GraphPath<Graph>::end() const {
    std::function<Segment<Edge> (size_t)> transformer = [this](size_t ind)->Segment<Edge> {
        return operator[](ind);
    };
    return {{size()}, {size()}, transformer};
}

template<class Graph>
size_t ag::GraphPath<Graph>::len() const {
    if(!valid())
        return 0;
    size_t res = start().size();
    for(Edge &edge : edges()) {
        res += edge.truncSize();
    }
    return res - cutLeft() - cutRight();
}

template<class Traits>
ag::GraphPath<Traits> &ag::GraphPath<Traits>::fastExtend(const Sequence &seq) {
    size_t pos = cutRight();
    cut_right = 0;
    while(pos < seq.size()) {
        Edge &e = finish().getOutgoing(seq[pos]);
        *this += e;
        pos += e.truncSize();
    }
    cut_right += pos - seq.size();
    return *this;
}

template<class Traits>
void ag::GraphPath<Traits>::simpleReroute(size_t left, size_t right, const std::vector<EdgeId> &edges) {
    VERIFY(right - left == edges.size());
    for(size_t i = 0; i < edges.size(); i++) {
        path[skip_left + left + i] = edges[i];
    }
}
