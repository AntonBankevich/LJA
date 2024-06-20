#pragma once

#include "assembly_graph/graph_paths.hpp"
#include "assembly_graph/assembly_graph_base.hpp"
#include "sequences/contigs.hpp"
namespace old {
    namespace ag {
        template<class Traits>
        class RAPathDirection;

        template<class Traits>
        class RAPathIterator;

        template<class Traits>
        class RAGraphPath {
        public:
            typedef typename Traits::Vertex Vertex;
            typedef typename Traits::Edge Edge;
            typedef typename Vertex::VertexId VertexId;
            typedef typename Edge::EdgeId EdgeId;

            friend class RAPathIterator<Traits>;

            friend class RAPathDirection<Traits>;

        private:
            VertexId start_;
            std::vector<EdgeId> path;
            size_t skip_left = 0;
            size_t skip_right = 0;
            size_t cut_left;
            size_t cut_right;

            void set(Edge &edge, size_t pos) { path[pos + skip_left] = edge.getId(); }

        public:
            typedef typename std::vector<EdgeId>::iterator iterator;
            typedef typename std::vector<EdgeId>::const_iterator const_iterator;
            typedef TransformingIterator<CountingIterator<size_t>, Vertex> vertex_iterator;
            typedef TransformingIterator<CountingIterator<size_t>, Edge> edge_iterator;
            typedef Generator<CountingIterator<size_t>, Segment<Edge>> segment_iterator;

            RAGraphPath(Vertex &_start, std::vector<EdgeId> _path, size_t cut_left, size_t cut_right,
                      size_t skip_left = 0, size_t skip_right = 0) :
                    start_(_start.getId()), path(std::move(_path)), cut_left(cut_left), cut_right(cut_right),
                    skip_left(skip_left), skip_right(skip_right) {}

            RAGraphPath(Vertex &_start, size_t cut_left = 0, size_t cut_right = 0) : start_(_start.getId()), // NOLINT(google-explicit-constructor)
                                                                                   cut_left(cut_left), cut_right(
                            cut_right) {} // NOLINT(google-explicit-constructor)
            RAGraphPath(Edge &edge, size_t cut_left = 0, size_t cut_right = 0) : start_(edge.getStart().getId()), // NOLINT(google-explicit-constructor)
                                                                               path({edge.getId()}), cut_left(cut_left),
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
            Edge &backEdge() const { return *path[path.size() - skip_right - 1]; }
            Edge &frontEdge() const { return *path[skip_left]; }
            Segment<Edge> back() const;
            Segment<Edge> front() const;
            size_t size() const { return path.size() - skip_left - skip_right; }
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
            IterableStorage<typename std::vector<typename Traits::Edge::EdgeId>::const_iterator> edgeIds() const &;
            segment_iterator begin() const;
            segment_iterator end() const;

            RAPathDirection<Traits> forward() { return {*this, false}; }
            RAPathDirection<Traits> backward() { return {*this, true}; }
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
            void setCutLeft(size_t value) {cut_left = value;}
            void setCutRight(size_t value) {cut_right = value;}
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

        template<class Traits>
        class RAPathIterator {
        private:
            RAGraphPath<Traits> *path;
            bool rc;
            int pos;
        public:
            RAPathIterator(RAGraphPath<Traits> &path, bool rc, int pos) : path(&path), rc(rc), pos(pos) {
            }

        public:
            typedef typename Traits::Edge Edge;
            typedef typename Traits::Edge &reference;
            typedef typename Traits::Edge *pointer;

            reference operator*() const { return rc ? path->path[pos]->rc() : *path->path[pos]; }
            pointer operator->() const { return rc ? &path->path[pos]->rc() : &(*path->path[pos]); }
            void set(Edge &edge) const { path->set(rc ? edge.rc() : edge, pos - path->skip_left); }
            RAPathIterator &operator++();
            RAPathIterator &operator--();

            RAPathIterator operator+(int d) const { return {*path, rc, rc ? pos - d : pos + d}; }
            RAPathIterator operator-(int d) const { return operator+(-d); }
            RAPathIterator operator++(int) { return *this + 1; }
            RAPathIterator operator--(int) { return *this - 1; }

            bool operator==(const RAPathIterator &other) const;
            bool operator!=(const RAPathIterator &other) const { return !(*this == other); }

            int getPos() const { return pos; }

            std::string str() const;

//        Careful! RC does not point to the same element. Instead it makes sure that begin->end and end->begin
            RAPathIterator RC() const { return {*path, !rc, pos + (rc ? 1 : -1)}; }
            RAPathIterator<Traits> SameElementRC() const { return {*path, !rc, pos}; }
        };

        template<class Traits>
        RAPathIterator<Traits> &RAPathIterator<Traits>::operator++() {
            VERIFY(pos >= -1);
            rc ? pos-- : pos++;
            return *this;
        }

        template<class Traits>
        RAPathIterator<Traits> &RAPathIterator<Traits>::operator--() {
            VERIFY(pos >= -1);
            rc ? pos++ : pos--;
            return *this;
        }

        template<class Traits>
        bool RAPathIterator<Traits>::operator==(const RAPathIterator<Traits> &other) const {
            return rc == other.rc && path == other.path && pos == other.pos;
        }

        template<class Traits>
        std::string RAPathIterator<Traits>::str() const {
            std::stringstream ss;
            ss << path->str() << ":" << pos << "(" << (rc ? "B" : "F") << ")";
            return ss.str();
        }

        template<class Traits>
        class RAPathDirection {
        private:
            RAGraphPath<Traits> *path;
            bool rc;

            void setCutLeft(size_t val) const { rc ? path->cut_right = val : path->cut_left = val; }
            void setCutRight(size_t val) const { rc ? path->cut_left = val : path->cut_right = val; }
        public:
            typedef typename Traits::Vertex Vertex;
            typedef typename Traits::Edge Edge;
            typedef typename Vertex::VertexId VertexId;
            typedef typename Edge::EdgeId EdgeId;

            RAPathDirection(RAGraphPath<Traits> &path, bool rc) : path(&path), rc(rc) {}
            Vertex &getStart() const { return rc ? path->getFinish().rc() : path->getStart(); }
            Vertex &getFinish() const { return rc ? path->getStart().rc() : path->getFinish(); }
            size_t cutRight() const { return rc ? path->leftCut() : path->rightCut(); }
            size_t cutLeft() const { return rc ? path->rightCut() : path->leftCut(); }
            void pop_front() const { return rc ? path->pop_back() : path->pop_front(); }
            void pop_back() const { return rc ? path->pop_front() : path->pop_back(); }
            bool isForward() const { return !rc; }
            size_t size() const { return path->size(); }
            RAGraphPath<Traits> &getPath() const { return *path; }
            RAPathIterator<Traits> end() const;
            RAPathIterator<Traits> begin() const;
            RAPathDirection RC() const { return {*path, !rc}; }
            bool valid() const { return path->valid(); }
            bool empty() const { return path->empty(); }
            bool isRC() const { return rc; }
            void invalidate() { path->invalidate(); }
            void rerouteSameSize(RAPathIterator<Traits> from, RAPathIterator<Traits> to, const ag::RAGraphPath<Traits> &alt) const;
            void cutFrontPrefix() const;
            void cutBackSuffix() const;
            void forceReroute(RAPathIterator<Traits> from, RAPathIterator<Traits> to, const RAGraphPath<Traits> &alt) const;
        };

        template<class Traits>
        void RAPathDirection<Traits>::rerouteSameSize(RAPathIterator<Traits> from, RAPathIterator<Traits> to,
                                                      const ag::RAGraphPath<Traits> &alt) const {
            if (from == begin())
                setCutLeft(alt.leftCut());
            if (to == end())
                setCutRight(alt.rightCut());
            size_t cur = 0;
            while (from != to && cur < alt.size()) {
                from.set(alt.getEdge(cur));
                cur++;
                ++from;
            }
            VERIFY(from == to);
            VERIFY(cur == alt.size());
            if (!empty())
                path->start_ = path->frontEdge().getStart().getId();
        }

        template<class Traits>
        void RAPathDirection<Traits>::cutFrontPrefix() const {
            VERIFY(!path->empty());
            if (!rc) {
                path->start = path->frontEdge().getFinish();
                path->path.pop_front();
            } else {
                path->path.pop_back();
            }
        }

        template<class Traits>
        RAPathIterator<Traits> RAPathDirection<Traits>::end() const {
            return {*path, rc, int(rc ? path->skip_left - 1 : path->path.size() - path->skip_right)};
        }

        template<class Traits>
        RAPathIterator<Traits> RAPathDirection<Traits>::begin() const {
            return {*path, rc, int(rc ? path->path.size() - path->skip_right - 1 : path->skip_left)};
        }

        template<class Traits>
        void RAPathDirection<Traits>::cutBackSuffix() const {
            RC().cutFrontPrefix();
        }

        template<class Traits>
        void RAPathDirection<Traits>::forceReroute(RAPathIterator<Traits> from, RAPathIterator<Traits> to,
                                                   const RAGraphPath<Traits> &alt) const {
            if (rc) {
                RAPathDirection rc_dir = RC();
                rc_dir.forceReroute(to.RC(), from.RC(), alt.RC());
            }
            if (from == begin() && to == end()) {
                *path = alt;
            } else if (from == begin()) {
                *path = alt + path->subPath(from.getPos());
            } else if (to == end()) {
                path->pop_back(path->size() - from.getPos());
                *path += alt;
            } else {
                *path = path->subPath(0, from.getPos()) + alt + path->subPath(to.getPos());
            }
        }
    }

    template<class Graph>
    ag::RAGraphPath<Graph> ag::RAGraphPath<Graph>::WalkForward(Edge &start) {
        RAGraphPath<Graph> res(start);
        Vertex *next = &start.getFinish();
        VERIFY(next != nullptr);
        while (*next != start.getStart() && *next != start.getStart().rc() && !next->isJunction()) {
            VERIFY(next != nullptr);
            res += next->front();
            next = &res.getFinish();
        }
        return std::move(res);
    }

    template<class Graph>
    ag::RAGraphPath<Graph> ag::RAGraphPath<Graph>::subPath(size_t from, size_t to) const {
        if (!valid()) {
            VERIFY(from == 0 && to == 0);
            return {};
        }
        if (from == to) {
            if ((from == 0 && leftCut() > 0) || (to == size() && rightCut() > 0)) {
                return {};
            } else {
                return RAGraphPath<Graph>(getVertex(from));
            }
        } else {
            return {getVertex(from),
                    std::vector<EdgeId>(path.begin() + skip_left + from, path.begin() + skip_left + to),
                    from == 0 ? leftCut() : 0, to == size() ? rightCut() : 0};
        }
    }

    template<class Graph>
    typename Graph::Vertex &ag::RAGraphPath<Graph>::getVertex(size_t i) const {
        VERIFY(valid());
        VERIFY(i <= size());
        if (i == 0)
            return *start_;
        else
            return path[skip_left + i - 1]->getFinish();
    }

    template<class Graph>
    size_t ag::RAGraphPath<Graph>::find(Edge &edge, size_t pos) const {
        while (pos < size() && edge != getEdge(pos))
            pos++;
        if (pos == size())
            return -1;
        return pos;
    }

    template<class Graph>
    size_t ag::RAGraphPath<Graph>::find(Vertex &v, size_t pos) const {
        while (pos <= size() && v != getVertex(pos))
            pos++;
        if (pos > size())
            return -1;
        return pos;
    }

    template<class Graph>
    size_t ag::RAGraphPath<Graph>::truncLen() const {
        size_t res = 0;
        for (Edge &edge: edges())
            res += edge.truncSize();
        if (res == 0)
            return 0;
        return res - cut_left - cut_right;
    }

    template<class Graph>
    IterableStorage<typename ag::RAGraphPath<Graph>::vertex_iterator>
    ag::RAGraphPath<Graph>::vertices(size_t from, size_t to) const &{
        if (to == size_t(-1))
            to = size() + 1;
        from += skip_left;
        to += skip_left;
        std::function<Vertex &(size_t)> transformer = [this](size_t ind) -> Vertex & { return getVertex(ind); };
        CountingIterator<size_t> end_it = CountingIterator<size_t>(valid() ? to : 0);
        vertex_iterator vbegin(CountingIterator<size_t>(from), end_it, transformer);
        vertex_iterator vend(end_it, end_it, transformer);
        return {vbegin, vend};
    }

    template<class Graph>
    IterableStorage<typename ag::RAGraphPath<Graph>::vertex_iterator> ag::RAGraphPath<Graph>::innerVertices() const &{
        return vertices(1, size());
    }

    template<class Graph>
    IterableStorage<typename ag::RAGraphPath<Graph>::edge_iterator> ag::RAGraphPath<Graph>::edges() const &{
        std::function<Edge &(size_t)> transformer = [this](size_t ind) -> Edge & { return *path[ind]; };
        CountingIterator<size_t> end_it = CountingIterator<size_t>(path.size() - skip_right);
        edge_iterator ebegin(CountingIterator<size_t>(skip_left), end_it, transformer);
        edge_iterator eend(end_it, end_it, transformer);
        return {ebegin, eend};
    }

    template<class Graph>
    ag::RAGraphPath<Graph> ag::RAGraphPath<Graph>::RC() const {
        if (!valid())
            return {};
        std::vector<EdgeId> res;
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            res.emplace_back((*it)->rc().getId());
        }
        return {getFinish().rc(), res, rightCut(), leftCut(), skip_right, skip_left};
    }

    template<class Graph>
    void ag::RAGraphPath<Graph>::invalidate() {
        start_ = {};
        path.clear();
        cut_left = 0;
        cut_right = 0;
        skip_left = 0;
        skip_right = 0;
    }

    template<class Graph>
    bool ag::RAGraphPath<Graph>::valid() const {
        VERIFY(start_.valid() || size() == 0);
        return start_.valid();
    }

    template<class Graph>
    ag::RAGraphPath<Graph> &ag::RAGraphPath<Graph>::cutBack(size_t l) {
        VERIFY(l <= len());
        size_t expected = len() - l;
        size_t cur_cut = 0;
        size_t cut = 0;
        l += cut_right;
        cut_right = 0;
        while (cur_cut < size() && l >= getEdge(size() - 1 - cur_cut).truncSize()) {
            if (getEdge(size() - 1 - cur_cut).truncSize() == 0) {
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
//ag::RAGraphPath<Graph> &ag::RAGraphPath<Graph>::uniqueExtendBack(size_t l) {
//    if(cut_right != 0) {
//        size_t tmp = std::min(l, cut_right);
//        l -= tmp;
//        cut_right -= tmp;
//    }
//    while(l > 0) {
//        VERIFY(getFinish().outDeg() == 1);
//        Edge &e = getFinish().front();
//        size_t tmp = std::min(e.truncSize(), l);
//        *this += e;
//        cutBack(e.truncSize() - tmp);
//    }
//    return *this;
//}


////TODO: Optimize
//template<class Graph>
//ag::RAGraphPath<Graph> &ag::RAGraphPath<Graph>::uniqueExtendFront(size_t l) {
//    *this = this->RC().uniqueExtendBack(l);
//    return *this;
//}

    template<class Graph>
    ag::RAGraphPath<Graph> &ag::RAGraphPath<Graph>::cutFront(size_t l) {
        VERIFY(l <= len());
        size_t expected = len() - l;
        size_t cur_cut = 0;
        size_t cut = 0;
        l += cut_left;
        cut_left = 0;
        while (cur_cut < size() && l >= getEdge(cur_cut).rc().truncSize()) {
            if (getEdge(cur_cut).rc().truncSize() == 0) {
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
    ag::RAGraphPath<Graph> &ag::RAGraphPath<Graph>::addStep() {
        cut_right -= 1;
        return *this;
    }

    template<class Graph>
    ag::RAGraphPath<Graph> &ag::RAGraphPath<Graph>::addStep(Edge &edge) {
        *this += Segment<Edge>(edge, 0, 1);
        return *this;
    }

    template<class Graph>
    ag::RAGraphPath<Graph> &ag::RAGraphPath<Graph>::extend(const Sequence &seq) {
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
    bool ag::RAGraphPath<Graph>::endClosed() const {
        return valid() && rightCut() == 0;
    }

    template<class Graph>
    bool ag::RAGraphPath<Graph>::startClosed() const {
        return valid() && leftCut() == 0;
    }

//template<class Graph>
//unsigned char ag::RAGraphPath<Graph>::lastNucl() const {
//    Segment<Edge> seg = back();
//    return seg.truncSeq()[seg.right - 1];
//}

    template<class Graph>
    size_t ag::RAGraphPath<Graph>::leftCut() const {
        return cut_left;
    }

    template<class Graph>
    size_t ag::RAGraphPath<Graph>::rightCut() const {
        return cut_right;
    }

    template<class Graph>
    std::vector<ag::RAGraphPath<Graph>> ag::RAGraphPath<Graph>::allSteps() {
        if (size() != 0 && cut_right > 0) {
            ag::RAGraphPath<Graph> copy = *this;
            return {std::move(copy.addStep())};
        }
        std::vector<ag::RAGraphPath<Graph>> res;
        Vertex &end = getFinish();
        for (Edge &edge: end) {
            RAGraphPath<Graph> copy = *this;
            res.emplace_back(std::move(copy.addStep(edge)));
        }
        return res;
    }

    template<class Graph>
    std::vector<ag::RAGraphPath<Graph>> ag::RAGraphPath<Graph>::allExtensions(size_t len) {
        std::vector<RAGraphPath<Graph>> res = {*this};
        size_t left = 0;
        size_t right = 1;
        for (size_t l = 0; l < len; l++) {
            for (size_t i = left; i < right; i++) {
                std::vector<RAGraphPath<Graph>> tmp = res[i].allSteps();
                res.insert(res.end(), tmp.begin(), tmp.end());
            }
            left = right;
            right = res.size();
        }
        return std::move(res);
    }

    template<class Graph>
    Sequence ag::RAGraphPath<Graph>::map(std::unordered_map<const Edge *, Sequence> &edge_map) {
        SequenceBuilder sb;
        bool start = true;
        for (Segment<Edge> seg: *this) {
            auto it = edge_map.find(&seg.contig());
            if (it == edge_map.end()) {
                if (start) {
                    sb.append((start_->getSeq() + seg.contig().truncSeq()).Subseq(seg.left,
                                                                                  seg.right + start_->getSeq().size()));
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
    Sequence ag::RAGraphPath<Graph>::Seq() const {
        if (!valid())
            return {};
        if (size() == 0) {
            return start_->getSeq().Subseq(leftCut(), start_->size() - rightCut());
        }
        Edge &edge = frontEdge();
        Sequence seq = edge.getSeq();
        if (size() == 1) {
            return seq.Subseq(leftCut(), seq.size() - rightCut());
        }

        SequenceBuilder sb;
        sb.append(seq.Subseq(leftCut()));
        for (size_t i = skip_left + 1; i + 1 < path.size() - skip_right; i++) {
            sb.append(operator[](i).truncSeq());
        }
        sb.append(backEdge().truncSeq().Subseq(0, backEdge().truncSize() - rightCut()));
        return sb.BuildSequence();
    }

    template<class Graph>
    Sequence ag::RAGraphPath<Graph>::truncSeq() const {
        SequenceBuilder sb;
        for (Segment<Edge> seg: *this) {
            sb.append(seg.truncSeq());
        }
        return sb.BuildSequence();
    }

    template<class Graph>
    ag::RAGraphPath<Graph>
    ag::RAGraphPath<Graph>::reroute(size_t left, size_t right, const RAGraphPath<Graph> &rerouting) const {
        VERIFY(left == 0 || getVertex(left) == rerouting.getStart());
        VERIFY(right == size() || getVertex(right) == rerouting.getFinish());
        RAGraphPath<Graph> res;
        res += subPath(0, left);
        res += rerouting;
        res += subPath(right, size());
        return std::move(res);
    }

    template<class Graph>
    void ag::RAGraphPath<Graph>::operator+=(const ag::RAGraphPath<Graph> &other) {
        if (other.size() == 0)
            return;
        if (!valid()) {
            *this = other;
            return;
        }
        for (Segment<Edge> al: other) {
            operator+=(al);
        }
    }

    template<class Graph>
    void ag::RAGraphPath<Graph>::operator+=(const Segment<Edge> &other) {
        if (!valid()) {
            *this = {other};
            return;
        }
        if (!empty() && backEdge() == other.contig() && cut_right + other.left == backEdge().truncSize()) {
            cut_right = other.cutRight();
        } else {
            VERIFY(cut_right == 0 && other.cutLeft() == 0);
            VERIFY(getFinish() == other.contig().getStart());
            if (skip_right == 0)
                path.emplace_back(other.contig().getId());
            else {
                path[path.size() - skip_right] = other.contig().getId();
                skip_right -= 1;
            }
            cut_right = other.cutRight();
        }
    }

    template<class Graph>
    void ag::RAGraphPath<Graph>::operator+=(Edge &other) {
        RAGraphPath<Graph>::operator+=(Segment<Edge>(other, 0, other.truncSize()));
    }

    template<class Graph>
    ag::RAGraphPath<Graph> ag::RAGraphPath<Graph>::operator+(const RAGraphPath<Graph> &other) const {
        RAGraphPath<Graph> res = *this;
        res += other;
        return std::move(res);
    }

    template<class Graph>
    ag::RAGraphPath<Graph> ag::RAGraphPath<Graph>::operator+(const Segment<Edge> &other) const {
        RAGraphPath<Graph> res = *this;
        res += other;
        return std::move(res);
    }

    template<class Graph>
    ag::RAGraphPath<Graph> ag::RAGraphPath<Graph>::operator+(Edge &other) const {
        RAGraphPath<Graph> res = *this;
        res += other;
        return std::move(res);
    }

    template<class Graph>
    std::string ag::RAGraphPath<Graph>::str() const {
        if (!valid())
            return "";
        std::stringstream ss;
        ss << skip_left << "_" << leftCut() << "[" << getStart().getInnerId() << "(" << getStart().size() << ")";
        for (const typename Graph::Edge &edge: edges()) {
            ss << " -> " << edge.getInnerId().eid << edge.firstNucl() << "(" << edge.rc().truncSize() << "|"
               << edge.truncSize() << ")->" << edge.getFinish().getInnerId() << "(" << edge.getFinish().size() << ")";
        }
        ss << "]" << rightCut() << "_" << skip_right;
        return ss.str();
    }

    template<class Graph>
    Segment<typename Graph::Edge> ag::RAGraphPath<Graph>::back() const {
        return {backEdge(), (size() == 1 ? leftCut() : 0), backEdge().truncSize() - rightCut()};
    }

    template<class Graph>
    Segment<typename Graph::Edge> ag::RAGraphPath<Graph>::front() const {
        return {frontEdge(), leftCut(), size() == 1 ? frontEdge().truncSize() - rightCut() : frontEdge().truncSize()};
    }

    template<class Graph>
    Segment<typename Graph::Edge> ag::RAGraphPath<Graph>::operator[](size_t i) const {
        return {getEdge(i), i == 0 ? leftCut() : 0,
                i == size() - 1 ? backEdge().truncSize() - rightCut() : getEdge(i).truncSize()};
    }

    template<class Graph>
    typename ag::RAGraphPath<Graph>::segment_iterator ag::RAGraphPath<Graph>::begin() const {
        std::function<Segment<Edge>(size_t)> transformer = [this](size_t ind) -> Segment<Edge> {
            return operator[](ind);
        };
        return {CountingIterator<size_t>(0), {size()}, transformer};
    }

    template<class Graph>
    typename ag::RAGraphPath<Graph>::segment_iterator ag::RAGraphPath<Graph>::end() const {
        std::function<Segment<Edge>(size_t)> transformer = [this](size_t ind) -> Segment<Edge> {
            return operator[](ind);
        };
        return {{size()}, {size()}, transformer};
    }

    template<class Graph>
    size_t ag::RAGraphPath<Graph>::len() const {
        if (!valid())
            return 0;
        size_t res = getStart().size();
        for (Edge &edge: edges()) {
            res += edge.truncSize();
        }
        return res - leftCut() - rightCut();
    }

    template<class Traits>
    ag::RAGraphPath<Traits> &ag::RAGraphPath<Traits>::fastExtend(const Sequence &seq) {
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
    void ag::RAGraphPath<Traits>::simpleReroute(size_t left, size_t right, const std::vector<EdgeId> &edges) {
        VERIFY(right - left == edges.size());
        for (size_t i = 0; i < edges.size(); i++) {
            path[skip_left + left + i] = edges[i];
        }
    }

    template<class Traits>
    IterableStorage<typename std::vector<typename Traits::Edge::EdgeId>::const_iterator>
    ag::RAGraphPath<Traits>::edgeIds() const &{
        return {path.begin() + skip_left, path.end() - skip_right};
    }

    template<class Traits>
    std::vector<typename Traits::Edge::EdgeId> ag::RAGraphPath<Traits>::asEdgeIds() const {return {edgeIds().begin(), edgeIds().end()};}

    template<class Traits>
    template<class Iterator>
    ag::RAGraphPath<Traits>::RAGraphPath(Iterator begin, Iterator end) : cut_left(0), cut_right(0) {
        while (begin != end) {
            *this += *begin;
            ++begin;
        }
        if (size() > 0) {
            start_ = &frontEdge().getStart();
        }
    }

    template<class Traits>
    void ag::RAGraphPath<Traits>::normalize() {
        if (!valid())
            return;
        while (!empty() && frontEdge().isPrefix())
            pop_front();
        while (!empty() && backEdge().isSuffix())
            pop_back();
        while (getStart().inDeg() == 1 && getStart().rc().front().isSuffix() &&
                leftCut() < getStart().rc().front().getFinish().size()) {
            if (skip_left == 0) {
                path.insert(path.begin(), getStart().rc().front().rc().getId());
            } else {
                skip_left -= 1;
                path[skip_left] = getStart().rc().front().rc().getId();
            }
            start_ = frontEdge().getStart().getId();
        }
        while (getFinish().outDeg() == 1 && getFinish().front().isSuffix() &&
                rightCut() < getFinish().front().getFinish().size()) {
            if (skip_right == 0) {
                path.push_back(getFinish().front().getId());
            } else {
                path[path.size() - skip_right] = getFinish().front().getId();
                skip_right -= 1;
            }
        }
    }

    template<class Traits>
    typename Traits::Edge &ag::RAGraphPath<Traits>::getEdge(size_t i) const {
        if (i >= size()) {
            std::cout << i << " " << path.size() << " " << skip_left << " " << skip_right << std::endl;
        }
        VERIFY(i < size());
        return *path[skip_left + i];
    }

    template<class Traits>
    void ag::RAGraphPath<Traits>::pop_back() {
        skip_right++;
        VERIFY(skip_left + skip_right <= path.size());
        size_t new_cut_size = path[path.size() - skip_right]->truncSize();
        path[path.size() - skip_right] = {};
        cut_right -= std::min(cut_right, new_cut_size);
    }

    template<class Traits>
    void ag::RAGraphPath<Traits>::pop_back(size_t len) {
        for (size_t i = 0; i < len; i++)
            pop_back();
    }

    template<class Traits>
    void ag::RAGraphPath<Traits>::pop_front() {
        skip_left++;
        VERIFY(skip_left + skip_right <= path.size());
        size_t new_cut_size = path[skip_left - 1]->rc().truncSize();
        start_ = path[skip_left - 1]->getFinish().getId();
        path[skip_left - 1] = {};
        cut_left -= std::min(new_cut_size, cut_left);
    }

    template<class Traits>
    void ag::RAGraphPath<Traits>::pop_front(size_t len) {
        for (size_t i = 0; i < len; i++)
            pop_front();
    }

    template<class Traits>
    bool ag::RAGraphPath<Traits>::operator==(const ag::RAGraphPath<Traits> &other) const {
        return start_ == other.start_ && cut_left == other.cut_left && cut_right == other.cut_right &&
               size() == other.size() &&
               std::equal(path.begin() + skip_left, path.end() - skip_right,
                          other.path.begin() + other.skip_left);
    }

    template<class Traits>
    bool ag::RAGraphPath<Traits>::operator!=(const ag::RAGraphPath<Traits> &other) const { return !operator==(other); }
}