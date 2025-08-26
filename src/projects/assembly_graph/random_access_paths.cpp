#include "random_access_paths.hpp"
using namespace ag;

Vertex &RAGraphPath::getVertex(size_t i) const {
    VERIFY(valid());
    VERIFY(i <= size());
    if (i == 0)
        return *start_;
    else
        return path[i - 1]->getFinish();
}


size_t RAGraphPath::find(Edge &edge, size_t pos) const {
    while (pos < size() && edge != getEdge(pos))
        pos++;
    if (pos == size())
        return -1;
    return pos;
}


size_t RAGraphPath::find(Vertex &v, size_t pos) const {
    while (pos <= size() && v != getVertex(pos))
        pos++;
    if (pos > size())
        return -1;
    return pos;
}


size_t RAGraphPath::truncLen() const {
    size_t res = 0;
    for (Edge &edge: edges())
        res += edge.truncSize();
    if (res == 0)
        return 0;
    return res - cut_left - cut_right;
}


IterableStorage<typename RAGraphPath::vertex_iterator>
RAGraphPath::vertices(size_t from, size_t to) const &{
    if (to == size_t(-1))
        to = size() + 1;
    std::function<Vertex &(size_t)> transformer = [this](size_t ind) -> Vertex & { return getVertex(ind); };
    CountingIterator<size_t> end_it = CountingIterator<size_t>(valid() ? to : 0);
    vertex_iterator vbegin(CountingIterator<size_t>(from), end_it, transformer);
    vertex_iterator vend(end_it, end_it, transformer);
    return {vbegin, vend};
}


IterableStorage<typename RAGraphPath::vertex_iterator> RAGraphPath::innerVertices() const &{
    return vertices(1, size());
}


IterableStorage<typename RAGraphPath::edge_iterator> RAGraphPath::edges() const &{
    std::function<Edge &(size_t)> transformer = [this](size_t ind) -> Edge & { return *path[ind]; };
    CountingIterator<size_t> end_it = CountingIterator<size_t>(path.size());
    edge_iterator ebegin(CountingIterator<size_t>(0), end_it, transformer);
    edge_iterator eend(end_it, end_it, transformer);
    return {ebegin, eend};
}


RAGraphPath RAGraphPath::RC() const {
    if (!valid())
        return {};
    std::vector<EdgeId> res;
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        res.emplace_back((*it)->rc().getId());
    }
    return {getFinish().rc(), res, rightCut(), leftCut()};
}


void RAGraphPath::invalidate() {
    start_ = {};
    path.clear();
    cut_left = 0;
    cut_right = 0;
}


bool RAGraphPath::valid() const {
    VERIFY(start_.valid() || size() == 0);
    return start_.valid();
}


RAGraphPath &RAGraphPath::cutBack(size_t l) {
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

RAGraphPath &RAGraphPath::cutFront(size_t l) {
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


RAGraphPath &RAGraphPath::addStep() {
    cut_right -= 1;
    return *this;
}


RAGraphPath &RAGraphPath::addStep(Edge &edge) {
    *this += Segment<Edge>(edge, 0, 1);
    return *this;
}


RAGraphPath &RAGraphPath::extend(const Sequence &seq) {
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


bool RAGraphPath::endClosed() const {
    return valid() && rightCut() == 0;
}


bool RAGraphPath::startClosed() const {
    return valid() && leftCut() == 0;
}

size_t RAGraphPath::leftCut() const {
    return cut_left;
}


size_t RAGraphPath::rightCut() const {
    return cut_right;
}


std::vector<RAGraphPath> RAGraphPath::allSteps() {
    if (size() != 0 && cut_right > 0) {
        RAGraphPath copy = *this;
        return {std::move(copy.addStep())};
    }
    std::vector<RAGraphPath> res;
    Vertex &end = getFinish();
    for (Edge &edge: end) {
        RAGraphPath copy = *this;
        res.emplace_back(std::move(copy.addStep(edge)));
    }
    return res;
}


std::vector<RAGraphPath> RAGraphPath::allExtensions(size_t len) {
    std::vector<RAGraphPath> res = {*this};
    size_t left = 0;
    size_t right = 1;
    for (size_t l = 0; l < len; l++) {
        for (size_t i = left; i < right; i++) {
            std::vector<RAGraphPath> tmp = res[i].allSteps();
            res.insert(res.end(), tmp.begin(), tmp.end());
        }
        left = right;
        right = res.size();
    }
    return std::move(res);
}


Sequence RAGraphPath::map(std::unordered_map<const Edge *, Sequence> &edge_map) {
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
            if (seg.left == 0 && seg.right == seg.contig().truncSize()) {
                right += sz;
            } else if (seg.left == 0) {
                right += std::min(sz, seg.right);
            } else if (seg.right == seg.contig().truncSize()) {
                left += sz - std::min(sz, seg.size());
                right += sz;
            } else {
                size_t l = seg.left * sz / seg.contig().truncSize();
                left += l;
                right += std::min(l + seg.size(), sz);
            }
            sb.append(it->second.Subseq(left, right));
            start = false;
        }
    }
    return sb.BuildSequence();
}


Sequence RAGraphPath::Seq() const {
    if(!valid())
        return {};
    if(empty()) {
        return getStart().getSeq().Subseq(cut_left, getStart().size() - cut_right);
    }
    SequenceBuilder sb;
    size_t first_right = size() == 1 ? cut_right : 0;
    if(cut_left >= getStart().size())
        sb.append(frontEdge().truncSeq().Subseq(cut_left - getStart().size(), frontEdge().truncSize() - first_right));
    else
        sb.append(getStart().getSeq().Subseq(cut_left) + frontEdge().truncSeq().Subseq(0, frontEdge().truncSize() - first_right));
    for (size_t i = 1; i < size(); i++) {
        Edge &edge = *path[i];
        size_t right = i + 1 == size() ? cut_right : 0;
        sb.append(edge.truncSeq().Subseq(0, edge.truncSize() - right));
    }
    return sb.BuildSequence();

}


Sequence RAGraphPath::truncSeq() const {
    SequenceBuilder sb;
    for (Segment<Edge> seg: *this) {
        sb.append(seg.truncSeq());
    }
    return sb.BuildSequence();
}


RAGraphPath
RAGraphPath::reroute(size_t left, size_t right, const RAGraphPath &rerouting) const {
    VERIFY(left == 0 || getVertex(left) == rerouting.getStart());
    VERIFY(right == size() || getVertex(right) == rerouting.getFinish());
    RAGraphPath res;
    res += subPath(0, left);
    res += rerouting;
    res += subPath(right, size());
    return std::move(res);
}


void RAGraphPath::operator+=(const RAGraphPath &other) {
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


void RAGraphPath::operator+=(const Segment<Edge> &other) {
    if (!valid()) {
        *this = {other};
        return;
    }
    if (!empty() && backEdge() == other.contig() && cut_right + other.left == backEdge().truncSize()) {
        cut_right = other.cutRight();
    } else {
        VERIFY(cut_right == 0 && other.cutLeft() == 0);
        VERIFY(getFinish() == other.contig().getStart());
        path.emplace_back(other.contig().getId());
        cut_right = other.cutRight();
    }
}


void RAGraphPath::operator+=(Edge &other) {
    RAGraphPath::operator+=(Segment<Edge>(other, 0, other.truncSize()));
}


RAGraphPath RAGraphPath::operator+(const RAGraphPath &other) const {
    RAGraphPath res = *this;
    res += other;
    return std::move(res);
}


RAGraphPath RAGraphPath::operator+(const Segment<Edge> &other) const {
    RAGraphPath res = *this;
    res += other;
    return std::move(res);
}


RAGraphPath RAGraphPath::operator+(Edge &other) const {
    RAGraphPath res = *this;
    res += other;
    return std::move(res);
}


std::string RAGraphPath::str() const {
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


Segment<Edge> RAGraphPath::back() const {
    return {backEdge(), (size() == 1 ? leftCut() : 0), backEdge().truncSize() - rightCut()};
}


Segment<Edge> RAGraphPath::front() const {
    return {frontEdge(), leftCut(), size() == 1 ? frontEdge().truncSize() - rightCut() : frontEdge().truncSize()};
}


Segment<Edge> RAGraphPath::operator[](size_t i) const {
    return {getEdge(i), i == 0 ? leftCut() : 0,
            i == size() - 1 ? backEdge().truncSize() - rightCut() : getEdge(i).truncSize()};
}


typename RAGraphPath::segment_iterator RAGraphPath::begin() const {
    std::function<Segment<Edge>(size_t)> transformer = [this](size_t ind) -> Segment<Edge> {
        return operator[](ind);
    };
    return {CountingIterator<size_t>(0), {size()}, transformer};
}


typename RAGraphPath::segment_iterator RAGraphPath::end() const {
    std::function<Segment<Edge>(size_t)> transformer = [this](size_t ind) -> Segment<Edge> {
        return operator[](ind);
    };
    return {{size()}, {size()}, transformer};
}


size_t RAGraphPath::len() const {
    if (!valid())
        return 0;
    size_t res = getStart().size();
    for (Edge &edge: edges()) {
        res += edge.truncSize();
    }
    return res - leftCut() - rightCut();
}


RAGraphPath &RAGraphPath::fastExtend(const Sequence &seq) {
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


void RAGraphPath::simpleReroute(size_t left, size_t right, const std::vector<EdgeId> &edges) {
    VERIFY(right - left == edges.size());
    for (size_t i = 0; i < edges.size(); i++) {
        path[left + i] = edges[i];
    }
}


const std::deque<EdgeId> RAGraphPath::edgeIds() const &{
    return path;
}


std::vector<EdgeId> RAGraphPath::asEdgeIds() const {return {edgeIds().begin(), edgeIds().end()};}


template<class Iterator>
RAGraphPath::RAGraphPath(Iterator begin, Iterator end) : cut_left(0), cut_right(0) {
    while (begin != end) {
        *this += *begin;
        ++begin;
    }
    if (size() > 0) {
        start_ = frontEdge().getStart().getId();
    }
}


Edge &RAGraphPath::getEdge(size_t i) const {
    VERIFY(i < size());
    return *path[i];
}


void RAGraphPath::pop_back() {
    size_t last_edge_size = path.back()->truncSize();
    path.pop_back();
    cut_right -= std::min(cut_right, last_edge_size);
}


void RAGraphPath::pop_back(size_t len) {
    for (size_t i = 0; i < len; i++)
        pop_back();
}


void RAGraphPath::pop_front() {
    size_t first_edge_size = path.front()->rc().truncSize();
    start_ = path.front()->getFinish().getId();
    path.pop_front();
    cut_left -= std::min(first_edge_size, cut_left);
}


void RAGraphPath::pop_front(size_t len) {
    for (size_t i = 0; i < len; i++)
        pop_front();
}


bool RAGraphPath::operator==(const RAGraphPath &other) const {
    return start_ == other.start_ && cut_left == other.cut_left && cut_right == other.cut_right &&
           size() == other.size() && std::equal(path.begin(), path.end(), other.path.begin());
}


bool RAGraphPath::operator!=(const RAGraphPath &other) const { return !operator==(other); }

RAPathIterator &RAPathIterator::operator++() {
    VERIFY(pos >= -1);
    rc ? pos-- : pos++;
    return *this;
}


RAPathIterator &RAPathIterator::operator--() {
    VERIFY(pos >= -1);
    rc ? pos++ : pos--;
    return *this;
}


bool RAPathIterator::operator==(const RAPathIterator &other) const {
    return rc == other.rc && path == other.path && pos == other.pos;
}


std::string RAPathIterator::str() const {
    std::stringstream ss;
    ss << path->str() << ":" << pos << "(" << (rc ? "B" : "F") << ")";
    return ss.str();
}

ag::RAPathDirection ag::RAGraphPath::forward() { return {*this, false}; }

ag::RAPathDirection ag::RAGraphPath::backward() { return {*this, true}; }



void RAPathDirection::rerouteSameSize(RAPathIterator from, RAPathIterator to,
                                      const ag::RAGraphPath &alt) const {
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


void RAPathDirection::cutFrontPrefix() const {
    VERIFY(!path->empty());
    if (!rc) {
        path->start_ = path->frontEdge().getFinish().getId();
        path->path.pop_front();
    } else {
        path->path.pop_back();
    }
}


RAPathIterator RAPathDirection::end() const {
    return {*path, rc, int(rc ? - 1 : path->path.size())};
}


RAPathIterator RAPathDirection::begin() const {
    return {*path, rc, int(rc ? path->path.size() - 1 : 0)};
}


void RAPathDirection::cutBackSuffix() const {
    RC().cutFrontPrefix();
}


void RAPathDirection::forceReroute(RAPathIterator from, RAPathIterator to,
                                   const RAGraphPath &alt) const {
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

ag::RAGraphPath ag::RAGraphPath::WalkForward(Edge &start) {
    RAGraphPath res(start);
    VertexId next = start.getFinish().getId();
    VERIFY(next.valid());
    while (*next != start.getStart() && *next != start.getStart().rc() && !next->isJunction()) {
        VERIFY(next.valid());
        VERIFY(next->outDeg() == 1);
        res += next->front();
        next = res.getFinish().getId();
    }
    return std::move(res);
}


ag::RAGraphPath ag::RAGraphPath::subPath(size_t from, size_t to) const {
    if (!valid()) {
        VERIFY(from == 0 && to == 0);
        return {};
    }
    if (from == to) {
        if ((from == 0 && leftCut() > 0) || (to == size() && rightCut() > 0)) {
            return {};
        } else {
            return RAGraphPath(getVertex(from));
        }
    } else {
        return {getVertex(from),
                std::vector<EdgeId>(path.begin() + from, path.begin() + to),
                from == 0 ? leftCut() : 0, to == size() ? rightCut() : 0};
    }
}
