#include "graph_paths.hpp"

using namespace ag;

PathDirection GraphPath::forward() {return {*this, false};}

PathDirection GraphPath::backward() {return {*this, true};}

ConstPathDirection GraphPath::forward() const {return {*this, false};}

ConstPathDirection GraphPath::backward() const {return {*this, true};}

IterableStorage<PathVertexIterator> ConstPathDirection::vertices() {
    return {{*this, firstPosition()}, {*this, endPosition()}};
}

PathIterator ConstPathDirection::begin() const {return {firstPosition()};}

PathIterator ConstPathDirection::end() const {return {lastPosition()};}

Segment<Edge> ConstPathDirection::getSegment(PathPosition position) const {
    Edge &next = position.nextEdge();
    size_t cut_left = firstPosition() == position ? cutLeft() : 0;
    size_t cut_right = lastPosition() == position + 1 ? cutRight() : 0;
    return {next, cut_left, next.truncSize() - cut_right};
}

PathVertexIterator &PathVertexIterator::operator++() {
    if(position == dir.lastPosition())
        position = dir.endPosition();
    else
        ++position;
    return *this;
}

PathVertexIterator PathVertexIterator::operator++(int) { PathVertexIterator res = *this; ++res; return std::move(res); }

bool PathVertexIterator::operator==(const PathVertexIterator &other) const {return dir == other.dir && position == other.position;}

bool PathVertexIterator::operator!=(const PathVertexIterator &other) const {return !(*this == other);}

IterableStorage<PathVertexIterator> GraphPath::vertices() const &{
    if(valid())
        return {{forward(), firstPosition()}, {forward(), {endPosition()}}};
    else
        return {{forward(), endPosition()}, {forward(), {endPosition()}}};
}

IterableStorage<PathVertexIterator> GraphPath::innerVertices() const & {
    if(empty())
        return {{forward(), endPosition()},{forward(), endPosition()}};
    return {{forward(), firstPosition() + 1},{forward(), lastPosition()}};
}

IterableStorage<PathIterator> GraphPath::edges() const & {
    if (empty())
        return {{firstPosition()}, firstPosition()};
    return {{firstPosition()}, {lastPosition()}};
}

SegmentIterator GraphPath::begin() const {return {forward(), forward().firstPosition()};}

SegmentIterator GraphPath::end() const {return {forward(), forward().lastPosition()};}

Sequence GraphPath::truncSeq() const {
    SequenceBuilder sb;
    PathPosition first = firstPosition();
    PathPosition last = lastPosition();
    for (PathPosition pp = first; pp != last;) {
        Edge &edge = pp.nextEdge();
        size_t left = pp == first ? cut_left : 0;
        ++pp;
        size_t right = pp == last ? cut_right : 0;
        sb.append(edge.truncSeq().Subseq(left, edge.truncSize() - right));
    }
    return sb.BuildSequence();
}

Sequence GraphPath::Seq() const {
    if(!valid())
        return {};
    PathPosition left = firstPosition();
    PathPosition right = lastPosition();
    size_t left_cut = cut_left;
    size_t right_cut = cut_right;
    while(left != right && left.nextEdge().rc().truncSize() <= left_cut) {
        left_cut -= left.nextEdge().rc().truncSize();
        ++left;
    }
    while(left != right && right.prevEdge().truncSize() <= right_cut) {
        right_cut -= right.prevEdge().truncSize();
        --right;
    }
    if(left == right) {
        return left.getVertex().getSeq().Subseq(left_cut, left.getVertex().size() - right_cut);
    }
    if(left + 1 == right) {
        if(left_cut >= left.getVertex().size())
            return left.nextEdge().truncSeq().Subseq(left_cut - left.getVertex().size(), left.nextEdge().truncSize() - right_cut);
        else
            return left.getVertex().getSeq().Subseq(left_cut) + left.nextEdge().truncSeq().Subseq(0, left.nextEdge().truncSize() - right_cut);
    }
    SequenceBuilder sb;
    sb.append(left.nextEdge().rc().truncSeq().rc().Subseq(left_cut));
    ++left;
    sb.append(left.getVertex().getSeq());
    --right;
    for(;left != right; ++left) sb.append(left.nextEdge().truncSeq());
    sb.append(right.nextEdge().truncSeq().Subseq(0, right.nextEdge().truncSize() - right_cut));
    return sb.BuildSequence();
}

size_t GraphPath::truncLen() const {
    return valid() ? len() - start->size() : 0;
}

size_t GraphPath::len() const {
    size_t res = start->size();
    for (Edge &edge : edges())
        res += edge.truncSize();
    return res - cut_left - cut_right;
}

Segment<Edge> GraphPath::back() const {
    return {backEdge(), (isSingleton() ? leftCut() : 0), backEdge().truncSize() - rightCut()};
}

Segment<Edge> GraphPath::front() const {
    return {frontEdge(), leftCut(), isSingleton() ? frontEdge().truncSize() - rightCut() : frontEdge().truncSize()};
}

PathPosition GraphPath::firstPosition() const {
    return valid() ? PathPosition(start, fsplits.begin(), rsplits.end()) : endPosition();
}

PathPosition GraphPath::lastPosition() const {
    return valid() ? PathPosition(getFinish(), fsplits.end(), rsplits.begin()) : endPosition();
}

PathPosition GraphPath::endPosition() const {
    return {fsplits.end(), rsplits.begin()};
}

PathPosition GraphPath::rcEndPosition() const {
    return {rsplits.end(), fsplits.begin()};
}

GraphPath GraphPath::subPath(PathPosition from) const { return subPath(from, lastPosition()); }

GraphPath GraphPath::subPath(PathPosition from, PathPosition to) const {
//        Handles the case when the path is a vertex segment
    if(from == firstPosition() && to == lastPosition())
        return {*this};
    if(!valid()) {
        VERIFY(from == endPosition());
        VERIFY(to == endPosition());
        return {};
    }
    GraphPath res(from.getVertex());
    PathPosition cur = from;
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

void GraphPath::operator+=(const GraphPath &other) {
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

void GraphPath::operator+=(const Segment<Edge> &other) {
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

void GraphPath::operator+=(Edge &other) {
    VERIFY(cut_right == 0);
    VERIFY(!valid() || getFinish() == other.getStart());
    if(!valid())
        start = other.getStart().getId();
    fsplits.push_back(other.getCode());
    rsplits.push_front(other.rc().getCode());
    rc_start = other.getFinish().rc().getId();
}

GraphPath GraphPath::operator+(const GraphPath &other) const {
    GraphPath res = *this;
    res += other;
    return std::move(res);
}

GraphPath GraphPath::operator+(const Segment<Edge> &other) const {
    GraphPath res = *this;
    res += other;
    return std::move(res);
}

GraphPath GraphPath::operator+(Edge &other) const {
    GraphPath res = *this;
    res += other;
    return std::move(res);
}

void GraphPath::pop_back(size_t len) {
    for (size_t i = 0; i < len; i++)
        pop_back();
}

void GraphPath::pop_front(size_t len) {
    for (size_t i = 0; i < len; i++)
        pop_front();
}

GraphPath &GraphPath::cutBack(size_t l) {
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

GraphPath &GraphPath::cutFront(size_t l) {
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

PathIterator PathDirection::begin() const {
    return {firstPosition()};
}

PathIterator PathDirection::end() const {
    return {lastPosition()};
}

IterableStorage<PathVertexIterator> PathDirection::vertices() {
    return {{*this, firstPosition()}, {*this, endPosition()}};
}

void PathDirection::setCutLeft(size_t val) const {
    if(rc)
        path->setCutRight(val);
    else
        path->setCutLeft(val);
}

void PathDirection::setCutRight(size_t val) const {
    if(rc)
        path->setCutLeft(val);
    else
        path->setCutRight(val);
}

Segment<Edge> PathDirection::getSegment(PathPosition position) const {
    Edge &next = position.nextEdge();
    size_t cut_left = firstPosition() == position ? leftCut() : 0;
    size_t cut_right = lastPosition() == position + 1 ? rightCut() : 0;
    return {next, cut_left, next.truncSize() - cut_right};
}

PathDirection &PathDirection::operator=(GraphPath &&other) {
    if(rc)
        *path = other.RC();
    else
        *path = std::move(other);
    return *this;
}

PathDirection &PathDirection::operator=(const GraphPath &other) {
    if(rc)
        *path = other.RC();
    else
        *path = other;
    return *this;
}

void PathDirection::forcePushBack(Edge &edge) const {
    if(rc)
        path->forcePushFront(edge.rc());
    else
        path->forcePushBack(edge);
}

void PathDirection::forcePushFront(Edge &edge) const {
    if(rc)
        path->forcePushBack(edge.rc());
    else
        path->forcePushFront(edge);
}

bool PathDirection::isSingleton() const {
    return !empty() && getFSplits().size() == frontEdge().getCode().size() && getRSplits().size() == frontEdge().rc().getCode().size();
}

std::string GraphPath::str() const {
    if (!valid())
        return "";
    std::stringstream ss;
    if (isLegacy()) {
        ss << "Legacy:" << leftCut()  << "[" << getStart().getInnerId() << "(" << getStart().size() << ")]" <<rightCut();
        return ss.str();
    }
    ss << leftCut() << "[" << getStart().getInnerId() << "(" << getStart().size() << ")";
    for (const Edge &edge: edges()) {
        ss << "->" << edge.rc().getCode() << edge.rc().truncSize() << "(" << edge.getInnerId().eid << "|" << edge.getCoverage()<< "|" <<
           edge.rc().getInnerId().eid << ")" << edge.getCode() << edge.truncSize() << "->" << edge.getFinish().getInnerId() << "(" <<
           edge.getFinish().size() << ")";
    }
    ss << "]" << rightCut();
    return ss.str();
}

Segment<Edge> GraphPath::getSegment(PathPosition position) const {
    Edge &next = position.nextEdge();
    size_t seg_cut_left = firstPosition() == position ? leftCut() : 0;
    size_t seg_cut_right = lastPosition() == position + 1 ? rightCut() : 0;
    return {next, seg_cut_left, next.truncSize() - seg_cut_right};
}

void GraphPath::push_front(Edge &other) {
    VERIFY(cut_left == 0);
    VERIFY(!valid() || getStart() == other.getFinish());
    if(!valid())
        rc_start = other.getFinish().rc().getId();
    fsplits.push_front(other.getCode());
    rsplits.push_back(other.rc().getCode());
    start = other.getStart().getId();
}

void GraphPath::push_front(const Segment<Edge> &other) {
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

void GraphPath::push_front(const GraphPath &other) {
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

void GraphPath::normalize() {
    if(empty())
        return;
    while(frontEdge().isPrefix())
        pop_front();
    while(backEdge().isSuffix())
        pop_back();
}

void GraphPath::pop_back() {
    pop_back(backEdge());
}

void GraphPath::pop_back(Edge &edge) {
    VERIFY(edge.getFinish() == getFinish());
    VERIFY(fsplits.endsWith(edge.getCode()));
    VERIFY(rsplits.startsWith(edge.rc().getCode()));
    fsplits.pop_back(edge.getCode().size());
    rsplits.pop_front(edge.rc().getCode().size());
    size_t lost_len = edge.truncSize();
    cut_right -= std::min(cut_right, lost_len);
    rc_start = edge.getStart().rc().getId();
}

void GraphPath::pop_front() {
    pop_front(frontEdge());
}

void GraphPath::pop_front(Edge &edge) {
    VERIFY(edge.getStart() == getStart());
    VERIFY(fsplits.startsWith(edge.getCode()));
    VERIFY(rsplits.endsWith(edge.rc().getCode()));
    fsplits.pop_front(edge.getCode().size());
    rsplits.pop_back(edge.rc().getCode().size());
    size_t lost_len = edge.rc().truncSize();
    cut_left -= std::min(cut_left, lost_len);
    start = edge.getFinish().getId();
}


void GraphPath::setCutLeft(size_t value) {
    cut_left = value;
}

void GraphPath::setCutRight(size_t value) {
    cut_right = value;
}

GraphPath &GraphPath::extend(const Sequence &seq) {
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

GraphPath &GraphPath::fastExtend(const Sequence &seq) {
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

bool GraphPath::operator==(const GraphPath &other) const {
    return start == other.start && rc_start == other.rc_start &&
           cut_left == other.cut_left && cut_right == other.cut_right &&
           fsplits == other.fsplits && rsplits == other.rsplits;
}

GraphPath::GraphPath(const RAGraphPath &path, size_t cut_left, size_t cut_right) : GraphPath() { // NOLINT(google-explicit-constructor)
    for(Edge &e : path.edges()) {
        *this += e;
    }
    setCutLeft(cut_left);
    setCutRight(cut_right);
}

GraphPath::GraphPath(Vertex &start, const Sequence &extension) : GraphPath(start) {
    size_t pos = 0;
    while(pos < extension.size()) {
        Edge &next = getFinish().getOutgoing(extension[pos]);
        *this += next;
        pos += next.getCode().size();
    }
}

GraphPath GraphPath::Load(std::istream &os, const IdIndex<Vertex> &index) {
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

bool GraphPath::startsWith(const GraphPath &other) const {
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

bool GraphPath::endsWith(const GraphPath &other) const {
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

bool GraphPath::nonContradicts(const GraphPath &other) const {
    return startsWith(other) || other.startsWith(*this);
}

void GraphPath::forcePushBack(Edge &edge) {
    if(empty()) {
        start= edge.getStart().getId();
    }
    fsplits.push_back(edge.getCode());
    rsplits.push_front(edge.rc().getCode());
    rc_start = edge.getFinish().rc().getId();
}

void GraphPath::forcePushFront(Edge &edge) {
    if(empty()) {
        rc_start = edge.getFinish().rc().getId();
    }
    fsplits.push_front(edge.getCode());
    rsplits.push_back(edge.rc().getCode());
    start = edge.getStart().getId();
}

size_t GraphPath::calculateSize() const {
    size_t res = 0;
    for(Edge &edge : edges()) {
        res++;
    }
    return res;
}

void GraphPath::resetEdgeCodes() {
    if(!valid())
        return;
    NuclDeck newf;
    NuclDeck newr;
    for(Edge &edge : edges()) {
        if(!edge.isSuffix())
            newf.push_back(edge.firstNucl());
        if(!edge.isPrefix())
            newr.push_front(edge.rc().firstNucl());
    }
    fsplits = newf;
    rsplits = newr;
}

GraphPath GraphPath::operator*(size_t mult) const {
    VERIFY(getStart() == getFinish());
    GraphPath res(getStart());
    for(size_t i = 0; i < mult; i++)
        res += *this;
    return std::move(res);
}

RAGraphPath GraphPath::asRAPath() const {
    RAGraphPath res;
    for(Edge &edge : edges())
        res += edge;
    res.setCutLeft(leftCut());
    res.setCutRight(rightCut());
    return std::move(res);
}

void GraphPath::replace_back(Edge &edge) {
    Edge &old = backEdge();
    fsplits.replaceBack(old.getCode().size(), edge.getCode());
    rsplits.replaceFront(old.rc().getCode().size(), edge.rc().getCode());
    rc_start = edge.getFinish().rc().getId();
}

void GraphPath::replace_front(Edge &edge) {
    Edge &old = frontEdge();
    fsplits.replaceFront(old.getCode().size(), edge.getCode());
    rsplits.replaceBack(old.rc().getCode().size(), edge.rc().getCode());
    start = edge.getStart().getId();
}

bool GraphPath::isSingleton() const {
    return !empty() && fsplits.size() == frontEdge().getCode().size() && rsplits.size() == frontEdge().rc().getCode().size();
}

//    TODO: optimize to O(1)
GraphPath &GraphPath::shorten(PathPosition from, PathPosition to) {
    while(from != firstPosition())
        pop_front();
    while(to != lastPosition())
        pop_back();
    return *this;
}

GraphPath::GraphPath(Edge &edge, size_t cut_left, size_t cut_right) : GraphPath() { // NOLINT(google-explicit-constructor)
    *this += edge;
    this->cut_left = cut_left;
    this->cut_right = cut_right;
}
