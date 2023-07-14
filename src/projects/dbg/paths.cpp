#include "paths.hpp"
using namespace dbg;
GraphPath GraphPath::WalkForward(dbg::Edge &start) {
    GraphPath res(start);
    Vertex *next = &start.getFinish();
    VERIFY(next != nullptr);
    while(*next != start.getStart() && !next->isJunction()) {
        VERIFY(next != nullptr);
        res += next->front();
        next = &res.finish();
    }
    return std::move(res);
}

dbg::GraphPath dbg::GraphPath::subPath(size_t from, size_t to) {
    if (from == to)
        return GraphPath(getVertex(from));
    else
        return {getVertex(from), std::vector<Edge *>(path.begin() + from, path.begin() + to),
                from == 0 ? leftSkip() : 0, to == size() - 1 ? rightSkip() : 0};
}

dbg::Vertex &dbg::GraphPath::getVertex(size_t i) {
    VERIFY(i <= path.size());
    if (i == 0)
        return *start_;
    else
        return path[i - 1]->getFinish();
}

dbg::Vertex &dbg::GraphPath::getVertex(size_t i) const {
    VERIFY(i <= path.size());
    if (i == 0)
        return *start_;
    else
        return path[i - 1]->getFinish();
}

size_t dbg::GraphPath::find(dbg::Edge &edge, size_t pos) const {
    while(pos < size() && edge != *path[pos])
        pos++;
    if(pos == size())
        return -1;
    return pos;
}

size_t dbg::GraphPath::find(dbg::Vertex &v, size_t pos) const {
    while(pos <= size() && v != getVertex(pos))
        pos++;
    if(pos > size())
        return -1;
    return pos;
}

double dbg::GraphPath::minCoverage() const {
    double res = 100000000;
    for (const Edge *edge : path) {
        res = std::min(edge->getCoverage(), res);
    }
    return res;
}

size_t dbg::GraphPath::len() const {
    size_t res = 0;
    for (Edge *edge : path)
        res += edge->size();
    return res - cut_left - cut_right;
}

IterableStorage<GraphPath::vertex_iterator> GraphPath::vertices() const & {
    std::function<Vertex &(size_t)> transformer = [this](size_t ind)->Vertex &{return getVertex(ind);};
    CountingIterator<size_t> end_it = CountingIterator<size_t>(size() + 1);
    vertex_iterator vbegin (CountingIterator<size_t>(0), end_it, transformer);
    vertex_iterator vend(end_it, end_it, transformer);
    return {vbegin, vend};
}

//IterableStorage<GraphPath::const_vertex_iterator> GraphPath::vertices() const {
//    std::function<const Vertex &(size_t)> transformer = [this](size_t ind)->const Vertex &{return getVertex(ind);};
//    const_vertex_iterator vbegin (CountingIterator<size_t>(0), CountingIterator<size_t>(size() + 1), transformer);
//    const_vertex_iterator vend(CountingIterator<size_t>(size() + 1), CountingIterator<size_t>(size() + 1), transformer);
//    return {vbegin, vend};
//}

IterableStorage<GraphPath::edge_iterator> GraphPath::edges() const &{
    std::function<Edge &(size_t)> transformer = [this](size_t ind)->Edge &{return *path[ind];};
    CountingIterator<size_t> end_it = CountingIterator<size_t>(size());
    edge_iterator ebegin (CountingIterator<size_t>(0), end_it, transformer);
    edge_iterator eend(end_it, end_it, transformer);
    return {ebegin, eend};
}

//IterableStorage<GraphPath::const_edge_iterator> GraphPath::edges() const {
//    std::function<const Edge &(size_t)> transformer = [this](size_t ind)->const Edge &{return *path[ind];};
//    const_edge_iterator ebegin (CountingIterator<size_t>(0), CountingIterator<size_t>(size()), transformer);
//    const_edge_iterator eend(CountingIterator<size_t>(size() + 1), CountingIterator<size_t>(size()), transformer);
//    return {ebegin, eend};
//}

//IterableStorage<GraphPath::const_segment_iterator> GraphPath::segments() const {
//    std::function<Segment<const Edge> (size_t)> transformer = [this](size_t ind)->Segment<const Edge> {
//        Segment<const Edge> res(*path[ind]);
//        if(ind == 0)
//            res = res.shrinkLeftBy(cut_left);
//        if(ind + 1 == size()) {
//            res = res.shrinkRightBy(cut_right);
//        }
//        return res;
//    };
//    const_segment_iterator sbegin (CountingIterator<size_t>(0), CountingIterator<size_t>(size() + 1), transformer);
//    const_segment_iterator send(CountingIterator<size_t>(size() + 1), CountingIterator<size_t>(size() + 1), transformer);
//    return {sbegin, send};
//}
//
dbg::GraphPath dbg::GraphPath::RC() const {
    if(!valid())
        return {};
    std::vector<Edge *> res;
    for(auto it  = path.rbegin(); it != path.rend(); ++it) {
        Edge *e = *it;
        res.emplace_back(&e->rc());
    }
    return {finish().rc(), res, rightSkip(), leftSkip()};
}

void dbg::GraphPath::invalidate() {
    start_ = nullptr;
    path.clear();
    cut_left = 0;
    cut_right = 0;
}

bool dbg::GraphPath::valid() const {
    VERIFY(start_ != nullptr || size() == 0);
    return start_ != nullptr;
}

void dbg::GraphPath::cutBack(size_t l) {
    VERIFY(l <= len());
    while(l >= back().size()) {
        l -= back().size();
        pop_back();
        cut_right = 0;
    }
    VERIFY(size() > 0);
    cut_right += l;
}

void dbg::GraphPath::cutFront(size_t l) {
    VERIFY(l <= len());
    size_t cut = 0;
    while(cut < size() && l >= operator[](cut).size()) {
        l -= operator[](cut).size();
        cut++;
        cut_left = 0;
    }
    path.erase(path.begin(), path.begin() + cut);
    cut_left += l;
}

dbg::GraphPath dbg::GraphPath::subalignment(size_t left, size_t right) const {
    if (left == right)
        return {getVertex(left)};
    else {
        size_t left_skip = left == 0 ? leftSkip() : 0;
        size_t right_skip = right == size() ? rightSkip() : 0;
        return {getVertex(left), {path.begin() + left, path.begin() + right}, left_skip, right_skip};
    }
}

dbg::GraphPath &dbg::GraphPath::addStep() {
    cut_right -= 1;
    return *this;
}

dbg::GraphPath &dbg::GraphPath::addStep(dbg::Edge &edge) {
    *this += Segment<Edge>(edge, 0, 1);
    return *this;
}

dbg::GraphPath &dbg::GraphPath::extend(const Sequence &seq) {
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

bool dbg::GraphPath::endClosed() const {
    return valid() && rightSkip() == 0;
}

bool dbg::GraphPath::startClosed() const {
    return valid() && leftSkip() == 0;
}

unsigned char dbg::GraphPath::lastNucl() const {
    Segment<Edge> seg = back();
    return seg.truncSeq()[seg.right - 1];
}

size_t dbg::GraphPath::leftSkip() const {
    return cut_left;
}

size_t dbg::GraphPath::rightSkip() const {
    return cut_right;
}

std::vector<dbg::GraphPath> dbg::GraphPath::allSteps() {
    if (size() != 0 && cut_right > 0) {
        GraphPath copy = *this;
        return {std::move(copy.addStep())};
    }
    std::vector<GraphPath> res;
    Vertex &end = finish();
    for (Edge &edge : end) {
        GraphPath copy = *this;
        res.emplace_back(std::move(copy.addStep(edge)));
    }
    return res;
}

std::vector<dbg::GraphPath> dbg::GraphPath::allExtensions(size_t len) {
    std::vector<GraphPath> res = {*this};
    size_t left = 0;
    size_t right = 1;
    for (size_t l = 0; l < len; l++) {
        for (size_t i = left; i < right; i++) {
            std::vector<GraphPath> tmp = res[i].allSteps();
            res.insert(res.end(), tmp.begin(), tmp.end());
        }
        left = right;
        right = res.size();
    }
    return std::move(res);
}

Sequence dbg::GraphPath::map(std::unordered_map<const Edge *, Sequence> &edge_map) {
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

Sequence dbg::GraphPath::Seq() const {
    if (size() == 0)
        return {};
    SequenceBuilder sb;
    sb.append(front().fullSeq());
    for(size_t i = 1; i < size(); i++) {
        sb.append(operator[](i).truncSeq());
    }
    return sb.BuildSequence();
}

Sequence dbg::GraphPath::truncSeq() const {
    SequenceBuilder sb;
    for (Segment<Edge> seg : *this) {
        sb.append(seg.truncSeq());
    }
    return sb.BuildSequence();
}

Sequence dbg::GraphPath::truncSubseq(size_t start_position, size_t sz) const {
    SequenceBuilder sb;
    for (size_t i = start_position; i < this->size(); i++) {
//            std::cout << i << " " << sz << " " << size << std::endl;
//            std::cout << als[i].contig().size() << " " << als[i].left << " " << als[i].right << " " << als[i].size() << std::endl;
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

dbg::GraphPath dbg::GraphPath::reroute(size_t left, size_t right, const dbg::GraphPath &rerouting) const {
    VERIFY(getVertex(left) == rerouting.start());
    VERIFY(getVertex(right) == rerouting.finish());
    return subalignment(0, left) + rerouting + subalignment(right, size());
}

void dbg::GraphPath::operator+=(const dbg::GraphPath &other) {
    if(other.size() == 0)
        return;
    if(!valid()) {
        *this = other;
        return;
    }
    VERIFY(finish() == other.getVertex(0));
    for (Segment<Edge> al : other) {
        operator+=(al);
    }
}

void dbg::GraphPath::operator+=(const Segment<Edge> &other) {
    if(!valid()) {
        *this = {other};
        return;
    }
    if(cut_right == 0) {
        VERIFY(other.left == 0 && finish() == other.contig().getStart());
        path.emplace_back(&other.contig());
        cut_right = other.contig().size() - other.right;
    } else {
        VERIFY(cut_right == other.contig().size() - other.left && other.contig() == backEdge());
    }
    cut_right = other.contig().size() - other.right;
}

void dbg::GraphPath::operator+=(Edge &other) {
    dbg::GraphPath::operator+=(Segment<Edge>(other, 0, other.size()));
}

dbg::GraphPath dbg::GraphPath::operator+(const dbg::GraphPath &other) const {
    GraphPath res = *this;
    res += other;
    return std::move(res);
}

dbg::GraphPath dbg::GraphPath::operator+(const Segment<Edge> &other) const {
    GraphPath res = *this;
    res += other;
    return std::move(res);
}

dbg::GraphPath dbg::GraphPath::operator+(Edge &other) const {
    GraphPath res = *this;
    res += other;
    return std::move(res);
}

std::string dbg::GraphPath::str(bool show_coverage) const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << leftSkip() << " " << start().getInnerId();
    for(const Segment<Edge> &seg : *this) {
        ss << " " << seg.size() << "/" <<seg.contig().size() << "ACGT"[seg.contig().truncSeq()[0]] ;
        if(show_coverage) {
            ss << "(" << seg.contig().getCoverage() << ")";
        }
        ss << " " << seg.contig().getFinish().getInnerId();
    }
    ss << " " << rightSkip();
    return ss.str();
}

Segment<Edge> GraphPath::back() const {
    return {*path.back(), (size() == 1 ? leftSkip() : 0), path.back()->size() - rightSkip()};
}

Segment<Edge> GraphPath::front() const {
    return {*path.front(), leftSkip(), size() == 1 ? path.front()->size() - rightSkip() : path.front()->size()};
}

Segment<Edge> GraphPath::operator[](size_t i) const {
    return {*path[i], i == 0 ? leftSkip() : 0, i == size() - 1 ? path.back()->size() - rightSkip() : path[i]->size()};
}

GraphPath::segment_iterator GraphPath::begin() const {
    std::function<Segment<Edge> (size_t)> transformer = [this](size_t ind)->Segment<Edge> {
        return operator[](ind);
    };
    CountingIterator<size_t> end_it(size());
    return {CountingIterator<size_t>(0), end_it, transformer};
}

GraphPath::segment_iterator GraphPath::end() const {
    std::function<Segment<Edge> (size_t)> transformer = [this](size_t ind)->Segment<Edge> {
        return operator[](ind);
    };
    CountingIterator<size_t> end_it(size());
    return {end_it, end_it, transformer};
}

dbg::GraphPath dbg::GraphAligner::align(const Sequence &seq, dbg::Edge *edge_to, size_t pos_to) {
    size_t k = dbg.hasher().getK();
    size_t cur = k;
    GraphPath res;
    while(cur < seq.size()) {
        size_t len = std::min(seq.size() - cur, edge_to->size() - pos_to);
        res += Segment<Edge>(*edge_to, pos_to, pos_to + len);
        cur += len;
        if(cur < seq.size()) {
            edge_to = &edge_to->getFinish().getOutgoing(seq[cur]);
            pos_to = 0;
        }
    }
    return res;
}

dbg::GraphPath dbg::GraphAligner::align(const Sequence &seq, const std::string &name) const {
    std::vector<hashing::KWH> kmers = dbg.extractVertexPositions(seq, 1);
    size_t k = dbg.hasher().getK();
    GraphPath res;
    if (kmers.size() == 0) {
        hashing::KWH kwh(dbg.hasher(), seq, 0);
        while (true) {
            if (dbg.isAnchor(kwh.hash())) {
                EdgePosition pos = dbg.getAnchor(kwh);
                VERIFY(kwh.pos < pos.pos);
                VERIFY(pos.pos + seq.size() - kwh.pos <= pos.edge->size() + k);
                Segment<Edge> seg(*pos.edge, pos.pos - kwh.pos, pos.pos + seq.size() - kwh.pos - k);
                return {seg};
            }
            if (!kwh.hasNext()) {
#pragma omp critical
                {
                    std::cout << "Error: could not align sequence " << seq.size() << " " << name << std::endl;
                    std::cout << seq << std::endl;
                    abort();
                };
                return res;
            }
            kwh = kwh.next();
        }
    }
    Vertex *prestart = &dbg.getVertex(kmers.front());
    if (kmers.front().pos > 0) {
        Vertex &rcstart = prestart->rc();
        if (!rcstart.hasOutgoing(seq[kmers.front().pos - 1] ^ 3)) {
            std::cout << "No outgoing for getStart" << std::endl << seq << std::endl <<
                      kmers.front().pos << " " << seq[kmers.front().pos - 1] << std::endl
                      << kmers.front().getSeq() << std::endl;
            VERIFY(false);
        }
        Edge &rcedge = rcstart.getOutgoing(seq[kmers.front().pos - 1] ^ 3);
        Edge &edge = rcedge.rc();
        VERIFY(edge.size() >= kmers.front().pos);
        Segment<Edge> seg(edge, edge.size() - kmers.front().pos, edge.size());
        res += seg;
    }
    size_t cpos = kmers.front().pos + k;
    while(cpos < seq.size()) {
        if(!prestart->hasOutgoing(seq[cpos])) {
            std::cout << "No outgoing for middle\n" << seq << "\n" << cpos << " " << prestart->getInnerId() <<
                      " " << prestart->outDeg() << "\n" << seq.Subseq(cpos -k, cpos) << "\n" << prestart->getSeq() << "\n" <<
                      size_t(seq[cpos]) << std::endl;
            for(Edge &tmp : *prestart) {
                std::cout << tmp.getInnerId() << " " << tmp.size() << std::endl;
            }
            VERIFY(false);
        }
        Edge &next = prestart->getOutgoing(seq[cpos]);
        size_t len = std::min<size_t>(next.size(), seq.size() - cpos);
        res += Segment<Edge>(next, 0, len);
        cpos += len;
        prestart = &next.getFinish();
    }
    return std::move(res);
}

dbg::GraphPath dbg::GraphAligner::align(const dbg::EdgePosition &pos, const Sequence &seq) const {
    GraphPath res(Segment<Edge>(*pos.edge, pos.pos, pos.pos));
    Edge *cedge = pos.edge;
    size_t epos = pos.pos;
    for (size_t cpos = 0; cpos < seq.size(); cpos++) {
        unsigned char c = seq[cpos];
        if (epos == cedge->size()) {
            Vertex &v = cedge->getFinish();
            if (v.hasOutgoing(c)) {
                cedge = &v.getOutgoing(c);
                res.addStep(*cedge);
                epos = 1;
            } else {
                return {};
            }
        } else {
            if (cedge->truncSeq()[epos] == c) {
                res.addStep();
                epos += 1;
            } else {
                return {};
            }
        }
    }
    return std::move(res);
}

std::vector<dbg::PerfectAlignment<dbg::Edge, dbg::Edge>> dbg::GraphAligner::oldEdgeAlign(dbg::Edge &contig) const {
    Sequence seq = contig.getSeq();
    std::vector<PerfectAlignment<Edge, Edge>> res;
    hashing::KWH kwh(dbg.hasher(), seq, 0);
    size_t k = dbg.hasher().getK();
    while (true) {
        if (!kwh.hasNext())
            break;
        if (res.empty() || kwh.pos >= res.back().seg_from.right) {
            Edge *edge = nullptr;
            size_t pos = 0;
            if (dbg.containsVertex(kwh.hash())) {
                Vertex &start = dbg.getVertex(kwh);
                if (start.hasOutgoing(seq[kwh.pos + k]))
                    edge = &dbg.getVertex(kwh).getOutgoing(seq[kwh.pos + k]);
            }
            if (edge == nullptr && dbg.isAnchor(kwh.hash())) {
                EdgePosition gpos = dbg.getAnchor(kwh);
                if (gpos.edge->truncSeq()[gpos.pos] == seq[kwh.pos + k]) {
                    edge = gpos.edge;
                    pos = gpos.pos;
                }
            }
            if (edge != nullptr) {
                size_t len = std::min(contig.size() - kwh.pos, edge->size() - pos);
                res.emplace_back(Segment<Edge>(contig, kwh.pos, kwh.pos + len), Segment<Edge>(*edge, pos, pos + len));
            }
        }
        kwh = kwh.next();
    }
    return std::move(res);
}

std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> dbg::GraphAligner::carefulAlign(Contig &contig) const {
    Sequence seq = contig.getSeq();
    size_t k = dbg.hasher().getK();
    if(contig.size() < k) {
        return {};
    }
    std::vector<PerfectAlignment<Contig, Edge>> res;
    hashing::KWH kwh(dbg.hasher(), seq, 0);
    while (true) {
        if (res.empty() || kwh.pos >= res.back().seg_from.right) {
            if (dbg.containsVertex(kwh.hash())) {
                Vertex &vertex = dbg.getVertex(kwh);
                Vertex &rcVertex = vertex.rc();
                if ((res.empty() || kwh.pos > res.back().seg_from.right)
                    && kwh.pos > 0 && rcVertex.hasOutgoing(seq[kwh.pos - 1] ^ 3)) {
                    Edge &edge = rcVertex.getOutgoing(seq[kwh.pos - 1] ^ 3);
                    size_t len = 1;
                    while (len < edge.size() && len < kwh.pos && edge.truncSeq()[len] == (seq[kwh.pos - len - 1] ^ 3))
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos - len, kwh.pos),
                                     Segment<Edge>(edge.rc(), edge.size() - len, edge.size()));
                }
                if (kwh.pos + k < seq.size() && vertex.hasOutgoing(seq[kwh.pos + k])) {
                    Edge &edge = vertex.getOutgoing(seq[kwh.pos + k]);
                    size_t len = 1;
                    while (len < edge.size() && kwh.pos + k + len < seq.size() &&
                            edge.truncSeq()[len] == seq[kwh.pos + k + len])
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos, kwh.pos + len),
                                     Segment<Edge>(edge, 0, len));
                }
            } else if ((res.empty() || kwh.pos > res.back().seg_from.right) && dbg.isAnchor(kwh.hash())) {
                EdgePosition pos = dbg.getAnchor(kwh);
//                TODO replace this code with a call to expand method of PerfectAlignment class after each edge is marked by its full sequence
                Edge &edge = *pos.edge;
                Vertex &start = pos.edge->getStart();
                CompositeSequence edge_seq({start.getSeq(), edge.truncSeq()});
                size_t left_from = kwh.pos;
                size_t right_from = kwh.pos + k;
                size_t left_to = pos.pos;
                size_t right_to = pos.pos + k;
                while (left_from > 0 && left_to > 0 && edge_seq[left_to - 1] == seq[left_from - 1]) {
                    left_from -= 1;
                    left_to -= 1;
                }
                while (right_from < seq.size() && right_to < edge_seq.size() &&
                       seq[right_from] == edge_seq[right_to]) {
                    right_from += 1;
                    right_to += 1;
                }
                if (left_to - left_from > k) {
                    res.emplace_back(Segment<Contig>(contig, left_from, right_from - k),
                                     Segment<Edge>(edge, left_to, right_to - k));
                }
            }
        }
        if (!kwh.hasNext())
            break;
        kwh = kwh.next();
    }
    return std::move(res);
}

PerfectAlignment<Contig, dbg::Edge> bestExtension(const Vertex &vertex, const Segment<Contig> &seg) {
    PerfectAlignment<Contig, dbg::Edge> best({seg.contig(), seg.left, seg.left}, {});
    for(Edge &edge : vertex) {
        size_t len = 0;
        while(len < edge.size() && seg.left + vertex.size() + len < seg.contig().size()) {
            if(seg.contig()[seg.left + vertex.size() + len] != edge.truncSeq()[len])
                break;
            len++;
        }
//        std::cout << len << std::endl;
//        std::cout << Segment<Contig>(seg.contig(), seg.left + vertex.seq.size(),
//                                     std::min(seg.contig().size(), seg.left + vertex.seq.size() + 200)).seq() << std::endl;
//        std::cout << edge.seq.Subseq(0, std::min<size_t>(edge.seq.size(), 200)) << std::endl;
        if(len >= best.size()) {
            best = {{seg.contig(), seg.left, seg.left + len}, {edge, 0, len}};
        }
    }
    return best;
}

PerfectAlignment<Contig, dbg::Edge> bestExtension(Edge &edge, const Segment<Contig> &seg) {
    size_t len = 0;
    while(len < edge.size() && seg.left + edge.getStart().getSeq().size() + len < seg.contig().size()) {
        if(seg.contig()[seg.left + edge.getStart().getSeq().size() + len] != edge.truncSeq()[len])
            break;
        len++;
    }
    return {Segment<Contig>(seg.contig(), seg.left, seg.left + len), Segment<Edge>(edge, 0, len)};
}

PerfectAlignment<Contig, dbg::Edge> GraphAligner::extendLeft(const hashing::KWH &kwh, Contig &contig) const {
    size_t k = dbg.hasher().getK();
    PerfectAlignment<Contig, dbg::Edge> best({contig, kwh.pos, kwh.pos}, {});
    if(kwh.pos == 0) {
        return best;
    }
    Vertex &start = dbg.getVertex(kwh);
    Contig rc_contig = contig.RC();
    if(start.inDeg() == 0) {
        return best;
    }
    PerfectAlignment<Contig, Edge> start_al = bestExtension(start.rc(), Segment<Contig>(rc_contig, contig.size() - k - kwh.pos, contig.size() - k));
    return {Segment<Contig>(contig, contig.size() - k - start_al.seg_from.right, contig.size() - k - start_al.seg_from.left),
            Segment<Edge>(start_al.seg_to.contig().rc(),
                          start_al.seg_to.contig().size() - start_al.seg_to.right, start_al.seg_to.contig().size() - start_al.seg_to.left)};
}

PerfectAlignment<Contig, dbg::Edge> GraphAligner::extendRight(const hashing::KWH &kwh, Contig &contig) const {
    size_t k = dbg.hasher().getK();
    PerfectAlignment<Contig, dbg::Edge> best({contig, kwh.pos, kwh.pos}, {});
    if(kwh.pos + k == contig.size()) {
        return best;
    }
    Vertex &start = dbg.getVertex(kwh);
    if(start.outDeg() == 0) {
        return best;
    }
    return bestExtension(start, Segment<Contig>(contig, kwh.pos, contig.size() - k));
}

std::vector<PerfectAlignment<Contig, dbg::Edge>> GraphAligner::sparseAlign(Contig &contig) const {
    std::vector<hashing::KWH> vlist = dbg.extractVertexPositions(contig.getSeq());
    std::vector<PerfectAlignment<Contig, dbg::Edge>> result;
    size_t k = dbg.hasher().getK();
    if(vlist.empty())
        return result;
    for(hashing::KWH &kwh : vlist) {
        if(result.empty() || result.back().seg_from.right != kwh.pos) {
            PerfectAlignment<Contig, Edge> new_al = extendLeft(kwh, contig);
            if(new_al.size() > 0) {
                result.emplace_back(new_al);
            }
        }
        {
            PerfectAlignment<Contig, Edge> new_al = extendRight(kwh, contig);
            if(new_al.size() > 0) {
                result.emplace_back(new_al);
            }
        }
    }
    return std::move(result);
}