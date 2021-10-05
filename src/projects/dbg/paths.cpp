#include "paths.hpp"

dbg::Path dbg::Path::WalkForward(dbg::Edge &start) {
    Path res(*start.start());
    res += start;
    Vertex *next = start.end();
    while(next != nullptr && next != start.start() && !next->isJunction()) {
        res += (*next)[0];
        next = (*next)[0].end();
    }
    return std::move(res);
}

dbg::Path dbg::Path::subPath(size_t from, size_t to) {
    if (from == to)
        return Path(getVertex(from));
    else
        return Path(getVertex(from), std::vector<Edge *>(path.begin() + from, path.begin() + to));
}

dbg::Path dbg::Path::RC() {
    std::vector<Edge *> rcPath;
    for (size_t i = path.size(); i > 0; i--) {
        rcPath.emplace_back(&path[i - 1]->rc());
    }
    return Path(back().end()->rc(), rcPath);
}

dbg::Vertex &dbg::Path::getVertex(size_t i) {
    VERIFY(i <= path.size());
    if (i == 0)
        return *start_;
    else
        return *path[i - 1]->end();
}

dbg::Vertex &dbg::Path::getVertex(size_t i) const {
    VERIFY(i <= path.size());
    if (i == 0)
        return *start_;
    else
        return *path[i - 1]->end();
}

size_t dbg::Path::find(dbg::Edge &edge, size_t pos) const {
    while(pos < size() && edge != *path[pos])
        pos++;
    if(pos == size())
        return -1;
    return pos;
}

size_t dbg::Path::find(dbg::Vertex &v, size_t pos) const {
    while(pos <= size() && v != getVertex(pos))
        pos++;
    if(pos > size())
        return -1;
    return pos;
}

double dbg::Path::minCoverage() const {
    double res = 100000;
    for (const Edge *edge : path) {
        res = std::min(edge->getCoverage(), res);
    }
    return res;
}

Sequence dbg::Path::Seq() const {
    SequenceBuilder sb;
    sb.append(start().seq);
    for (const Edge *e : path) {
        sb.append(e->seq);
    }
    return sb.BuildSequence();
}

Sequence dbg::Path::truncSeq() const {
    SequenceBuilder sb;
    for (const Edge *e : path) {
        sb.append(e->seq);
    }
    return sb.BuildSequence();
}

size_t dbg::Path::len() const {
    size_t res = 0;
    for (Edge *edge : path)
        res += edge->size();
    return res;
}

dbg::Path dbg::Path::operator+(const dbg::Path &other) const {
    VERIFY(finish() == *other.start_);
    std::vector<Edge *> edges = path;
    edges.insert(edges.end(), other.path.begin(), other.path.end());
    return {start(), std::move(edges)};
}

void dbg::Path::operator+=(dbg::Edge &edge) {
    path.emplace_back(&edge);
}

dbg::GraphAlignment dbg::GraphAlignment::RC() const {
    if(!valid())
        return {};
    GraphAlignment res(finish().rc());
    for (size_t i = 0; i < als.size(); i++) {
        res += als[als.size() - 1 - i].RC();
    }
    return res;
}

void dbg::GraphAlignment::invalidate() {
    start_ = nullptr;
    als.clear();
}

bool dbg::GraphAlignment::valid() const {
    return start_ != nullptr;
}

void dbg::GraphAlignment::cutBack(size_t l) {
    VERIFY(l <= len());
    while (size() > 0 && als.back().size() <= l) {
        l -= als.back().size();
        pop_back();
    }
    VERIFY(size() > 0);
    if (l > 0) {
        VERIFY(als.back().right > l);
        als.back().right -= l;
    }
}

dbg::GraphAlignment dbg::GraphAlignment::subalignment(size_t left, size_t right) const {
    if (left == right)
        return GraphAlignment(getVertex(left));
    else
        return {&getVertex(left), {als.begin() + left, als.begin() + right}};
}

dbg::GraphAlignment &dbg::GraphAlignment::addStep() {
    als.back().right += 1;
    return *this;
}

dbg::GraphAlignment &dbg::GraphAlignment::addStep(dbg::Edge &edge) {
    als.emplace_back(edge, 0, 1);
    return *this;
}

dbg::GraphAlignment &dbg::GraphAlignment::extend(const Sequence &seq) {
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
            if (als.back().contig().seq[als.back().right] == c) {
                addStep();
            } else {
                invalidate();
                return *this;
            }
        }
    }
    return *this;
}

bool dbg::GraphAlignment::endClosed() const {
    return start_ != nullptr && (als.size() == 0 || als.back().right == als.back().contig().size());
}

bool dbg::GraphAlignment::startClosed() const {
    return start_ != nullptr && (als.size() == 0 || als.front().left == 0);
}

unsigned char dbg::GraphAlignment::lastNucl() const {
    VERIFY(als.back().size() > 0);
    return als.back().contig().seq[als.back().right - 1];
}

size_t dbg::GraphAlignment::leftSkip() const {
    return als.size() == 0 ? 0 : als.front().left;
}

size_t dbg::GraphAlignment::rightSkip() const {
    return als.size() == 0 ? 0 : als.back().contig().size() - als.back().right;
}

std::vector<dbg::GraphAlignment> dbg::GraphAlignment::allSteps() {
    if (als.size() != 0 && als.back().right < als.back().contig().size()) {
        GraphAlignment copy = *this;
        return {std::move(copy.addStep())};
    }
    std::vector<GraphAlignment> res;
    Vertex &end = als.size() == 0 ? *start_ : *back().contig().end();
    for (Edge &edge : end) {
        GraphAlignment copy = *this;
        res.emplace_back(std::move(copy.addStep(edge)));
    }
    return res;
}

std::vector<dbg::GraphAlignment> dbg::GraphAlignment::allExtensions(size_t len) {
    std::vector<GraphAlignment> res = {*this};
    size_t left = 0;
    size_t right = 1;
    for (size_t l = 0; l < len; l++) {
        for (size_t i = left; i < right; i++) {
            std::vector<GraphAlignment> tmp = res[i].allSteps();
            res.insert(res.end(), tmp.begin(), tmp.end());
        }
        left = right;
        right = res.size();
    }
    return std::move(res);
}

Sequence dbg::GraphAlignment::map(std::unordered_map<const Edge *, Sequence> &edge_map) {
    SequenceBuilder sb;
    bool start = true;
    for (Segment<Edge> &seg : als) {
        auto it = edge_map.find(&seg.contig());
        if (it == edge_map.end()) {
            if (start) {
                sb.append((start_->seq + seg.contig().seq).Subseq(seg.left, seg.right + start_->seq.size()));
                start = false;
            } else {
                sb.append(seg.seq());
            }
        } else {
            size_t left = start_->seq.size();
            if (start) {
                left = 0;
            }
            size_t right = start_->seq.size();
            size_t sz = it->second.size() - start_->seq.size();
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

Sequence dbg::GraphAlignment::Seq() const {
    if (als.size() == 0) {
        return {};
    }
    SequenceBuilder sb;
    size_t k = start_->seq.size();
    if (als[0].left >= k)
        sb.append(als[0].contig().seq.Subseq(als[0].left - k, als[0].right));
    else {
        sb.append(start_->seq.Subseq(als[0].left, k));
        sb.append(als[0].contig().seq.Subseq(0, als[0].right));
    }
    for (size_t i = 1; i < als.size(); i++) {
        sb.append(als[i].seq());
    }
    return sb.BuildSequence();
}

Sequence dbg::GraphAlignment::truncSeq() const {
    SequenceBuilder sb;
    for (size_t i = 0; i < als.size(); i++) {
        sb.append(als[i].seq());
    }
    return sb.BuildSequence();
}

Sequence dbg::GraphAlignment::truncSeq(size_t start_position, size_t size) const {
    SequenceBuilder sb;
    size_t sz = 0;
    for (size_t i = start_position; i < als.size(); i++) {
//            std::cout << i << " " << sz << " " << size << std::endl;
//            std::cout << als[i].contig().size() << " " << als[i].left << " " << als[i].right << " " << als[i].size() << std::endl;
        if (sz + als[i].size() >= size) {
            sb.append(als[i].seq().Subseq(0, size - sz));
            break;
        } else {
            sb.append(als[i].seq());
            sz += als[i].size();
        }
    }
    return sb.BuildSequence();
}

dbg::GraphAlignment dbg::GraphAlignment::reroute(size_t left, size_t right, const dbg::GraphAlignment &rerouting) const {
    VERIFY(getVertex(left) == rerouting.start());
    VERIFY(getVertex(right) == rerouting.finish());
    return subalignment(0, left) + rerouting + subalignment(right, size());
}

dbg::GraphAlignment dbg::GraphAlignment::reroute(size_t left, size_t right, const dbg::Path &rerouting) const {
    VERIFY(getVertex(left) == rerouting.start());
    VERIFY(getVertex(right) == rerouting.finish());
    return subalignment(0, left) + rerouting + subalignment(right, size());
}

void dbg::GraphAlignment::operator+=(const dbg::Path &other) {
    if(other.size() == 0)
        return;
    if(!valid()) {
        start_ = &other.start();
    }
    VERIFY(finish() == other.getVertex(0));
    for (Edge *edge : other) {
        operator+=(Segment<Edge>(*edge, 0, edge->size()));
    }
}

void dbg::GraphAlignment::operator+=(const dbg::GraphAlignment &other) {
    if(other.size() == 0)
        return;
    if(!valid()) {
        start_ = &other.start();
    }
    VERIFY(finish() == other.getVertex(0));
    for (const Segment<Edge> &al : other) {
        operator+=(al);
    }
}

void dbg::GraphAlignment::operator+=(const Segment<Edge> &other) {
    if(!valid()) {
        start_ = other.contig().start();
    }
    if (!als.empty() && als.back().right < als.back().contig().size()) {
        als.back() = als.back() + other;
    } else {
        VERIFY(als.empty() || other.left == 0);
        VERIFY(finish() == *other.contig().start());
        als.push_back(other);
    }
}

void dbg::GraphAlignment::operator+=(Edge &other) {
    dbg::GraphAlignment::operator+=(Segment<Edge>(other, 0, other.size()));
}

dbg::GraphAlignment dbg::GraphAlignment::operator+(const dbg::GraphAlignment &other) const {
    GraphAlignment res = *this;
    res += other;
    return std::move(res);
}

dbg::GraphAlignment dbg::GraphAlignment::operator+(const dbg::Path &other) const {
    GraphAlignment res = *this;
    res += other;
    return std::move(res);
}

dbg::GraphAlignment dbg::GraphAlignment::operator+(const Segment<Edge> &other) const {
    GraphAlignment res = *this;
    res += other;
    return std::move(res);
}

dbg::GraphAlignment dbg::GraphAlignment::operator+(Edge &other) const {
    GraphAlignment res = *this;
    res += other;
    return std::move(res);
}

double dbg::GraphAlignment::minCoverage() const {
    double res = 100000;
    for (const Segment<Edge> &seg : als) {
        res = std::min(seg.contig().getCoverage(), res);
    }
    return res;
}

dbg::Path dbg::GraphAlignment::path() {
    std::vector<Edge *> res;
    for (auto &seg : als) {
        res.push_back(&seg.contig());
    }
    return {*start_, res};
}

std::string dbg::GraphAlignment::str(bool show_coverage) const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << leftSkip() << " " << start().getId();
    for(const Segment<Edge> &seg : als) {
        ss << " " << seg.size() << "ACGT"[seg.contig().seq[0]];
        if(show_coverage) {
            ss << "(" << seg.contig().getCoverage() << ")";
        }
        ss << " " << seg.contig().end()->getId();
    }
    ss << " " << rightSkip();
    return ss.str();
}

size_t dbg::GraphAlignment::len() const {
    size_t res = 0;
    for (auto &seg : als) {
        res += seg.size();
    }
    return res;
}

size_t dbg::GraphAlignment::find(dbg::Edge &edge, size_t pos) const {
    while(pos < size() && edge != als[pos].contig())
        pos++;
    if(pos == size())
        return -1;
    return pos;
}

size_t dbg::GraphAlignment::find(dbg::Vertex &v, size_t pos) const {
    while(pos <= size() && v != getVertex(pos))
        pos++;
    if(pos > size())
        return -1;
    return pos;
}

dbg::GraphAlignment dbg::GraphAligner::align(const Sequence &seq, dbg::Edge *edge_to, size_t pos_to) {
    size_t k = dbg.hasher().getK();
    size_t cur = k;
    GraphAlignment res;
    while(cur < seq.size()) {
        size_t len = std::min(seq.size() - cur, edge_to->size() - pos_to);
        res += Segment<Edge>(*edge_to, pos_to, pos_to + len);
        cur += len;
        if(cur < seq.size()) {
            edge_to = &edge_to->end()->getOutgoing(seq[cur]);
            pos_to = 0;
        }
    }
    return res;
}

dbg::GraphAlignment dbg::GraphAligner::align(const Sequence &seq) const {
    std::vector<hashing::KWH> kmers = dbg.extractVertexPositions(seq, 1);
    size_t k = dbg.hasher().getK();
    GraphAlignment res;
    if (kmers.size() == 0) {
        hashing::KWH kwh(dbg.hasher(), seq, 0);
        while (true) {
            if (dbg.isAnchor(kwh.hash())) {
                EdgePosition pos = dbg.getAnchor(kwh);
                VERIFY(kwh.pos < pos.pos);
                VERIFY(pos.pos + seq.size() - kwh.pos <= pos.edge->size() + k);
                Segment<Edge> seg(*pos.edge, pos.pos - kwh.pos, pos.pos + seq.size() - kwh.pos - k);
                return {pos.edge->start(), std::vector<Segment<Edge>>({seg})};
            }
            if (!kwh.hasNext()) {
#pragma omp critical
                {
                    std::cout << "Error: could not align sequence " << seq.size() << std::endl;
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
            std::cout << "No outgoing for start" << std::endl << seq << std::endl <<
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
            std::cout << "No outgoing for middle\n" << seq << "\n" << cpos << " " << prestart->getId() <<
                " " << prestart->outDeg() << "\n" << seq.Subseq(cpos -k, cpos) << "\n" << prestart->seq << "\n" <<
                size_t(seq[cpos]) << std::endl;
            for(Edge &tmp : *prestart) {
                std::cout << tmp.getId() << " " << tmp.size() << std::endl;
            }
            VERIFY(false);
        }
        Edge &next = prestart->getOutgoing(seq[cpos]);
        size_t len = std::min<size_t>(next.size(), seq.size() - cpos);
        res += Segment<Edge>(next, 0, len);
        cpos += len;
        prestart = next.end();
    }
    return std::move(res);
}

dbg::GraphAlignment dbg::GraphAligner::align(const dbg::EdgePosition &pos, const Sequence &seq) const {
    GraphAlignment res(pos.edge->start(), {{*pos.edge, pos.pos, pos.pos}});
    Edge *cedge = pos.edge;
    size_t epos = pos.pos;
    for (size_t cpos = 0; cpos < seq.size(); cpos++) {
        unsigned char c = seq[cpos];
        if (epos == cedge->size()) {
            Vertex &v = *cedge->end();
            if (v.hasOutgoing(c)) {
                cedge = &v.getOutgoing(c);
                res.addStep(*cedge);
                epos = 1;
            } else {
                return {};
            }
        } else {
            if (cedge->seq[epos] == c) {
                res.addStep();
                epos += 1;
            } else {
                return {};
            }
        }
    }
    return std::move(res);
}

std::vector<dbg::PerfectAlignment<dbg::Edge, dbg::Edge>>
dbg::GraphAligner::oldEdgeAlign(dbg::Edge
                                &contig) const {
    Sequence seq = contig.start()->seq + contig.seq;
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
                if (gpos.edge->seq[gpos.pos] == seq[kwh.pos + k]) {
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
    Sequence seq = contig.seq;
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
                    while (len < edge.size() && len < kwh.pos && edge.seq[len] == (seq[kwh.pos - len - 1] ^ 3))
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos - len, kwh.pos),
                                     Segment<Edge>(edge.rc(), edge.size() - len, edge.size()));
                }
                if (kwh.pos + k < seq.size() && vertex.hasOutgoing(seq[kwh.pos + k])) {
                    Edge &edge = vertex.getOutgoing(seq[kwh.pos + k]);
                    size_t len = 1;
                    while (len < edge.size() && kwh.pos + k + len < seq.size() &&
                           edge.seq[len] == seq[kwh.pos + k + len])
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos, kwh.pos + len),
                                     Segment<Edge>(edge, 0, len));
                }
            } else if ((res.empty() || kwh.pos > res.back().seg_from.right) && dbg.isAnchor(kwh.hash())) {
                EdgePosition pos = dbg.getAnchor(kwh);
//                TODO replace this code with a call to expand method of PerfectAlignment class after each edge is marked by its full sequence
                Edge &edge = *pos.edge;
                Vertex &start = *pos.edge->start();
                CompositeSequence edge_seq({start.seq, edge.seq});
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
