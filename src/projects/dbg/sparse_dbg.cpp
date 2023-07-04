#include "sparse_dbg.hpp"
using namespace dbg;

Edge Edge::_fake = Edge();

bool IsMarkerCorrect(EdgeMarker marker) {
    return marker == EdgeMarker::correct || marker == EdgeMarker::unique || marker == EdgeMarker::repeat;
}

size_t Edge::updateTipSize() const {
    size_t new_val = 0;
    if(extraInfo == size_t(-1) && finish->inDeg() == 1) {
        for (const Edge & other : *finish) {
            other.finish->lock();
            new_val = std::max(new_val, other.extraInfo);
            other.finish->unlock();
        }
        if(new_val != size_t(-1))
            new_val += size();
        finish->lock();
        extraInfo = new_val;
        finish->unlock();
    }
    return new_val;
}

Edge &Edge::rc() const {
    VERIFY(!getStart()->getSeq().empty());
    Vertex &vend = finish->rc();
    unsigned char c;
    size_t k = vend.getSeq().size();
    if(size() > k) {
        c = (!getSeq())[k];
    } else {
        c = (!start->getSeq())[k - size()];
    }
    return vend.getOutgoing(c);
}

Edge &Edge::sparseRcEdge() const {
    Vertex &vend = getFinish()->rc();
    VERIFY(!start->getSeq().empty());
    for(Edge &candidate : vend) {
        if(*candidate.getFinish() == start->rc() && candidate.size() == size() &&
           (size() <= start->getSeq().size() || candidate.getSeq().startsWith((!getSeq()).Subseq(start->getSeq().size())))) {
            return candidate;
        }
    }
    std::cout << start->getSeq() + getSeq() << std::endl;
    std::cout << getSeq() << std::endl;
    std::cout << vend.getSeq() << std::endl;
    for(Edge &candidate : vend) {
        if(*candidate.getFinish() == start->rc() && candidate.size() == size() &&
           (size() <= getSeq().size() || candidate.getSeq().startsWith((!getSeq()).Subseq(start->getSeq().size())))) {
            std::cout << vend.getSeq() + candidate.getSeq() << std::endl;
        }
    }
    VERIFY(false);
    return vend.front();
}

void Edge::bindTip(Vertex &start, Vertex &end) {
    VERIFY(finish == nullptr);
    finish = &end;
    Sequence rcseq = !(start.getSeq() + getSeq());
    end.rc().addEdgeLockFree(Edge(end.rc(), start.rc(), rcseq.Subseq(start.getSeq().size())));
}

size_t Edge::getTipSize() const {
    return extraInfo;
}

void Edge::setTipSize(size_t val) const {
    extraInfo = val;
}

Vertex *Edge::getStart() const {
    return start;
}

Vertex *Edge::getFinish() const {
    return finish;
}

size_t Edge::size() const {
    return getSeq().size();
}

double Edge::getCoverage() const {
    VERIFY(size() != 0);
    return double(cov) / size();
}

size_t Edge::intCov() const {
    return cov;
}

void Edge::incCov(size_t val) const {
#pragma omp atomic
    cov += val;
}

bool Edge::operator==(const Edge &other) const {
    return this == &other;
}

bool Edge::operator!=(const Edge &other) const {
    return this != &other;
}

bool Edge::operator<(const Edge &other) const {
    if(this == &other)
        return false;
    if(start != other.start)
        return *start < *other.start;
    return this->getSeq() < other.getSeq();
}

bool Edge::operator>(const Edge &other) const {
    if(this == &other)
        return false;
    if(start != other.start)
        return *start > *other.start;
    return other.getSeq() < getSeq();
}

bool Edge::operator<=(const Edge &other) const {
    return *this == other || *this < other;
}

std::string Edge::str() const {
    std::stringstream ss;
    const dbg::Vertex &v = *getStart();
    ss << v.hash() << v.isCanonical() << "ACGT"[getSeq()[0]];
    return ss.str();
}

Sequence Edge::kmerSeq(size_t pos) const {
    VERIFY(pos <= getSeq().size());
    size_t k = start->getSeq().size();
    if (pos >= k)
        return getSeq().Subseq(pos - k, pos);
    else {
        return getStart()->getSeq().Subseq(pos) + getSeq().Subseq(0, pos);
    }
}

Sequence Edge::suffix(size_t pos) const {
    VERIFY(pos <= getSeq().size());
    size_t k = start->getSeq().size();
    if (pos >= k)
        return getSeq().Subseq(pos - k, getSeq().size());
    else {
        return getStart()->getSeq().Subseq(pos) + getSeq().Subseq(0, getSeq().size());
    }
}

std::string Edge::getId() const {
    return start->getId() + "ACGT"[getSeq()[0]];
}

std::string Edge::oldId() const {
    return start->oldId() + "ACGT"[getSeq()[0]];
}


std::string Edge::getShortId() const {
    return start->getShortId() + "ACGT"[getSeq()[0]];
}

Sequence Edge::firstNucl() const {
    return getSeq().Subseq(0, 1);
}

//std::ostream &operator<<(std::ostream &os, const Edge &edge) {
//    os << edge.getShortId();
//    return os;
//}

Vertex::Vertex(hashing::htype hash, Vertex *_rc) : hash_(hash), rc_(_rc), canonical(false) {
    omp_init_lock(&writelock);
}

bool Vertex::isCanonical() const {
    return canonical;
}

size_t Vertex::coverage() const {
    return coverage_;
}

bool Vertex::isCanonical(const Edge &edge) const {
    const Vertex &other = edge.getFinish()->rc();
    if(hash() != other.hash())
        return hash() < other.hash();
    if (isCanonical() != other.isCanonical())
        return isCanonical();
    const Edge &rc_edge = edge.rc();
    return edge.getSeq() <= rc_edge.getSeq();
}

void Vertex::checkConsistency() const {
    for (const Edge &edge : outgoing_) {
        if (edge.getFinish() != nullptr) {
            if (edge.rc().getFinish() != &(this->rc())) {
                std::cout << this << " " << getSeq() << " " << edge.getSeq() << " " << edge.rc().getFinish() << " "
                          << &(this->rc()) << std::endl;
            }
            VERIFY(edge.rc().getFinish() == &(this->rc()));
            VERIFY(edge.start == this);
        }
    }
}

std::string Vertex::getId() const {
    std::stringstream ss;
    if(!isCanonical())
        ss << "-";
    ss << hash();
    return ss.str();
}

std::string Vertex::oldId() const {
    std::stringstream ss;
    ss << hash() << isCanonical();
    return ss.str();
}

std::string Vertex::getShortId() const {
    std::stringstream ss;
    if(!isCanonical())
        ss << "-";
    ss << hash() % 1000000000;
    return ss.str();
}

void Vertex::incCoverage() {
#pragma omp atomic update
    coverage_ += 1;
#pragma omp atomic update
    rc().coverage_ += 1;
}

void Vertex::setSeq(Sequence _seq) {
    lock();
    if (getSeq().empty()) {
        seq = std::move(_seq);
        unlock();
        rc_->lock();
        rc_->seq = !getSeq();
        rc_->unlock();
    } else {
        unlock();
    }
}

//void Vertex::clearSequence() {
//    if (!seq.empty()) {
//        seq = Sequence();
//        rc_->seq = Sequence();
//    }
//}

Edge &Vertex::addEdgeLockFree(const Edge &edge) {
    VERIFY(this == edge.getStart());
    for (Edge &e : outgoing_) {
        if (edge.size() <= e.size()) {
            if (edge.getSeq() == e.getSeq().Subseq(0, edge.size())) {
                return e;
            }
        } else if (edge.getSeq().Subseq(0, e.size()) == e.getSeq()) {
            e = edge;
            return e;
        }
    }
    outgoing_.emplace_back(edge);
    return outgoing_.back();
}

void Vertex::addEdge(const Edge &e) {
    omp_set_lock(&writelock);
    addEdgeLockFree(e);
    omp_unset_lock(&writelock);
}

Edge &Vertex::getOutgoing(unsigned char c) const {
    for (Edge &edge : outgoing_) {
        if (edge.getSeq()[0] == c) {
            return edge;
        }
    }
    std::cout << getSeq() << std::endl;
    std::cout << size_t(c) << std::endl;
    for (const Edge &edge : outgoing_) {
        std::cout << edge.getSeq() << std::endl;
    }
    VERIFY(false);
    return outgoing_.front();
}

bool Vertex::hasOutgoing(unsigned char c) const {
    for (const Edge &edge : outgoing_) {
        if (edge.getSeq()[0] == c) {
            return true;
        }
    }
    return false;
}

bool Vertex::operator<(const Vertex &other) const {
    return hash_ < other.hash_ || (hash_ == other.hash_ && canonical && !other.canonical);
}

bool Vertex::operator>(const Vertex &other) const {
    return hash_ > other.hash_ || (hash_ == other.hash_ && !canonical && other.canonical);
}

void Vertex::clear() {
    outgoing_.clear();
    rc_->outgoing_.clear();
}

void Vertex::clearOutgoing() {
    outgoing_.clear();
}

Vertex::Vertex(hashing::htype hash) : hash_(hash), rc_(new Vertex(hash, this)), canonical(true) {
    omp_init_lock(&writelock);
}

Vertex::~Vertex() {
    if (rc_ != nullptr) {
        rc_->rc_ = nullptr;
        delete rc_;
    }
    rc_ = nullptr;
}

void Vertex::sortOutgoing() {
    outgoing_.sort();
//    std::sort(outgoing_.begin(), outgoing_.end());
}

bool Vertex::isJunction() const {
    return outDeg() != 1 || inDeg() != 1;
}

bool Vertex::operator==(const Vertex &other) const {
    return this == &other;
}

bool Vertex::operator!=(const Vertex &other) const {
    return this != &other;
}

void SparseDBG::checkSeqFilled(size_t threads, logging::Logger &logger) {
    logger.trace() << "Checking vertex sequences" << std::endl;
    std::function<void(size_t, std::pair<const hashing::htype, Vertex> &)> task =
            [&logger](size_t pos, std::pair<const hashing::htype, Vertex> &pair) {
                const Vertex &vert = pair.second;
                if (vert.getSeq().empty() || vert.rc().getSeq().empty()) {
                    logger.trace() << "Sequence not filled " << pair.first << std::endl;
                    VERIFY(false);
                }
                if (!vert.isCanonical()) {
                    logger.trace() << "Canonical vertex marked not canonical " << pair.first << std::endl;
                    VERIFY(false);
                }
                if (vert.rc().isCanonical()) {
                    logger.trace() << "Noncanonical vertex marked canonical " << pair.first << std::endl;
                    VERIFY(false);
                }
            };
    processObjects(v.begin(), v.end(), logger, threads, task);
    logger.trace() << "Vertex sequence check success" << std::endl;
}

SparseDBG SparseDBG::Subgraph(std::vector<Segment<Edge>> &pieces) {
    SparseDBG res(hasher_);
    for(auto &it : v) {
        res.addVertex(it.second.getSeq());
    }
    for(Segment<Edge> &seg : pieces) {
        Vertex *left;
        Vertex *right;
        if (seg.left == 0) {
            left = &res.getVertex(seg.contig().getStart()->getSeq());
        } else {
            left = &res.addVertex(seg.contig().kmerSeq(seg.left));
        }
        Segment<Edge> rcSeg = seg.RC();
        if (rcSeg.left == 0) {
            right = &res.getVertex(rcSeg.contig().getStart()->getSeq());
        } else if(seg == rcSeg) {
            right = left;
        } else {
            right = &res.addVertex(rcSeg.contig().kmerSeq(rcSeg.left));
        }
        left->addEdge(Edge(*left, right->rc(), seg.seq()));
        right->addEdge(Edge(*right, left->rc(), rcSeg.seq()));
    }
    return std::move(res);
}

SparseDBG SparseDBG::SplitGraph(const std::vector<EdgePosition> &breaks) {
    SparseDBG res(hasher_);
    for(Vertex &it : verticesUnique()) {
        res.addVertex(it);
    }
    std::unordered_set<Edge *> broken_edges;
    for(const EdgePosition &epos : breaks) {
        if(!epos.isBorder()) {
            res.addVertex(epos.kmerSeq());
            broken_edges.emplace(epos.edge);
            broken_edges.emplace(&epos.edge->rc());
        }
    }
    for(Edge &edge : edges()) {
        if(broken_edges.find(&edge) == broken_edges.end()) {
            Vertex &start = res.getVertex(*edge.getStart());
            Vertex &end = res.getVertex(*edge.getFinish());
            Edge new_edge(start, end, edge.getSeq());
            start.addEdge(new_edge);
        } else {
            Vertex &newVertex = res.getVertex(*edge.getStart());
            res.processEdge(newVertex, edge.getSeq());
        }
    }
    return std::move(res);
}

SparseDBG SparseDBG::AddNewSequences(logging::Logger &logger, size_t threads, const std::vector<Sequence> &new_seqs) {
    SparseDBG res(hasher_);
    for(Vertex &it : verticesUnique()) {
        res.addVertex(it);
    }
    for(const Sequence &seq: new_seqs) {
        hashing::KWH kwh(hasher_, seq, 0);
        while(true) {
            if(!res.containsVertex(kwh.hash())) {
                res.addVertex(kwh);
            }
            if(!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
    }
    std::function<void(size_t, Edge &)> task = [&res](size_t num, Edge &edge) {
        res.processEdge(edge);
    };
    omp_set_num_threads(threads);
    processObjects(this->edgesUnique().begin(), this->edgesUnique().end(), logger, threads, task);
    std::function<void(size_t, const Sequence &)> task1 = [&res](size_t num, const Sequence &seq) {
        res.processRead(seq);
    };
    ParallelProcessor<const Sequence>(task1, logger, threads).processObjects(new_seqs.begin(), new_seqs.end(), 1024);
//    processObjects<std::vector<Sequence>::const_iterator>(new_seqs.begin(), new_seqs.end(), logger, threads, task1);
    return std::move(res);
}

void SparseDBG::checkConsistency(size_t threads, logging::Logger &logger) {
    logger.trace() << "Checking consistency" << std::endl;
    std::function<void(size_t, std::pair<const hashing::htype, Vertex> &)> task =
            [this](size_t pos, std::pair<const hashing::htype, Vertex> &pair) {
                const Vertex &vert = pair.second;
                vert.checkConsistency();
                vert.rc().checkConsistency();
            };
    processObjects(v.begin(), v.end(), logger, threads, task);
    logger.trace() << "Consistency check success" << std::endl;
}

void SparseDBG::checkDBGConsistency(size_t threads, logging::Logger &logger) {
    logger.trace() << "Checking kmer index" << std::endl;
    std::function<void(size_t, Edge &)> task =
            [this](size_t pos, Edge &edge) {
                hashing::KWH kwh(hasher(), edge.getStart()->getSeq() + edge.getSeq(), 0);
                while (true) {
                    if(this->containsVertex(kwh.hash())) {
                        VERIFY_OMP((kwh.pos == 0 && this->getVertex(kwh) == *edge.getStart()) || (kwh.pos == edge.size() && this->getVertex(kwh) == *edge.getFinish()), "Vertex kmer index corruption");
                    }
                    if(this->isAnchor(kwh.hash())) {
                        EdgePosition ep = getAnchor(kwh);
                        VERIFY_OMP(ep.edge == &edge && ep.pos == kwh.pos, "Anchor kmer index corruption " + itos(ep.pos) + " " +
                                itos(ep.edge->size()));
                    }
                    if(!kwh.hasNext())
                        break;
                    kwh = kwh.next();
                }
            };
    processObjects(edges().begin(), edges().end(), logger, threads, task);
    logger.trace() << "Index check success" << std::endl;
    size_t sz = 0;
    for(Edge &edge : edgesUnique()) {
        sz += edge.size();
    }
    for(Vertex &vertex: vertices()) {
        std::vector<char> out;
        for(Edge &edge : vertex) {
            out.emplace_back(edge.getSeq()[0]);
        }
        std::sort(out.begin(), out.end());
        for(size_t i = 0; i + 1 < out.size(); i++) {
            VERIFY(out[i] != out[i + 1]);
        }
    }
    if(sz > 10000000)
        return;
    ParallelRecordCollector<hashing::htype> hashs(threads);
    std::function<void(size_t, Edge &)> task1 =
            [this, &hashs](size_t pos, Edge &edge) {
                hashing::KWH kwh(hasher(), edge.getStart()->getSeq() + edge.getSeq(), 1);
                for(size_t i = 1; i < edge.size(); i++) {
                    hashs.emplace_back(kwh.hash());
                    kwh = kwh.next();
                }
            };
    processObjects(edgesUnique().begin(), edgesUnique().end(), logger, threads, task1);
    for(Vertex &vertex : verticesUnique()) {
        hashs.emplace_back(vertex.hash());
    }
    std::vector<hashing::htype> res = hashs.collect();
    __gnu_parallel::sort(res.begin(), res.end());
    bool ok = true;
    for(size_t i = 0; i + 1 < res.size(); i++) {
        ok &= res[i] != res[i + 1];
    }
    if(!ok) {
        logger.trace() << "Duplicated k-mers in the graph" << std::endl;
    }
}

Vertex &SparseDBG::addVertex(const hashing::KWH &kwh) {
    Vertex &newVertex = innerAddVertex(kwh.hash());
    Vertex &res = kwh.isCanonical() ? newVertex : newVertex.rc();
    res.setSeq(kwh.getSeq());
    return res;
}

Vertex &SparseDBG::addVertex(const Sequence &seq) {
    return addVertex(hashing::KWH(hasher_, seq, 0));
}

Vertex &SparseDBG::addVertex(const Vertex &other_graph_vertex) {
    Vertex &newVertex = innerAddVertex(other_graph_vertex.hash());
    if(other_graph_vertex.isCanonical()) {
        newVertex.setSeq(other_graph_vertex.getSeq());
        return newVertex;
    } else {
        newVertex.setSeq(!other_graph_vertex.getSeq());
        return newVertex.rc();
    }
}

Vertex &SparseDBG::bindTip(Vertex &start, Edge &tip) {
    Sequence seq = start.getSeq() + tip.getSeq();
    Vertex &end = addVertex(seq.Subseq(seq.size() - hasher().getK()));
    tip.bindTip(start, end);
    return end;
}

Vertex &SparseDBG::getVertex(const hashing::KWH &kwh) {
    auto it = v.find(kwh.hash());
    VERIFY(it != v.end());
    if (kwh.isCanonical()) {
        return it->second;
    } else {
        return it->second.rc();
    }
}

Vertex &SparseDBG::getVertex(const Sequence &seq) {
    return getVertex(hashing::KWH(hasher_, seq, 0));
}

Vertex &SparseDBG::getVertex(const Vertex &other_graph_vertex) {
    if(other_graph_vertex.isCanonical())
        return v.find(other_graph_vertex.hash())->second;
    else
        return v.find(other_graph_vertex.hash())->second.rc();
}

std::array<Vertex *, 2> SparseDBG::getVertices(hashing::htype hash) {
    Vertex &res = v.find(hash)->second;
    return {&res, &res.rc()};
}

void SparseDBG::fillAnchors(size_t w, logging::Logger &logger, size_t threads) {
    logger.trace() << "Adding anchors from long edges for alignment" << std::endl;
    ParallelRecordCollector<std::pair<const hashing::htype, EdgePosition>> res(threads);
    std::function<void(size_t, Edge &)> task = [&res, w, this](size_t pos, Edge &edge) {
        Vertex &vertex = *edge.getStart();
        if (edge.size() > w) {
            Sequence seq = vertex.getSeq() + edge.getSeq();
//                    Does not run for the first and last kmers.
            for (hashing::KWH kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                if (kmer.pos % w == 0) {
                    EdgePosition ep(edge, kmer.pos);
                    if (kmer.isCanonical())
                        res.emplace_back(kmer.hash(), ep);
                    else {
                        res.emplace_back(kmer.hash(), ep.RC());
                    }
                }
            }
        }
    };
    processObjects(edges().begin(), edges().end(), logger, threads, task);
    for (auto &tmp : res) {
        anchors.emplace(tmp);
    }
    anchors_filled = true;
    logger.trace() << "Added " << anchors.size() << " anchors" << std::endl;
}

void SparseDBG::fillAnchors(size_t w, logging::Logger &logger, size_t threads,
                            const std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> &to_add) {
    logger.trace() << "Adding anchors from long edges for alignment" << std::endl;
    ParallelRecordCollector<std::pair<const hashing::htype, EdgePosition>> res(threads);
    std::function<void(size_t, Edge &)> task = [&res, w, this, &to_add](size_t pos, Edge &edge) {
        Vertex &vertex = *edge.getStart();
        if (edge.size() > w || !to_add.empty()) {
            Sequence seq = vertex.getSeq() + edge.getSeq();
//                    Does not run for the first and last kmers.
            for (hashing::KWH kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                if (kmer.pos % w == 0 || to_add.find(kmer.hash()) != to_add.end()) {
                    EdgePosition ep(edge, kmer.pos);
                    if (kmer.isCanonical())
                        res.emplace_back(kmer.hash(), ep);
                    else {
                        res.emplace_back(kmer.hash(), ep.RC());
                    }
                }
            }
        }
    };
    processObjects(edges().begin(), edges().end(), logger, threads, task);
    for (auto &tmp : res) {
        anchors.emplace(tmp);
    }
    anchors_filled = true;
    logger.trace() << "Added " << anchors.size() << " anchors" << std::endl;
}

EdgePosition SparseDBG::getAnchor(const hashing::KWH &kwh) {
    if (kwh.isCanonical())
        return anchors.find(kwh.hash())->second;
    else
        return anchors.find(kwh.hash())->second.RC();
}

std::vector<hashing::KWH> SparseDBG::extractVertexPositions(const Sequence &seq, size_t max) const {
    std::vector<hashing::KWH> res;
    hashing::KWH kwh(hasher(), seq, 0);
    while (true) {
        if (containsVertex(kwh.hash())) {
            res.emplace_back(kwh);
        }
        if (!kwh.hasNext() || res.size() == max)
            break;
        kwh = kwh.next();
    }
    return std::move(res);
}

void SparseDBG::printFastaOld(const std::experimental::filesystem::path &out) {
    std::ofstream os;
    os.open(out);
    for(Edge &edge : edges()) {
        os << ">" << edge.getStart()->hash() << edge.getStart()->isCanonical() << "ACGT"[edge.getSeq()[0]] << "\n" <<
           edge.getStart()->getSeq() << edge.getSeq() << "\n";
    }
    os.close();
}

void SparseDBG::processRead(const Sequence &seq) {
    std::vector<hashing::KWH> kmers = extractVertexPositions(seq);
    if (kmers.size() == 0) {
        std::cout << seq << std::endl;
    }
    VERIFY(kmers.size() > 0);
    std::vector<Vertex *> vertices;
    for (size_t i = 0; i < kmers.size(); i++) {
        vertices.emplace_back(&getVertex(kmers[i]));
        if (i == 0 || vertices[i] != vertices[i - 1]) {
            vertices.back()->setSeq(kmers[i].getSeq().copy());
            vertices.back()->incCoverage();
        }
    }
    for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
        VERIFY(kmers[i].pos + hasher_.getK() <= seq.size())
        if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
            (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
            kmers[i + 1].pos - kmers[i].pos < hasher_.getK()) {
            continue;
        }
        vertices[i]->addEdge(Edge(*vertices[i], *vertices[i + 1], Sequence(seq.Subseq(kmers[i].pos + hasher_.getK(),
                                                                                    kmers[i + 1].pos +
                                                                                    hasher_.getK()).str())));
        vertices[i + 1]->rc().addEdge(Edge(vertices[i + 1]->rc(), vertices[i]->rc(),
                                           !Sequence(seq.Subseq(kmers[i].pos, kmers[i + 1].pos).str())));
    }
    if (kmers.front().pos > 0) {
        vertices.front()->rc().addSequence(!(seq.Subseq(0, kmers[0].pos)));
    }
    if (kmers.back().pos + hasher_.getK() < seq.size()) {
        vertices.back()->addSequence(seq.Subseq(kmers.back().pos + hasher_.getK(), seq.size()));
    }
}

//This method does not add rc edges so should be run for both edge and its rc
void SparseDBG::processEdge(Vertex &vertex, Sequence old_seq) {
    Sequence seq = vertex.getSeq() + old_seq;
    std::vector<hashing::KWH> kmers = extractVertexPositions(seq);
    VERIFY(kmers.front().pos == 0 && kmers.back().pos == old_seq.size());
    std::vector<Vertex *> vertices;
    for (size_t i = 0; i < kmers.size(); i++) {
        vertices.emplace_back(&getVertex(kmers[i]));
    }
    for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
        VERIFY(kmers[i].pos + hasher_.getK() <= seq.size())
        if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
            (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
            kmers[i + 1].pos - kmers[i].pos < hasher_.getK()) {
            continue;
        }
        vertices[i]->addEdge(
                Edge(*vertices[i], *vertices[i + 1], old_seq.Subseq(kmers[i].pos, kmers[i + 1].pos)));
    }
}

void SparseDBG::processEdge(Edge &other_graph_edge) {
    Sequence seq = other_graph_edge.getStart()->getSeq() + other_graph_edge.getSeq();
    std::vector<hashing::KWH> kmers = extractVertexPositions(seq);
    VERIFY(kmers.front().pos == 0 && kmers.back().pos == other_graph_edge.size());
    std::vector<Vertex *> vertices;
    for (size_t i = 0; i < kmers.size(); i++) {
        vertices.emplace_back(&getVertex(kmers[i]));
    }
    for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
        VERIFY(kmers[i].pos + hasher_.getK() <= seq.size())
        if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
            (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
            kmers[i + 1].pos - kmers[i].pos < hasher_.getK()) {
            continue;
        }
        vertices[i]->addEdge(Edge(*vertices[i], *vertices[i + 1], other_graph_edge.getSeq().Subseq(kmers[i].pos, kmers[i + 1].pos)));
        vertices[i + 1]->rc().addEdge(Edge(vertices[i + 1]->rc(), vertices[i]->rc(),
                  other_graph_edge.rc().getSeq().Subseq(other_graph_edge.size() - kmers[i + 1].pos, other_graph_edge.size() - kmers[i].pos)));
    }
}

IterableStorage<ApplyingIterator<SparseDBG::vertex_iterator_type, Vertex, 2>> SparseDBG::vertices(bool unique) {
    std::function<std::array<Vertex*, 2>(std::pair<const hashing::htype, Vertex> &)> apply =
            [unique](std::pair<const hashing::htype, Vertex> &it) -> std::array<Vertex*, 2> {
                if(unique)
                    return {&it.second};
                else
                    return {&it.second, &it.second.rc()};
            };
    ApplyingIterator<vertex_iterator_type, Vertex, 2> begin(v.begin(), v.end(), apply);
    ApplyingIterator<vertex_iterator_type, Vertex, 2> end(v.end(), v.end(), apply);
    return {begin, end};
}

IterableStorage<ApplyingIterator<SparseDBG::vertex_iterator_type, Vertex, 2>> SparseDBG::verticesUnique() {
    return vertices(true);
}

IterableStorage<ApplyingIterator<SparseDBG::vertex_iterator_type, Edge, 8>> SparseDBG::edges(bool unique) {
    std::function<std::array<Edge*, 8>(const std::pair<const hashing::htype, Vertex> &)> apply = [unique](const std::pair<const hashing::htype, Vertex> &it) {
        std::array<Edge*, 8> res = {};
        size_t cur = 0;
        for(Edge &edge : it.second) {
            if(!unique || edge <= edge.rc()) {
                res[cur] = &edge;
                cur++;
            }
        }
        for(Edge &edge : it.second.rc()) {
            if(!unique || edge <= edge.rc()) {
                res[cur] = &edge;
                cur++;
            }
        }
        return res;
    };
    ApplyingIterator<vertex_iterator_type, Edge, 8> begin(v.begin(), v.end(), apply);
    ApplyingIterator<vertex_iterator_type, Edge, 8> end(v.end(), v.end(), apply);
    return {begin, end};
}

IterableStorage<ApplyingIterator<SparseDBG::vertex_iterator_type, Edge, 8>> SparseDBG::edgesUnique() {
    return edges(true);
}

void SparseDBG::removeIsolated() {
    vertex_map_type newv;
    std::vector<hashing::htype> todelete;
    for (auto it = v.begin(); it != v.end();) {
        if (it->second.outDeg() == 0 && it->second.inDeg() == 0) {
            it = v.erase(it);
        } else {
            ++it;
        }
    }
}

void SparseDBG::removeMarked() {
    for (auto it = v.begin(); it != v.end();) {
        if (it->second.marked() || (it->second.inDeg() == 0 && it->second.outDeg() == 0)) {
            it = v.erase(it);
        } else {
            ++it;
        }
    }
}

void SparseDBG::resetMarkers() {
    for(dbg::Edge &edge: edges()) {
        edge.mark(dbg::EdgeMarker::common);
    }
    for(dbg::Vertex &vertex : vertices()) {
        vertex.unmark();
    }
}

//const Vertex &SparseDBG::getVertex(const hashing::KWH &kwh) const {
//    auto it = v.find(kwh.hash());
//    VERIFY(it != v.end());
//    if (kwh.isCanonical()) {
//        return it->second;
//    } else {
//        return it->second.rc();
//    }
//}
std::vector<EdgePosition> EdgePosition::step() const {
    if (pos == edge->size()) {
        std::vector<EdgePosition> res;
        Vertex &v = *edge->getFinish();
        for (Edge &next : v) {
            res.emplace_back(next, 1);
        }
        return std::move(res);
    } else {
        return {{*edge, pos + 1}};
    }
}
