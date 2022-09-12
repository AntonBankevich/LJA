#include "sparse_dbg.hpp"
using namespace dbg;

Edge Edge::_fake = Edge(nullptr, nullptr, Sequence());

bool IsMarkerCorrect(EdgeMarker marker) {
    return marker == EdgeMarker::correct || marker == EdgeMarker::unique || marker == EdgeMarker::repeat;
}

size_t Edge::updateTipSize() const {
    size_t new_val = 0;
    if(extraInfo == size_t(-1) && end_->inDeg() == 1) {
        for (const Edge & other : *end_) {
            other.end_->lock();
            new_val = std::max(new_val, other.extraInfo);
            other.end_->unlock();
        }
        if(new_val != size_t(-1))
            new_val += size();
        end_->lock();
        extraInfo = new_val;
        end_->unlock();
    }
    return new_val;
}

Edge &Edge::rc() const {
    VERIFY(!start()->seq.empty());
    Vertex &vend = end_->rc();
    unsigned char c;
    size_t k = vend.seq.size();
    if(size() > k) {
        c = (!seq)[k];
    } else {
        c = (!start_->seq)[k - size()];
    }
    return vend.getOutgoing(c);
}

Edge &Edge::sparseRcEdge() const {
    Vertex &vend = end()->rc();
    VERIFY(!start_->seq.empty());
    for(Edge &candidate : vend) {
        if(*candidate.end() == start_->rc() && candidate.size() == size() &&
           (size() <= start_->seq.size() || candidate.seq.startsWith((!seq).Subseq(start_->seq.size())))) {
            return candidate;
        }
    }
    std::cout << start_->seq + seq << std::endl;
    std::cout << seq << std::endl;
    std::cout << vend.seq << std::endl;
    for(Edge &candidate : vend) {
        if(*candidate.end() == start_->rc() && candidate.size() == size() &&
           (size() <= seq.size() || candidate.seq.startsWith((!seq).Subseq(start_->seq.size())))) {
            std::cout << vend.seq + candidate.seq << std::endl;
        }
    }
    VERIFY(false);
    return vend[0];
}

void Edge::bindTip(Vertex &start, Vertex &end) {
    VERIFY(end_ == nullptr);
    end_ = &end;
    Sequence rcseq = !(start.seq + seq);
    end.rc().addEdgeLockFree(Edge(&end.rc(), &start.rc(), rcseq.Subseq(start.seq.size())));
}

size_t Edge::getTipSize() const {
    return extraInfo;
}

void Edge::setTipSize(size_t val) const {
    extraInfo = val;
}

Vertex *Edge::start() const {
    return start_;
}

Vertex *Edge::end() const {
    return end_;
}

size_t Edge::size() const {
    return seq.size();
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
    if(start_ != other.start_)
        return *start_ < *other.start_;
    return this->seq < other.seq;
}

bool Edge::operator>(const Edge &other) const {
    if(this == &other)
        return false;
    if(start_ != other.start_)
        return *start_ > *other.start_;
    return other.seq < seq;
}

bool Edge::operator<=(const Edge &other) const {
    return *this == other || *this < other;
}

std::string Edge::str() const {
    std::stringstream ss;
    const dbg::Vertex &v = *start();
    ss << v.hash() << v.isCanonical() << "ACGT"[seq[0]];
    return ss.str();
}

Sequence Edge::kmerSeq(size_t pos) const {
    VERIFY(pos <= seq.size());
    size_t k = start_->seq.size();
    if (pos >= k)
        return seq.Subseq(pos - k, pos);
    else {
        return start()->seq.Subseq(pos) + seq.Subseq(0, pos);
    }
}

Sequence Edge::suffix(size_t pos) const {
    VERIFY(pos <= seq.size());
    size_t k = start_->seq.size();
    if (pos >= k)
        return seq.Subseq(pos - k, seq.size());
    else {
        return start()->seq.Subseq(pos) + seq.Subseq(0, seq.size());
    }
}

std::string Edge::getId() const {
    return start_->getId() + "ACGT"[seq[0]];
}

std::string Edge::oldId() const {
    return start_->oldId() + "ACGT"[seq[0]];
}


std::string Edge::getShortId() const {
    return start_->getShortId() + "ACGT"[seq[0]];
}

Sequence Edge::firstNucl() const {
    return seq.Subseq(0, 1);
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
    const Vertex &other = edge.end()->rc();
    if(hash() != other.hash())
        return hash() < other.hash();
    if (isCanonical() != other.isCanonical())
        return isCanonical();
    const Edge &rc_edge = edge.rc();
    return edge.seq <= rc_edge.seq;
}

void Vertex::checkConsistency() const {
    for (const Edge &edge : outgoing_) {
        if (edge.end() != nullptr) {
            if (edge.rc().end() != &(this->rc())) {
                std::cout << this << " " << seq << " " << edge.seq << " " << edge.rc().end() << " "
                          << &(this->rc()) << std::endl;
            }
            VERIFY(edge.rc().end() == &(this->rc()));
            VERIFY(edge.start_ == this);
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

void Vertex::setSequence(const Sequence &_seq) {
    lock();
    if (seq.empty()) {
        if (seq.empty()) {
            seq = Sequence(_seq.str());
            unlock();
            rc_->lock();
            rc_->seq = !_seq;
            rc_->unlock();
        } else {
            unlock();
        }
    } else {
//            if(seq != _seq) {
//                std::cout << seq << std::endl << _seq << std::endl;
//                VERIFY(false);
//            }
//            VERIFY(_seq == seq);
        unlock();
    }
}

void Vertex::clearSequence() {
    if (!seq.empty()) {
        seq = Sequence();
        rc_->seq = Sequence();
    }
}

Edge &Vertex::addEdgeLockFree(const Edge &edge) {
    VERIFY(this == edge.start());
    for (Edge &e : outgoing_) {
        if (edge.size() <= e.size()) {
            if (edge.seq == e.seq.Subseq(0, edge.size())) {
                return e;
            }
        } else if (edge.seq.Subseq(0, e.size()) == e.seq) {
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
        if (edge.seq[0] == c) {
            return edge;
        }
    }
    std::cout << seq << std::endl;
    std::cout << size_t(c) << std::endl;
    for (const Edge &edge : outgoing_) {
        std::cout << edge.seq << std::endl;
    }
    VERIFY(false);
    return outgoing_[0];
}

bool Vertex::hasOutgoing(unsigned char c) const {
    for (const Edge &edge : outgoing_) {
        if (edge.seq[0] == c) {
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
    std::sort(outgoing_.begin(), outgoing_.end());
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
                if (vert.seq.empty() || vert.rc().seq.empty()) {
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
        res.addVertex(it.second.seq);
    }
    for(Segment<Edge> &seg : pieces) {
        Vertex *left;
        Vertex *right;
        if (seg.left == 0) {
            left = &res.getVertex(seg.contig().start()->seq);
        } else {
            left = &res.addVertex(seg.contig().kmerSeq(seg.left));
        }
        Segment<Edge> rcSeg = seg.RC();
        if (rcSeg.left == 0) {
            right = &res.getVertex(rcSeg.contig().start()->seq);
        } else if(seg == rcSeg) {
            right = left;
        } else {
            right = &res.addVertex(rcSeg.contig().kmerSeq(rcSeg.left));
        }
        left->addEdge(Edge(left, &right->rc(), seg.seq()));
        right->addEdge(Edge(right, &left->rc(), rcSeg.seq()));
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
            Vertex &start = res.getVertex(*edge.start());
            Vertex &end = res.getVertex(*edge.end());
            Edge new_edge(&start, &end, edge.seq);
            start.addEdge(new_edge);
        } else {
            Vertex &newVertex = res.getVertex(*edge.start());
            res.processEdge(newVertex, edge.seq);
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
                hashing::KWH kwh(hasher(), edge.start()->seq + edge.seq, 0);
                while (true) {
                    if(this->containsVertex(kwh.hash())) {
                        VERIFY_OMP((kwh.pos == 0 && this->getVertex(kwh) == *edge.start()) || (kwh.pos == edge.size() && this->getVertex(kwh) == *edge.end()), "Vertex kmer index corruption");
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
            out.emplace_back(edge.seq[0]);
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
                hashing::KWH kwh(hasher(), edge.start()->seq + edge.seq, 1);
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
    res.setSequence(kwh.getSeq());
    return res;
}

Vertex &SparseDBG::addVertex(const Sequence &seq) {
    return addVertex(hashing::KWH(hasher_, seq, 0));
}

Vertex &SparseDBG::addVertex(const Vertex &other_graph_vertex) {
    Vertex &newVertex = innerAddVertex(other_graph_vertex.hash());
    if(other_graph_vertex.isCanonical()) {
        newVertex.setSequence(other_graph_vertex.seq);
        return newVertex;
    } else {
        newVertex.setSequence(!other_graph_vertex.seq);
        return newVertex.rc();
    }
}

Vertex &SparseDBG::bindTip(Vertex &start, Edge &tip) {
    Sequence seq = start.seq + tip.seq;
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
        Vertex &vertex = *edge.start();
        if (edge.size() > w) {
            Sequence seq = vertex.seq + edge.seq;
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
    logger.trace() << "Added " << anchors.size() << " anchors" << std::endl;
}

void SparseDBG::fillAnchors(size_t w, logging::Logger &logger, size_t threads,
                            const std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> &to_add) {
    logger.trace() << "Adding anchors from long edges for alignment" << std::endl;
    ParallelRecordCollector<std::pair<const hashing::htype, EdgePosition>> res(threads);
    std::function<void(size_t, Edge &)> task = [&res, w, this, &to_add](size_t pos, Edge &edge) {
        Vertex &vertex = *edge.start();
        if (edge.size() > w || !to_add.empty()) {
            Sequence seq = vertex.seq + edge.seq;
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
        os << ">" << edge.start()->hash() << edge.start()->isCanonical() << "ACGT"[edge.seq[0]] << "\n" <<
           edge.start()->seq << edge.seq << "\n";
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
            vertices.back()->setSequence(kmers[i].getSeq());
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
        vertices[i]->addEdge(Edge(vertices[i], vertices[i + 1], Sequence(seq.Subseq(kmers[i].pos + hasher_.getK(),
                                                                                    kmers[i + 1].pos +
                                                                                    hasher_.getK()).str())));
        vertices[i + 1]->rc().addEdge(Edge(&vertices[i + 1]->rc(), &vertices[i]->rc(),
                                           !Sequence(seq.Subseq(kmers[i].pos, kmers[i + 1].pos).str())));
    }
    if (kmers.front().pos > 0) {
        vertices.front()->rc().addEdge(Edge(&vertices.front()->rc(), nullptr, !(seq.Subseq(0, kmers[0].pos))));
    }
    if (kmers.back().pos + hasher_.getK() < seq.size()) {
        vertices.back()->addEdge(
                Edge(vertices.back(), nullptr, seq.Subseq(kmers.back().pos + hasher_.getK(), seq.size())));
    }
}

//This method does not add rc edges so should be run for both edge and its rc
void SparseDBG::processEdge(Vertex &vertex, Sequence old_seq) {
    Sequence seq = vertex.seq + old_seq;
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
                Edge(vertices[i], vertices[i + 1], old_seq.Subseq(kmers[i].pos, kmers[i + 1].pos)));
    }
}

void SparseDBG::processEdge(Edge &other_graph_edge) {
    Sequence seq = other_graph_edge.start()->seq + other_graph_edge.seq;
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
        vertices[i]->addEdge(Edge(vertices[i], vertices[i + 1], other_graph_edge.seq.Subseq(kmers[i].pos, kmers[i + 1].pos)));
        vertices[i + 1]->rc().addEdge(Edge(&vertices[i + 1]->rc(), &vertices[i]->rc(),
                  other_graph_edge.rc().seq.Subseq(other_graph_edge.size() - kmers[i + 1].pos, other_graph_edge.size() - kmers[i].pos)));
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
    for(dbg::Edge &edge: edges())
        edge.mark(dbg::EdgeMarker::common);
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
        Vertex &v = *edge->end();
        for (Edge &next : v) {
            res.emplace_back(next, 1);
        }
        return std::move(res);
    } else {
        return {{*edge, pos + 1}};
    }
}
