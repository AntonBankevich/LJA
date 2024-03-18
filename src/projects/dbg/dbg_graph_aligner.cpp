#include "dbg_graph_aligner.hpp"

using namespace dbg;

PerfectAlignment<Contig, dbg::Edge> bestExtension(const Vertex &vertex, const Segment<Contig> &seg) {
    PerfectAlignment<Contig, dbg::Edge> best({seg.contig(), seg.left, seg.left}, {});
    for (Edge &edge: vertex) {
        size_t len = 0;
        while (len < edge.truncSize() && seg.left + vertex.size() + len < seg.contig().truncSize()) {
            if (seg.contig()[seg.left + vertex.size() + len] != edge.truncSeq()[len])
                break;
            len++;
        }
//        std::cout << len << std::endl;
//        std::cout << Segment<Contig>(seg.contig(), seg.left + vertex.seq.size(),
//                                     std::min(seg.contig().size(), seg.left + vertex.seq.size() + 200)).seq() << std::endl;
//        std::cout << edge.seq.Subseq(0, std::min<size_t>(edge.seq.size(), 200)) << std::endl;
        if (len >= best.size()) {
            best = {{seg.contig(), seg.left, seg.left + len},
                    {edge,         0,        len}};
        }
    }
    return best;
}

PerfectAlignment<Contig, dbg::Edge> bestExtension(Edge &edge, const Segment<Contig> &seg) {
    size_t len = 0;
    while (len < edge.truncSize() && seg.left + edge.getStart().getSeq().size() + len < seg.contig().truncSize()) {
        if (seg.contig()[seg.left + edge.getStart().getSeq().size() + len] != edge.truncSeq()[len])
            break;
        len++;
    }
    return {Segment<Contig>(seg.contig(), seg.left, seg.left + len), Segment<Edge>(edge, 0, len)};
}

dbg::GraphPath KmerIndex::align(const Sequence &seq, dbg::Edge *edge_to, size_t pos_to) {
    VERIFY(alignmentReady());
    size_t k = hasher().getK();
    size_t cur = k;
    dbg::GraphPath res;
    while(cur < seq.size()) {
        size_t len = std::min(seq.size() - cur, edge_to->truncSize() - pos_to);
        res += Segment<dbg::Edge>(*edge_to, pos_to, pos_to + len);
        cur += len;
        if(cur < seq.size()) {
            edge_to = &edge_to->getFinish().getOutgoing(seq[cur]);
            pos_to = 0;
        }
    }
    return res;
}

dbg::GraphPath KmerIndex::align(const Sequence &seq, const std::string &name) const {
    VERIFY(alignmentReady());
    std::vector<hashing::MovingKWH> kmers = extractVertexPositions(seq, 1);
    size_t k = hasher().getK();
    dbg::GraphPath res;
    if (kmers.empty()) {
        for(const hashing::MovingKWH &kwh : hasher().kmers(seq)) {
            if (isAnchor(kwh.hash())) {
                dbg::EdgePosition pos = getAnchor(kwh);
                VERIFY(kwh.getPos() < pos.pos);
                VERIFY(pos.pos + seq.size() - kwh.getPos() <= pos.edge->fullSize());
                Segment<dbg::Edge> seg(*pos.edge, pos.pos - kwh.getPos(), pos.pos + seq.size() - kwh.getPos() - k);
                return {seg};
            }
        }
#pragma omp critical
        {
            std::cout << "Error: could not find alignment anchors in sequence " << name << " of size " << seq.size() << std::endl;
            std::cout << seq << std::endl;
            abort();
        };
        return {};
    }
    dbg::Vertex *prestart = &getVertex(kmers.front());
    if (kmers.front().getPos() > 0) {
        dbg::Vertex &rcstart = prestart->rc();
        if (!rcstart.hasOutgoing(seq[kmers.front().getPos() - 1] ^ 3)) {
            std::cout << "No incoming for start vertex" << std::endl << seq << std::endl <<
                      kmers.front().getPos() << " " << seq[kmers.front().getPos() - 1] << std::endl
                      << kmers.front().getSeq() << std::endl;
            VERIFY(false);
        }
        dbg::Edge &rcedge = rcstart.getOutgoing(seq[kmers.front().getPos() - 1] ^ 3);
        dbg::Edge &edge = rcedge.rc();
        VERIFY(edge.truncSize() >= kmers.front().getPos());
        Segment<dbg::Edge> seg(edge, edge.truncSize() - kmers.front().getPos(), edge.truncSize());
        res += seg;
    }
    size_t cpos = kmers.front().getPos() + k;
    while(cpos < seq.size()) {
        if(seq.Subseq(cpos -k, cpos) != prestart->getSeq() || !prestart->hasOutgoing(seq[cpos])) {
            std::cout << "No outgoing for middle\n" << seq << "\n" << cpos << " " << prestart->getInnerId() <<
                      " " << prestart->outDeg() << "\n" << seq.Subseq(cpos -k, cpos) << "\n" << prestart->getSeq() << "\n" <<
                      size_t(seq[cpos]) << std::endl;
            std::cout << (seq.Subseq(cpos -k, cpos) != prestart->getSeq()) << " " <<  !prestart->hasOutgoing(seq[cpos]) << std::endl;
            std::cout << hashing::MovingKWH(hasher(), seq.Subseq(cpos - k, cpos), 0).hash() << " " <<
                      hashing::MovingKWH(hasher(), prestart->getSeq(), 0).hash() << std::endl;
            for(dbg::Edge &tmp : *prestart) {
                std::cout << tmp.getInnerId() << " " << tmp.truncSize() << std::endl;
            }
            VERIFY(false);
        }
        dbg::Edge &next = prestart->getOutgoing(seq[cpos]);
        size_t len = std::min<size_t>(next.truncSize(), seq.size() - cpos);
        res += Segment<dbg::Edge>(next, 0, len);
        cpos += len;
        prestart = &next.getFinish();
    }
    return std::move(res);
}

dbg::GraphPath KmerIndex::align(const dbg::EdgePosition &pos, const Sequence &seq) const {
    VERIFY(alignmentReady());
    dbg::GraphPath res(Segment<dbg::Edge>(*pos.edge, pos.pos, pos.pos));
    dbg::Edge *cedge = pos.edge;
    size_t epos = pos.pos;
    for (size_t cpos = 0; cpos < seq.size(); cpos++) {
        unsigned char c = seq[cpos];
        if (epos == cedge->truncSize()) {
            dbg::Vertex &vertex = cedge->getFinish();
            if (vertex.hasOutgoing(c)) {
                cedge = &vertex.getOutgoing(c);
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

std::vector<dbg::PerfectAlignment<dbg::Edge, dbg::Edge>> KmerIndex::oldEdgeAlign(dbg::Edge &contig) const {
    VERIFY(alignmentReady());
    Sequence seq = contig.getSeq();
    std::vector<PerfectAlignment < dbg::Edge, dbg::Edge>> res;
    size_t k = hasher().getK();
    for(const hashing::MovingKWH &kwh : hasher().kmers(seq, 0, seq.size() - hasher().getK())) {
        if (res.empty() || kwh.getPos() >= res.back().seg_from.right) {
            dbg::Edge *edge = nullptr;
            size_t pos = 0;
            if (containsVertex(kwh.hash())) {
                dbg::Vertex &start = getVertex(kwh);
                if (start.hasOutgoing(seq[kwh.getPos() + k]))
                    edge = &getVertex(kwh).getOutgoing(seq[kwh.getPos() + k]);
            }
            if (edge == nullptr && isAnchor(kwh.hash())) {
                dbg::EdgePosition gpos = getAnchor(kwh);
                if (gpos.edge->truncSeq()[gpos.pos] == seq[kwh.getPos() + k]) {
                    edge = gpos.edge;
                    pos = gpos.pos;
                }
            }
            if (edge != nullptr) {
                size_t len = std::min(contig.truncSize() - kwh.getPos(), edge->truncSize() - pos);
                res.emplace_back(Segment<dbg::Edge>(contig, kwh.getPos(), kwh.getPos() + len), Segment<dbg::Edge>(*edge, pos, pos + len));
            }
        }
    }
    return std::move(res);
}

std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> dbg::KmerIndex::carefulAlign(Contig &contig) const {
    VERIFY(alignmentReady());
    Sequence seq = contig.getSeq();
    size_t k = hasher().getK();

    if(contig.truncSize() < k) {
        return {};
    }
    std::vector<PerfectAlignment<Contig, Edge>> res;
    for(const hashing::MovingKWH &kwh : hasher().kmers(seq)) {
        if (res.empty() || kwh.getPos() >= res.back().seg_from.right) {
            if (containsVertex(kwh.hash())) {
                Vertex &vertex = getVertex(kwh);
                Vertex &rcVertex = vertex.rc();
                if ((res.empty() || kwh.getPos() > res.back().seg_from.right)
                    && kwh.getPos() > 0 && rcVertex.hasOutgoing(seq[kwh.getPos() - 1] ^ 3)) {
                    Edge &edge = rcVertex.getOutgoing(seq[kwh.getPos() - 1] ^ 3);
                    size_t len = 1;
                    while (len < edge.truncSize() && len < kwh.getPos() && edge.truncSeq()[len] == (seq[kwh.getPos() - len - 1] ^ 3))
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.getPos() - len, kwh.getPos()),
                                     Segment<Edge>(edge.rc(), edge.truncSize() - len, edge.truncSize()));
                }
                if (kwh.getPos() + k < seq.size() && vertex.hasOutgoing(seq[kwh.getPos() + k])) {
                    Edge &edge = vertex.getOutgoing(seq[kwh.getPos() + k]);
                    size_t len = 1;
                    while (len < edge.truncSize() && kwh.getPos() + k + len < seq.size() &&
                           edge.truncSeq()[len] == seq[kwh.getPos() + k + len])
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.getPos(), kwh.getPos() + len),
                                     Segment<Edge>(edge, 0, len));
                }
            } else if ((res.empty() || kwh.getPos() > res.back().seg_from.right) && isAnchor(kwh.hash())) {
                EdgePosition pos = getAnchor(kwh);
//                TODO replace this code with a call to expand method of PerfectAlignment class after each edge is marked by its full sequence
                Edge &edge = *pos.edge;
                Vertex &start = pos.edge->getStart();
                CompositeSequence edge_seq({start.getSeq(), edge.truncSeq()});
                size_t left_from = kwh.getPos();
                size_t right_from = kwh.getPos() + k;
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
    }
    return std::move(res);
}

PerfectAlignment<Contig, dbg::Edge> KmerIndex::extendLeft(const hashing::MovingKWH &kwh, Contig &contig) const {
    size_t k = hasher().getK();
    PerfectAlignment<Contig, dbg::Edge> best({contig, kwh.getPos(), kwh.getPos()}, {});
    if(kwh.getPos() == 0) {
        return best;
    }
    Vertex &start = getVertex(kwh);
    Contig rc_contig = contig.RC();
    if(start.inDeg() == 0) {
        return best;
    }
    PerfectAlignment<Contig, Edge> start_al = bestExtension(start.rc(), Segment<Contig>(rc_contig, contig.truncSize() - k - kwh.getPos(), contig.truncSize() - k));
    return {Segment<Contig>(contig, contig.truncSize() - k - start_al.seg_from.right, contig.truncSize() - k - start_al.seg_from.left),
            Segment<Edge>(start_al.seg_to.contig().rc(),
                          start_al.seg_to.contig().truncSize() - start_al.seg_to.right,
                          start_al.seg_to.contig().truncSize() - start_al.seg_to.left)};
}

PerfectAlignment<Contig, dbg::Edge> KmerIndex::extendRight(const hashing::MovingKWH &kwh, Contig &contig) const {
    size_t k = hasher().getK();
    PerfectAlignment<Contig, dbg::Edge> best({contig, kwh.getPos(), kwh.getPos()}, {});
    if(kwh.getPos() + k == contig.truncSize()) {
        return best;
    }
    Vertex &start = getVertex(kwh);
    if(start.outDeg() == 0) {
        return best;
    }
    return bestExtension(start, Segment<Contig>(contig, kwh.getPos(), contig.truncSize() - k));
}

std::vector<PerfectAlignment<Contig, dbg::Edge>> KmerIndex::sparseAlign(Contig &contig) const {
    VERIFY(alignmentReady());
    std::vector<hashing::MovingKWH> vlist = extractVertexPositions(contig.getSeq());
    std::vector<PerfectAlignment<Contig, dbg::Edge>> result;
    if(vlist.empty())
        return result;
    for(hashing::MovingKWH &kwh : vlist) {
        if(result.empty() || result.back().seg_from.right != kwh.getPos()) {
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

Vertex &KmerIndex::getVertex(const hashing::KWH &kwh) const {
    return getVertex(kwh.hash(), kwh.getSeq(hasher_.getK()).isCanonical());
}

Vertex &KmerIndex::getVertex(const Sequence &seq) const {
    return getVertex(hashing::MovingKWH(hasher_, seq, 0));
}

Vertex &KmerIndex::getVertex(const Vertex &other_graph_vertex) const {
    return getVertex(other_graph_vertex.getHash(), other_graph_vertex.isCanonical());
}

void KmerIndex::fillAnchors(logging::Logger &logger, size_t threads, SparseDBG &dbg, size_t _w) {
    w = _w;
    logger.trace() << "Adding anchors from long edges for alignment" << std::endl;
    ParallelRecordCollector<std::pair<const hashing::htype, EdgePosition>> res(threads);
    std::function<void(size_t, Edge &)> task = [&res, this](size_t pos, Edge &edge) {
        Vertex &vertex = edge.getStart();
        if (edge.truncSize() > w) {
            Sequence seq = vertex.getSeq() + edge.truncSeq();
//                    Does not run for the first and last kmers.
            for(const hashing::MovingKWH &kmer : hasher().innerKmers(seq)) {
                if (kmer.getPos() % w == 0) {
                    EdgePosition ep(edge, kmer.getPos());
                    if (kmer.isCanonical())
                        res.emplace_back(kmer.hash(), ep);
                    else {
                        res.emplace_back(kmer.hash(), ep.RC());
                    }
                }
            }
        }
    };
    processObjects(dbg.edges().begin(), dbg.edges().end(), logger, threads, task);
    for (auto &tmp : res) {
        anchors.emplace(tmp);
    }
    anchors_filled = true;
    logger.trace() << "Added " << anchors.size() << " anchors" << std::endl;
}

void KmerIndex::fillAnchors(logging::Logger &logger, size_t threads, SparseDBG &dbg, size_t _w,
                            const std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> &to_add) {
    w = _w;
    logger.trace() << "Adding anchors from long edges for alignment" << std::endl;
    ParallelRecordCollector<std::pair<const hashing::htype, EdgePosition>> res(threads);
    std::function<void(size_t, Edge &)> task = [&res, this, &to_add](size_t pos, Edge &edge) {
        Vertex &vertex = edge.getStart();
        if (edge.truncSize() > w || !to_add.empty()) {
            Sequence seq = edge.getSeq();
//                    Does not run for the first and last kmers.
            for(const hashing::MovingKWH &kmer : hasher().innerKmers(seq)) {
                VERIFY_OMP(v.find(kmer.hash()) == v.end());
                if (kmer.getPos() % w == 0 || to_add.find(kmer.hash()) != to_add.end()) {
                    EdgePosition ep(edge, kmer.getPos());
                    VERIFY_MSG(!containsVertex(kmer.hash()), "Kmer is present both as vertex and as edge anchor");
                    if (kmer.isCanonical())
                        res.emplace_back(kmer.hash(), ep);
                    else {
                        res.emplace_back(kmer.hash(), ep.RC());
                    }
                }
            }
        }
    };
    processObjects(dbg.edges().begin(), dbg.edges().end(), logger, threads, task);
    for (auto &tmp : res) {
        anchors.emplace(tmp);
    }
    anchors_filled = true;
    logger.trace() << "Added " << anchors.size() << " anchors" << std::endl;
}

EdgePosition KmerIndex::getAnchor(const hashing::KWH &kwh) const {
    if (kwh.isCanonical())
        return anchors.find(kwh.hash())->second;
    else
        return anchors.find(kwh.hash())->second.RC();
}

KmerIndex::KmerIndex(SparseDBG &dbg) : hasher_(dbg.hasher()) {
    for(Vertex &vertex : dbg.verticesUnique()) {
        VERIFY(vertex.isCanonical());
        VERIFY(v.find(vertex.getHash()) == v.end());
        VERIFY(vertex.getSeq().empty() || hashing::MovingKWH(hasher(), vertex.getSeq(), 0).hash() == vertex.getHash());
        v[vertex.getHash()] = &vertex;
    }
}
