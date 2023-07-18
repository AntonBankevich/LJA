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

DBGGraphPath GraphAligner::align(const Sequence &seq, dbg::Edge *edge_to, size_t pos_to) {
    size_t k = dbg.hasher().getK();
    size_t cur = k;
    DBGGraphPath res;
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

DBGGraphPath GraphAligner::align(const Sequence &seq, const std::string &name) const {
    std::vector<hashing::KWH> kmers = dbg.extractVertexPositions(seq, 1);
    size_t k = dbg.hasher().getK();
    DBGGraphPath res;
    if (kmers.size() == 0) {
        hashing::KWH kwh(dbg.hasher(), seq, 0);
        while (true) {
            if (dbg.isAnchor(kwh.hash())) {
                dbg::EdgePosition pos = dbg.getAnchor(kwh);
                VERIFY(kwh.pos < pos.pos);
                VERIFY(pos.pos + seq.size() - kwh.pos <= pos.edge->truncSize() + k);
                Segment<dbg::Edge> seg(*pos.edge, pos.pos - kwh.pos, pos.pos + seq.size() - kwh.pos - k);
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
    dbg::Vertex *prestart = &dbg.getVertex(kmers.front());
    if (kmers.front().pos > 0) {
        dbg::Vertex &rcstart = prestart->rc();
        if (!rcstart.hasOutgoing(seq[kmers.front().pos - 1] ^ 3)) {
            std::cout << "No outgoing for getStart" << std::endl << seq << std::endl <<
                      kmers.front().pos << " " << seq[kmers.front().pos - 1] << std::endl
                      << kmers.front().getSeq() << std::endl;
            VERIFY(false);
        }
        dbg::Edge &rcedge = rcstart.getOutgoing(seq[kmers.front().pos - 1] ^ 3);
        dbg::Edge &edge = rcedge.rc();
        VERIFY(edge.truncSize() >= kmers.front().pos);
        Segment<dbg::Edge> seg(edge, edge.truncSize() - kmers.front().pos, edge.truncSize());
        res += seg;
    }
    size_t cpos = kmers.front().pos + k;
    while(cpos < seq.size()) {
        if(!prestart->hasOutgoing(seq[cpos])) {
            std::cout << "No outgoing for middle\n" << seq << "\n" << cpos << " " << prestart->getInnerId() <<
                      " " << prestart->outDeg() << "\n" << seq.Subseq(cpos -k, cpos) << "\n" << prestart->getSeq() << "\n" <<
                      size_t(seq[cpos]) << std::endl;
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

DBGGraphPath GraphAligner::align(const dbg::EdgePosition &pos, const Sequence &seq) const {
    DBGGraphPath res(Segment<dbg::Edge>(*pos.edge, pos.pos, pos.pos));
    dbg::Edge *cedge = pos.edge;
    size_t epos = pos.pos;
    for (size_t cpos = 0; cpos < seq.size(); cpos++) {
        unsigned char c = seq[cpos];
        if (epos == cedge->truncSize()) {
            dbg::Vertex &v = cedge->getFinish();
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

std::vector<dbg::PerfectAlignment<dbg::Edge, dbg::Edge>> GraphAligner::oldEdgeAlign(dbg::Edge &contig) const {
    Sequence seq = contig.getSeq();
    std::vector<PerfectAlignment < dbg::Edge, dbg::Edge>> res;
    hashing::KWH kwh(dbg.hasher(), seq, 0);
    size_t k = dbg.hasher().getK();
    while (true) {
        if (!kwh.hasNext())
            break;
        if (res.empty() || kwh.pos >= res.back().seg_from.right) {
            dbg::Edge *edge = nullptr;
            size_t pos = 0;
            if (dbg.containsVertex(kwh.hash())) {
                dbg::Vertex &start = dbg.getVertex(kwh);
                if (start.hasOutgoing(seq[kwh.pos + k]))
                    edge = &dbg.getVertex(kwh).getOutgoing(seq[kwh.pos + k]);
            }
            if (edge == nullptr && dbg.isAnchor(kwh.hash())) {
                dbg::EdgePosition gpos = dbg.getAnchor(kwh);
                if (gpos.edge->truncSeq()[gpos.pos] == seq[kwh.pos + k]) {
                    edge = gpos.edge;
                    pos = gpos.pos;
                }
            }
            if (edge != nullptr) {
                size_t len = std::min(contig.truncSize() - kwh.pos, edge->truncSize() - pos);
                res.emplace_back(Segment<dbg::Edge>(contig, kwh.pos, kwh.pos + len), Segment<dbg::Edge>(*edge, pos, pos + len));
            }
        }
        kwh = kwh.next();
    }
    return std::move(res);
}

std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> dbg::GraphAligner::carefulAlign(Contig &contig) const {
    Sequence seq = contig.getSeq();
    size_t k = dbg.hasher().getK();
    if(contig.truncSize() < k) {
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
                    while (len < edge.truncSize() && len < kwh.pos && edge.truncSeq()[len] == (seq[kwh.pos - len - 1] ^ 3))
                        len += 1;
                    res.emplace_back(Segment<Contig>(contig, kwh.pos - len, kwh.pos),
                                     Segment<Edge>(edge.rc(), edge.truncSize() - len, edge.truncSize()));
                }
                if (kwh.pos + k < seq.size() && vertex.hasOutgoing(seq[kwh.pos + k])) {
                    Edge &edge = vertex.getOutgoing(seq[kwh.pos + k]);
                    size_t len = 1;
                    while (len < edge.truncSize() && kwh.pos + k + len < seq.size() &&
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
    PerfectAlignment<Contig, Edge> start_al = bestExtension(start.rc(), Segment<Contig>(rc_contig, contig.truncSize() - k - kwh.pos, contig.truncSize() - k));
    return {Segment<Contig>(contig, contig.truncSize() - k - start_al.seg_from.right, contig.truncSize() - k - start_al.seg_from.left),
            Segment<Edge>(start_al.seg_to.contig().rc(),
                          start_al.seg_to.contig().truncSize() - start_al.seg_to.right,
                          start_al.seg_to.contig().truncSize() - start_al.seg_to.left)};
}

PerfectAlignment<Contig, dbg::Edge> GraphAligner::extendRight(const hashing::KWH &kwh, Contig &contig) const {
    size_t k = dbg.hasher().getK();
    PerfectAlignment<Contig, dbg::Edge> best({contig, kwh.pos, kwh.pos}, {});
    if(kwh.pos + k == contig.truncSize()) {
        return best;
    }
    Vertex &start = dbg.getVertex(kwh);
    if(start.outDeg() == 0) {
        return best;
    }
    return bestExtension(start, Segment<Contig>(contig, kwh.pos, contig.truncSize() - k));
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
