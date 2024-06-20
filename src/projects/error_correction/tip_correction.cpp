#include "tip_correction.hpp"
#include "common/logging.hpp"

using namespace dbg;
void MakeUnreliable(Edge &e) {
    e.is_reliable = false;
    for(Edge &edge : e.getFinish()) {
        if(edge.is_reliable) {
            edge.is_reliable = false;
            MakeUnreliable(edge);
        }
    }
}

inline void FillReliableTips(logging::Logger &logger, dbg::SparseDBG &sdbg, double reliable_threshold) {
    logger.info() << "Remarking reliable edges" << std::endl;
    for(Edge &edge : sdbg.edges()) edge.is_reliable = true;
    size_t infty = 1000000000;
    std::unordered_map<VertexId, size_t> max_tip;
    std::vector<EdgeId> queue;
    for(Edge &edge : sdbg.edges()) {
        if(edge.getFinish().outDeg() == 0 && edge.getFinish().inDeg() == 1 && edge.truncSize() < 15000 &&
                edge.getCoverage() < reliable_threshold) {
            max_tip[edge.getFinish().getId()] = 0;
            queue.emplace_back(edge.getId());
        }
    }
    while(!queue.empty()) {
        Edge &new_edge = *queue.back();
        queue.pop_back();
        Vertex &v = new_edge.getStart();
        bool good = true;
        size_t val = 0;
        EdgeId best;
        for(Edge &out : v) {
            if(max_tip.find(out.getFinish().getId()) == max_tip.end()) {
                good = false;
                break;
            } else {
                if(val < max_tip[out.getFinish().getId()] + out.truncSize()) {
                    val = max_tip[out.getFinish().getId()] + out.truncSize();
                    best = out.getId();
                }
            }
        }
        if(good && val < 30000) {
            max_tip[v.getId()] = val;
            for(Edge &out : v) {
                if(out.getId() == best) {
                    out.is_reliable = true;
                } else {
                    MakeUnreliable(out);
                }
            }
            if(v.inDeg() == 1 && v.rc().begin()->truncSize() < 10000 &&
                    v.rc().begin()->getCoverage() < reliable_threshold) {
                queue.emplace_back(v.rc().begin()->rc().getId());
            }
        }
    }
}

inline dbg::GraphPath ReliablePath(Vertex &v, size_t max_size = 1000000000) {
    dbg::GraphPath path(v);
    size_t len = 0;
    while(len < max_size) {
        EdgeId next;
        for (Edge &edge : path.getFinish()) {
            if (edge.is_reliable) {
                if(next.valid()) {
                    next = EdgeId();
                    break;
                }
                next = edge.getId();
            }
        }
        if(!next.valid())
            break;
        path += *next;
        len += next->truncSize();
    }
    return path;
}



inline dbg::GraphPath CorrectSuffix(const dbg::GraphPath &al) {
    PathPosition first_unreliable = al.lastPosition();
    size_t bad_end_size = 0;
    while(first_unreliable != al.firstPosition() && !first_unreliable.prevEdge().is_reliable) {
        bad_end_size += first_unreliable.prevEdge().truncSize();
        --first_unreliable;
    }
    if(first_unreliable == al.lastPosition() || first_unreliable == al.firstPosition()) {
        return al;
    }
    size_t max_len = bad_end_size  * 11/10 + 100;
    dbg::GraphPath alternative = ReliablePath(first_unreliable.getVertex(), max_len);
    if(alternative.getFinish().outDeg() != 0 && alternative.truncLen() + 100 < bad_end_size) {
        return al;
    }
    Sequence tip = al.subPath(first_unreliable).truncSeq();
    Sequence alt = alternative.truncSeq();
    Sequence projection = alt;
    if(alt.size() > tip.size())
        projection = alt.Subseq(0, bestPrefix(tip, alt).first);
    dbg::GraphPath res = al.subPath(al.firstPosition(), first_unreliable);
    res.extend(projection);
    return res;
}

size_t CorrectTips(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                 std::vector<DBGAlignedReadStorage *> storages) {
    logger.info() << "Correcting tips using reliable edge marks" << std::endl;
    omp_set_num_threads(threads);
    ParallelCounter cnt(threads);
    for(dbg::DBGAlignedReadStorage *storageIt : storages) {
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storageIt, cnt)
        for (size_t i = 0; i < storageIt->getReads().size(); i++) {
            ag::AlignedRead<DBGTraits> &read = storageIt->getReads()[i];
            if (!read.valid())
                continue;
            dbg::GraphPath al = read.getPath();
            dbg::GraphPath al1 = CorrectSuffix(al);
            dbg::GraphPath al2 = CorrectSuffix(al1.RC()).RC();
            if (al != al2 && al2.truncLen() > 500) {
                cnt += 1;
                storageIt->getReads().rerouteRead(read, al2, "Tip corrected");
            }
        }
        storageIt->getReads().applyCorrections(logger, threads);
    }
    return cnt.get();
}

void TipCorrectionPipeline(logging::Logger &logger, dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads, size_t threads,
                           double reliable_threshold) {
    FillReliableTips(logger, dbg, reliable_threshold);
    size_t cnt = CorrectTips(logger, threads, dbg, {&reads});
    for(Edge &edge : dbg.edges())
        edge.is_reliable = false;
    logger.info() << "Corrected " << cnt << " reads using algorithm " << "TipCorrector" << std::endl;
}
