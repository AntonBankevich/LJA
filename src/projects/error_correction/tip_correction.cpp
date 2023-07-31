#include "tip_correction.hpp"
#include "common/logging.hpp"

using namespace dbg;
void MakeUnreliable(Edge &e) {
    e.getData().is_reliable = false;
    for(Edge &edge : e.getFinish()) {
        if(edge.getData().is_reliable) {
            edge.getData().is_reliable = false;
            MakeUnreliable(edge);
        }
    }
}

inline void FillReliableTips(logging::Logger &logger, dbg::SparseDBG &sdbg, double reliable_threshold) {
    logger.info() << "Remarking reliable edges" << std::endl;
    for(auto &vit : sdbg.verticesUnique()) {
        for(Vertex * vp : {&vit, &vit.rc()}) {
            Vertex &v = *vp;
            for(Edge &edge : v) {
                edge.getData().is_reliable = true;
            }
        }
    }
    size_t infty = 1000000000;
    std::unordered_map<Vertex *, size_t> max_tip;
    std::vector<Edge*> queue;
    for(Edge &edge : sdbg.edges()) {
        if(edge.getFinish().outDeg() == 0 && edge.getFinish().inDeg() == 1 && edge.truncSize() < 15000 &&
                edge.getData().getCoverage() < reliable_threshold) {
            max_tip[&edge.getFinish()] = 0;
            queue.emplace_back(&edge);
        }
    }
    while(!queue.empty()) {
        Edge &new_edge = *queue.back();
        queue.pop_back();
        Vertex &v = new_edge.getStart();
        bool good = true;
        size_t val = 0;
        Edge *best = nullptr;
        for(Edge &out : v) {
            if(max_tip.find(&out.getFinish()) == max_tip.end()) {
                good = false;
                break;
            } else {
                if(val < max_tip[&out.getFinish()] + out.truncSize()) {
                    val = max_tip[&out.getFinish()] + out.truncSize();
                    best = &out;
                }
            }
        }
        if(good && val < 30000) {
            max_tip[&v] = val;
            for(Edge &out : v) {
                if(&out == best) {
                    out.getData().is_reliable = true;
                } else {
                    MakeUnreliable(out);
                }
            }
            if(v.inDeg() == 1 && v.rc().begin()->truncSize() < 10000 &&
                    v.rc().begin()->getData().getCoverage() < reliable_threshold) {
                queue.emplace_back(&(v.rc().begin()->rc()));
            }
        }
    }
}

inline DBGGraphPath ReliablePath(Vertex &v, size_t max_size = 1000000000) {
    DBGGraphPath path(v);
    size_t len = 0;
    while(len < max_size) {
        Edge *next = nullptr;
        for (Edge &edge : path.finish()) {
            if (edge.getData().is_reliable) {
                if(next != nullptr) {
                    next = nullptr;
                    break;
                }
                next = &edge;
            }
        }
        if(next == nullptr)
            break;
        path += *next;
        len += next->truncSize();
    }
    return path;
}



inline DBGGraphPath CorrectSuffix(const DBGGraphPath &al) {
    size_t first_unreliable = al.size();
    size_t bad_end_size = 0;
    while(first_unreliable > 0 && !al[first_unreliable - 1].contig().getData().is_reliable) {
        first_unreliable--;
        bad_end_size += al[first_unreliable].size();
    }
    if(first_unreliable == al.size() || first_unreliable == 0) {
        return al;
    }
    size_t max_len = bad_end_size  * 11/10 + 100;
    DBGGraphPath alternative = ReliablePath(al.getVertex(first_unreliable), max_len);
    if(alternative.finish().outDeg() != 0 && alternative.len() + 2000 < bad_end_size) {
        return al;
    }
    Sequence tip = al.truncSubseq(first_unreliable);
    Sequence alt = alternative.truncSeq();
    Sequence projection = alt;
    if(alt.size() > tip.size())
        projection = alt.Subseq(0, bestPrefix(tip, alt).first);
    DBGGraphPath res = al.subPath(0, first_unreliable);
    res.extend(projection);
    return res;
}

void CorrectTips(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                        const std::vector<RecordStorage *> &storages) {
    logger.info() << "Correcting tips using reliable getEdge marks" << std::endl;
    omp_set_num_threads(threads);
    ParallelCounter cnt(threads);
    for(RecordStorage *storageIt : storages) {
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storageIt, cnt)
        for (size_t i = 0; i < storageIt->size(); i++) {
            AlignedRead &read = storageIt->operator[](i);
            if (!read.valid())
                continue;
            DBGGraphPath al = read.path.getAlignment();
            DBGGraphPath al1 = CorrectSuffix(al);
            DBGGraphPath al2 = CorrectSuffix(al1.RC()).RC();
            if (al != al2) {
                cnt += 1;
                storageIt->reroute(read, al, al2, "Tip corrected");
            }
        }
        storageIt->applyCorrections(logger, threads);
    }
}

void TipCorrectionPipeline(logging::Logger &logger, dbg::SparseDBG &dbg, RecordStorage &reads, size_t threads,
                           double reliable_threshold) {
    FillReliableTips(logger, dbg, reliable_threshold);
    CorrectTips(logger, threads, dbg, {&reads});
}
