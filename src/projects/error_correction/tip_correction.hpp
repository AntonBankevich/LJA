#pragma once

#include "edit_distance.hpp"

inline void MakeUnreliable(Edge &e) {
    e.is_reliable = false;
    for(Edge &edge : *e.end()) {
        if(edge.is_reliable) {
            edge.is_reliable = false;
            MakeUnreliable(edge);
        }
    }
}

inline void FillReliableTips(logging::Logger &logger, dbg::SparseDBG &sdbg, double reliable_threshold) {
    logger.info() << "Remarking reliable edges" << std::endl;
    for(auto &vit : sdbg) {
        for(Vertex * vp : {&vit.second, &vit.second.rc()}) {
            Vertex &v = *vp;
            for(Edge &edge : v) {
                edge.is_reliable = true;
            }
        }
    }
    size_t infty = 1000000000;
    std::unordered_map<Vertex *, size_t> max_tip;
    std::vector<Edge*> queue;
    for(Edge &edge : sdbg.edges()) {
        if(edge.end()->outDeg() == 0 && edge.end()->inDeg() == 1 && edge.size() < 10000 && edge.getCoverage() < reliable_threshold) {
            max_tip[edge.end()] = 0;
            queue.emplace_back(&edge);
        }
    }
    while(!queue.empty()) {
        Edge &new_edge = *queue.back();
        queue.pop_back();
        Vertex &v = *new_edge.start();
        bool good = true;
        size_t val = 0;
        Edge *best = nullptr;
        for(Edge &out : v) {
            if(max_tip.find(out.end()) == max_tip.end()) {
                good = false;
                break;
            } else {
                if(val < max_tip[out.end()] + out.size()) {
                    val = max_tip[out.end()] + out.size();
                    best = &out;
                }
            }
        }
        if(good && val < 30000) {
            max_tip[&v] = val;
            for(Edge &out : v) {
                if(&out == best) {
                    out.is_reliable = true;
                } else {
                    MakeUnreliable(out);
                }
            }
            if(v.inDeg() == 1 && v.rc().begin()->size() < 10000 && v.rc().begin()->getCoverage() < reliable_threshold) {
                queue.emplace_back(&(v.rc().begin()->rc()));
            }
        }
    }
}

inline Path ReliablePath(Vertex &v, size_t max_size = 1000000000) {
    Path path(v);
    size_t len = 0;
    while(len < max_size) {
        Edge *next = nullptr;
        for (Edge &edge : path.finish()) {
            if (edge.is_reliable) {
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
        len += next->size();
    }
    return path;
}



inline GraphAlignment CorrectSuffix(const GraphAlignment &al) {
    size_t first_unreliable = al.size();
    size_t bad_end_size = 0;
    while(first_unreliable > 0 && !al[first_unreliable - 1].contig().is_reliable) {
        first_unreliable--;
        bad_end_size += al[first_unreliable].size();
    }
    if(first_unreliable == al.size()) {
        return al;
    }
    size_t max_len = bad_end_size  * 11/10 + 100;
    Path alternative = ReliablePath(al.getVertex(first_unreliable), max_len);
    if(alternative.finish().outDeg() != 0 && alternative.len() + 10 < bad_end_size) {
        return al;
    }
    Sequence tip = al.truncSeq(first_unreliable);
    Sequence alt = alternative.truncSeq();
    Sequence projection = alt.Subseq(0, bestPrefix(tip, alt).first);
    GraphAlignment res = al.subalignment(0, first_unreliable);
    res.extend(projection);
    return res;
}

inline void CorrectTips(logging::Logger &logger, size_t threads, SparseDBG &dbg,
                        const std::vector<RecordStorage *> &storages) {
    logger.info() << "Correcting tips using reliable edge marks" << std::endl;
    omp_set_num_threads(threads);
    ParallelCounter cnt(threads);
    for(RecordStorage *storageIt : storages) {
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(storageIt, cnt)
        for (size_t i = 0; i < storageIt->size(); i++) {
            AlignedRead &read = storageIt->operator[](i);
            if (!read.valid())
                continue;
            GraphAlignment al = read.path.getAlignment();
            GraphAlignment al1 = CorrectSuffix(al);
            GraphAlignment al2 = CorrectSuffix(al1.RC()).RC();
            if (al != al2) {
                cnt += 1;
                storageIt->reroute(read, al, al2, "Tip corrected");
            }
        }
        storageIt->applyCorrections(logger, threads);
    }
    logger.info() << "Corrected tips for " << cnt.get() << " reads" << std::endl;
}

inline void TipCorrectionPipeline(logging::Logger &logger, SparseDBG &dbg, RecordStorage &reads,
                           size_t threads,double reliable_threshold) {
    FillReliableTips(logger, dbg, reliable_threshold);
    CorrectTips(logger, threads, dbg, {&reads});
}
