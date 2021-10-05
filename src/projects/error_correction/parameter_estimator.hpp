#pragma once
#include "dbg/sparse_dbg.hpp"

class DatasetParameters {
public:
    bool diploid;
    double coverage_median;
    double coverage5;
    double coverage95;
};

std::vector<const dbg::Edge *> GetOutgoing(const dbg::Vertex &start, double min_cov) {
    std::vector<const dbg::Edge *> res;
    for(const dbg::Edge &edge : start) {
        if(edge.getCoverage() >= min_cov)
            res.emplace_back(&edge);
    }
    return std::move(res);
}

std::vector<const dbg::Edge *> CoveredPath(const dbg::Edge &start, double min_cov, size_t max_len) {
    std::vector<const dbg::Edge *> res = {&start};
    size_t len = start.size();
    while(len < max_len) {
        std::vector<const dbg::Edge *> out = GetOutgoing(*res.back()->end(), min_cov);
        if(out.size() == 1) {
            res.emplace_back(out[0]);
            len += out[0]->size();
        }
    }
    return std::move(res);
}

const dbg::Edge *SeekCovered(const dbg::Vertex &start, double min_cov) {
    const dbg::Vertex *tmp = &start;
    size_t max_steps = 20;
    for(size_t i = 0; i < max_steps; i++) {
        for(const dbg::Edge &edge : *tmp) {
            if (edge.getCoverage() >= min_cov)
                return &edge;
        }
        for(const dbg::Edge &edge : *tmp) {
            if (edge.end()->outDeg() > 0) {
                tmp = edge.end();
                break;
            }
        }
    }
    return nullptr;
}

//DatasetParameters EstimateDatasetParameters(const dbg::SparseDBG &dbg) {
//    size_t num_votes = 1000;
//    size_t hap_vote_len = 20000;
//    double min_cov = 2;
//    size_t hap_votes = 0;
//    size_t dip_votes = 0;
//    size_t failed_votes = 0;
//    size_t cnt = 0;
//    for(auto & rec : dbg) {
//        if(hap_votes + dip_votes + failed_votes >= num_votes)
//            break;
//        const dbg::Vertex &v = rec.second;
//        const dbg::Edge *first_covered =SeekCovered(v, min_cov);
//        if(first_covered == nullptr) {
//            failed_votes += 1;
//            continue;
//        }
//        size_t clen = 0;
//        const dbg::Vertex *cur = first_covered->start();
//        bool is_hap = true;
//        while(clen < hap_vote_len) {
//            std::vector<const dbg::Edge *> out = GetOutgoing(*cur, min_cov);
//            if(out.size() > 2 || out.empty())
//                break;
//            if(out.size() == 1) {
//                std::vector<const dbg::Edge *> path = CoveredPath(*out[0], min_cov, hap_vote_len);
//                for(const dbg::Edge * edge : path) {
//                    clen += edge->size();
//                }
//                cur = path.back()->end();
//            } else {
//                std::vector<const dbg::Edge *> path1 = CoveredPath(*out[0], min_cov, hap_vote_len);
//                std::vector<const dbg::Edge *> path2 = CoveredPath(*out[1], min_cov, hap_vote_len);
//                if(path1.back()->end() == path2.back()->end()) {
//                    is_hap = false;
//
//                }
//            }
//        }
//    }
//    VERIFY(failed_votes * 4 < num_votes);
//}