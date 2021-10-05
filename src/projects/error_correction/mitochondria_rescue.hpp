#pragma once
#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_alignment_storage.hpp"

inline bool CheckCov(const Component &component, double &d) {
    size_t bad_cnt = 0;
    size_t good_cnt = 0;
    for(Edge &edge :component.edgesUnique()) {
        if(edge.getCoverage() <= d)
            bad_cnt++;
        else
            good_cnt++;
    }
    return good_cnt >= bad_cnt * 3 / 2;
//    TODO Need more elaborate check
}

inline void MRescue(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                    RecordStorage &reads_storage, size_t unique_length, double error_fraction = 0.05) {
    logger.info() << "Attempting to rescue small circular highly covered components" << std::endl;
    std::unordered_set<Edge const *> bad_edges;
    size_t cnt = 0;
    for(const Component &component : CCSplitter().splitGraph(dbg)) {
        if(component.size() > 100)
            continue;
        bool ok = true;
        std::vector<double> covs;
        for(Edge &edge : component.edgesUnique()) {
            if(edge.size() > unique_length) {
                ok = false;
                break;
            }
            covs.emplace_back(edge.getCoverage());
        }
        if(!ok)
            continue;
        std::sort(covs.begin(), covs.end());
        for(size_t i = 0; i + 1 < covs.size(); i++) {
            if(covs[i] < covs[i + 1] * error_fraction) {
                if(CheckCov(component, covs[i])) {
                    size_t sz = 0;
                    for(Edge &edge :component.edges()) {
                        if(edge.getCoverage() <= covs[i]) {
                            bad_edges.emplace(&edge);
                        } else {
                            sz += edge.size();
                        }
                    }
                    logger.trace() << "Rescued component of size " << (sz / 2) << std::endl;
                    cnt++;
                    break;
                }
            }
        }
    }
    std::function<bool(const Edge&)> is_bad = [&bad_edges](const Edge &edge) {
        return bad_edges.find(&edge) != bad_edges.end();
    };
    reads_storage.invalidateBad(logger, threads, is_bad, "after_mitres");
    logger.info() << "Rescued " << cnt << " circular highly covered components" << std::endl;
}
