#pragma once
#include "dbg/sparse_dbg.hpp"

class DatasetParameters {
private:
    double avg_cov;
    double sigma2;
    std::vector<double> density;
    size_t len_step;
    size_t observation_step;
public:
    DatasetParameters(const std::vector<size_t> &observations, size_t len_step) :  len_step(len_step) {
        avg_cov = 0;
        size_t max_val = 0;
        sigma2 = 0;
        for(size_t val : observations) {
            avg_cov += val;
            sigma2 += val * val;
            max_val = std::max(val, max_val);
        }
        avg_cov /= observations.size();
        sigma2 /= observations.size();
        sigma2 -= avg_cov * avg_cov;
        observation_step = std::max<size_t>(1, avg_cov / 20);
        max_val = std::min<size_t>(max_val, max_val / observation_step * 100);
        density = std::vector<double>(max_val / observation_step + 1, 0);
        for(size_t val : observations) {
            density[std::min<size_t>(density.size() - 1, val / observation_step)] += 1;
        }
        for(double &val : density) {
            val /= observations.size();
        }
        VERIFY_MSG(observations.empty() || sigma2 >= -1e-6, sigma2);
    }

    void Print(std::ostream &os) {
        os << "Average: " << avg_cov << " variance2: " << sigma2 << " using length step: " << len_step << "\n";
        for(double cov : {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18.,  20., 25., 30., 40., 50., 60., 80., 100.}) {
            os << "Penalties for " << cov << "\n";
            for (size_t i = 0; i < 20; i++)
                os << getPoissonLogPDiff(len_step, cov, i);
            os << "\n";
            for (size_t i = 0; i < 20; i++)
                os << getNormalLogPDiff(len_step, cov, i);
            os << "\n";
        }
    }

    double getP(double cov) const {
        size_t cell = std::min<size_t>(density.size() - 1, size_t(cov) / observation_step);
        return density[cell];
    }

    double getPoissonLogPDiff(size_t length, double cov, size_t mult) const {
        size_t events = std::max<size_t>(1ull, length / len_step);
        auto observed_coverage = size_t(cov * events);
        double dmult = mult;
        return (log(1 + 1. / double(dmult))* cov - avg_cov) * events;
    }

    double getNormalLogPDiff(size_t length, double cov, size_t mult) const {
        size_t events = std::max<size_t>(1ull, length / len_step);
        double avg1 = events * avg_cov * mult;
        double avg2 = events * avg_cov * (mult + 1);
        double sigma21 = sigma2 * mult * events;
        double sigma22 = sigma2 * (mult + 1) * events;
        double observed_coverage = cov * events;
//        return log(sigma21 / sigma22) / 2 - ((observed_coverage - avg2) * (observed_coverage - avg2) / sigma22 - (observed_coverage - avg1) * (observed_coverage - avg1) / sigma21) / 2;
        return log(mult / (mult + 1)) / 2 - ((cov - avg_cov * (mult + 1)) * (cov - avg_cov * (mult + 1)) / (mult + 1) - (cov - avg_cov * mult) * (cov - avg_cov * mult)  / mult) *events / 2 / sigma2;

    }
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

DatasetParameters EstimateDatasetParameters(dbg::SparseDBG &dbg, const RecordStorage &recordStorage, bool diploid) {
    std::vector<size_t> observations;
    for(dbg::Vertex &v : dbg.vertices()) {
        if(v.inDeg() != 1 || v.outDeg() != 2)
            continue;
        double min_cov= std::min(v[0].getCoverage(), v[1].getCoverage());
        if(min_cov < 4 || min_cov * 10 < v.rc()[0].getCoverage())
            continue;
        observations.emplace_back(size_t(min_cov));
    }
    return {observations, 5000};
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
