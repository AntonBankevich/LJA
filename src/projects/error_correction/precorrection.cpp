#include <assembly_graph/paths.hpp>
#include <dbg/graph_alignment_storage.hpp>
#include "precorrection.hpp"
#include "correction_utils.hpp"
#include "dbg/sparse_dbg.hpp"

dbg::GraphPath FindOnlyPathForward(dbg::Vertex &start, double reliable_coverage, size_t max_size, dbg::Vertex *finish = nullptr) {
    dbg::GraphPath res(start);
    size_t sz = 0;
    while(sz < max_size) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.finish()) {
            if(edge.getCoverage() > 1 && edge.getCoverage() < reliable_coverage) {
                next = nullptr;
                break;
            } else if(edge.getCoverage() >= reliable_coverage) {
                if(next != nullptr) {
                    next = nullptr;
                    break;
                } else {
                    next = &edge;
                }
            }
        }
        if(next == nullptr)
            break;
        size_t len = std::min(max_size - sz, next->truncSize());
        res += Segment<dbg::Edge>(*next, 0, len);
        sz += res.back().size();
        if(&res.finish() == finish)
            break;
    }
    return std::move(res);
}

dbg::GraphPath PrecorrectTip(const Segment<dbg::Edge> &seg, double reliable_coverage) {
    dbg::GraphPath res = FindOnlyPathForward(seg.contig().getStart(), reliable_coverage, seg.size());
    if(res.truncLen() >= seg.size()) {
        res.cutBack(res.truncLen() - seg.size());
        return std::move(res);
    } else {
        return {seg};
    }
}

dbg::GraphPath PrecorrectBulge(dbg::Edge &bulge, double reliable_coverage) {
    dbg::GraphPath res = FindOnlyPathForward(bulge.getStart(), reliable_coverage, bulge.truncSize() + 20,
                                           &bulge.getFinish());
    if(res.finish() == bulge.getFinish() && res.endClosed() && res.truncLen() + 20 > bulge.truncSize()) {
        return std::move(res);
    } else {
        res = FindOnlyPathForward(bulge.getFinish().rc(), reliable_coverage, bulge.truncSize() + 20, &bulge.getStart().rc()).RC();
        if(res.start() == bulge.getStart() && res.startClosed() && res.truncLen() + 20 > bulge.truncSize())
            return std::move(res);
        else {
            std::vector<dbg::GraphPath> candidates = FindPlausibleBulgeAlternatives(dbg::GraphPath(bulge), 10, reliable_coverage);
            if(candidates.size() == 1 && candidates[0].truncLen() + 20 > bulge.truncSize() && candidates[0].truncLen() <
                                                                                            bulge.truncSize() + 20) {
                return std::move(candidates[0]);
            }
            return {bulge};
        }
    }
}


std::string Precorrector::correctRead(dbg::GraphPath &path) {
    if(path.size() == 1)
        return "";
    dbg::GraphPath corrected_path;
    size_t ncor = 0;
    std::vector<std::string> message;
    for(size_t i = 0; i < path.size(); i++) {
        if(path[i].contig().getCoverage() != 1 ||
           (i > 0 && path[i - 1].contig().getCoverage() < reliable_threshold) ||
           (i + 1 < path.size() && path[i + 1].contig().getCoverage() < reliable_threshold)) {
            corrected_path += path[i];
            continue;
        }
        dbg::GraphPath correction;
        std::string m = "";
        if(i == 0) {
            correction = PrecorrectTip(path[i].RC(), reliable_threshold).RC();
            m = "pit";
        } else if(i + 1 == path.size()) {
            correction = PrecorrectTip(path[i], reliable_threshold);
            m = "pot";
        } else {
            correction = PrecorrectBulge(path[i].contig(), reliable_threshold);
            m = "pb";
        }
        if(correction.size() != 1 || correction[0].contig() != path[i].contig()) {
            ncor += 1;
            message.emplace_back(m);
        }
        corrected_path += correction;
    }
    if(!message.empty()) {
        path = corrected_path;
    }
    return join("_", message);
}
