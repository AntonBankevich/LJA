#include <assembly_graph/random_access_paths.hpp>
#include <dbg/dbg_read_alignment_storage.hpp>
#include "precorrection.hpp"
#include "correction_utils.hpp"
#include "dbg/sparse_dbg.hpp"

dbg::GraphPath FindOnlyPathForward(dbg::Vertex &start, double reliable_coverage, size_t max_size, dbg::Vertex *finish = nullptr) {
    dbg::GraphPath res(start);
    size_t sz = 0;
    while(sz < max_size) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.getFinish()) {
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
        if(&res.getFinish() == finish)
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
bool isSimplestBulge(dbg::Vertex &start, dbg::Vertex &finish) {
    return start.outDeg() == 2 && finish.inDeg() == 2 && start.front().getFinish() == finish &&
        start.back().getFinish() == finish && start.front().getCoverage() == 1 && start.back().getCoverage() == 1;
}

dbg::GraphPath PrecorrectBulge(dbg::Edge &bulge, double reliable_coverage) {
    dbg::GraphPath res = FindOnlyPathForward(bulge.getStart(), reliable_coverage, bulge.truncSize() + 20,
                                           &bulge.getFinish());
    if(res.getFinish() == bulge.getFinish() && res.endClosed() && res.truncLen() + 20 > bulge.truncSize()) {
        return std::move(res);
    } else {
        res = FindOnlyPathForward(bulge.getFinish().rc(), reliable_coverage, bulge.truncSize() + 20, &bulge.getStart().rc()).RC();
        if(res.getStart() == bulge.getStart() && res.startClosed() && res.truncLen() + 20 > bulge.truncSize())
            return std::move(res);
        else {
            std::vector<dbg::GraphPath> candidates = FindPlausibleBulgeAlternatives(dbg::GraphPath(bulge), 10, reliable_coverage);
            if(candidates.size() == 1 && candidates[0].truncLen() + 20 > bulge.truncSize() && candidates[0].truncLen() <
                                                                                            bulge.truncSize() + 20) {
                return std::move(candidates[0]);
            }
            return dbg::GraphPath(bulge);
        }
    }
}


std::string Precorrector::correctRead(const std::string &name, dbg::GraphPath &path) {
    if(path.isSingleton())
        return "";
    dbg::GraphPath corrected_path;
    size_t ncor = 0;
    std::vector<std::string> message;
    for(dbg::PathPosition pp = path.firstPosition(); pp != path.lastPosition(); ++pp) {
        dbg::PathPosition ppp1 = pp + 1;
        if(pp.nextEdge().getCoverage() != 1 ||
           (pp != path.firstPosition() && pp.prevEdge().getCoverage() < reliable_threshold) ||
           (ppp1 != path.lastPosition() && ppp1.nextEdge().getCoverage() < reliable_threshold)) {
            corrected_path += path.getSegment(pp);
            continue;
        }
        dbg::GraphPath correction;
        std::string m = "";
        if(pp == path.firstPosition()) {
            correction = PrecorrectTip(path.front().RC(), reliable_threshold).RC();
            m = "pit";
        } else if(ppp1 == path.lastPosition()) {
            correction = PrecorrectTip(path.back(), reliable_threshold);
            m = "pot";
        } else {
            if(isSimplestBulge(pp.getVertex(), ppp1.getVertex())) {
                dbg::Edge &other = pp.getVertex().front() == pp.nextEdge() ? pp.getVertex().back() : pp.getVertex().front();
                dbg::EdgeId other_canonical = other.rc().getId() < other.getId() ? other.rc().getId() : other.getId();
                dbg::EdgeId cur = pp.nextEdge().getId();
                if(cur->rc().getId() < cur) cur = cur->rc().getId();
                if(other_canonical < cur) {
                    correction = {other};
                    m = "p1b";
                } else {
                    correction = {pp.nextEdge()};
                }
            } else {
                correction = PrecorrectBulge(pp.nextEdge(), reliable_threshold);
                m = "pb";
            }
        }
        if(!correction.isSingleton() || correction.front() != pp.nextEdge()) {
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
