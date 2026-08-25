#include <assembly_graph/random_access_paths.hpp>
#include <dbg/dbg_read_alignment_storage.hpp>
#include "precorrection.hpp"
#include "correction_utils.hpp"
#include "dbg/sparse_dbg.hpp"

ag::GraphPath FindOnlyPathForward(dbg::Vertex &start, const std::function<bool(const ag::Edge &)> &isReliable,
                                   const std::function<bool(const ag::Edge &)> &isSuspicious,
                                   size_t max_size, dbg::Vertex *finish = nullptr) {
    ag::GraphPath res(start);
    size_t sz = 0;
    while(sz < max_size) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.getFinish()) {
            if(!isSuspicious(edge) && !isReliable(edge)) {
                next = nullptr;
                break;
            } else if(isReliable(edge)) {
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

ag::GraphPath PrecorrectTip(const Segment<dbg::Edge> &seg, const std::function<bool(const ag::Edge &)> &isReliable,
                             const std::function<bool(const ag::Edge &)> &isSuspicious) {
    ag::GraphPath res = FindOnlyPathForward(seg.contig().getStart(), isReliable, isSuspicious, seg.size());
    if(res.truncLen() >= seg.size()) {
        res.cutBack(res.truncLen() - seg.size());
        return std::move(res);
    } else {
        return {seg};
    }
}
bool isSimplestBulge(dbg::Vertex &start, dbg::Vertex &finish, const std::function<bool(const ag::Edge &)> &isSuspicious) {
    return start.outDeg() == 2 && finish.inDeg() == 2 && start.front().getFinish() == finish &&
        start.back().getFinish() == finish && isSuspicious(start.front()) && isSuspicious(start.back());
}

ag::GraphPath PrecorrectBulge(dbg::Edge &bulge, const std::function<bool(const ag::Edge &)> &isReliable,
                               const std::function<bool(const ag::Edge &)> &isSuspicious) {
    ag::GraphPath res = FindOnlyPathForward(bulge.getStart(), isReliable, isSuspicious, bulge.truncSize() + 20,
                                           &bulge.getFinish());
    if(res.getFinish() == bulge.getFinish() && res.endClosed() && res.truncLen() + 20 > bulge.truncSize()) {
        return std::move(res);
    } else {
        res = FindOnlyPathForward(bulge.getFinish().rc(), isReliable, isSuspicious, bulge.truncSize() + 20, &bulge.getStart().rc()).RC();
        if(res.getStart() == bulge.getStart() && res.startClosed() && res.truncLen() + 20 > bulge.truncSize())
            return std::move(res);
        else {
            std::vector<ag::GraphPath> candidates = dbg::FindPlausibleBulgeAlternatives(ag::GraphPath(bulge), 10, isReliable);
            if(candidates.size() == 1 && candidates[0].truncLen() + 20 > bulge.truncSize() && candidates[0].truncLen() <
                                                                                            bulge.truncSize() + 20) {
                return std::move(candidates[0]);
            }
            return ag::GraphPath(bulge);
        }
    }
}


std::string Precorrector::correctRead(const std::string &name, ag::GraphPath &path) {
    if(path.isSingleton())
        return "";
    ag::GraphPath corrected_path;
    size_t ncor = 0;
    std::vector<std::string> message;
    for(ag::PathPosition pp = path.firstPosition(); pp != path.lastPosition(); ++pp) {
        ag::PathPosition ppp1 = pp + 1;
        if(!isSuspicious(pp.nextEdge()) ||
           (pp != path.firstPosition() && !isReliable(pp.prevEdge())) ||
           (ppp1 != path.lastPosition() && !isReliable(ppp1.nextEdge()))) {
            corrected_path += path.getSegment(pp);
            continue;
        }
        ag::GraphPath correction;
        std::string m = "";
        if(pp == path.firstPosition()) {
            correction = PrecorrectTip(path.front().RC(), isReliable, isSuspicious).RC();
            m = "pit";
        } else if(ppp1 == path.lastPosition()) {
            correction = PrecorrectTip(path.back(), isReliable, isSuspicious);
            m = "pot";
        } else {
            if(isSimplestBulge(pp.getVertex(), ppp1.getVertex(), isSuspicious)) {
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
                correction = PrecorrectBulge(pp.nextEdge(), isReliable, isSuspicious);
                m = "pb";
            }
        }
        if(!correction.isSingleton() || correction.frontEdge() != pp.nextEdge()) {
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
