#include "correction_utils.hpp"

std::unordered_map<dbg::Vertex *, size_t> findReachable(dbg::Vertex &start, double min_cov, size_t max_dist) {
    typedef std::pair<size_t, dbg::Vertex*> StoredValue;
    std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
    std::unordered_map<dbg::Vertex *, size_t> res;
    queue.emplace(0, &start);
    while(!queue.empty()) {
        StoredValue next = queue.top();
        queue.pop();
        if(res.find(next.second) == res.end()) {
            res[next.second] = next.first;
            for(dbg::Edge &edge : *next.second) {
                size_t new_len = next.first + edge.truncSize();
                if((edge.getData().getCoverage() >= min_cov || edge.getData().is_reliable) && new_len <= max_dist) {
                    queue.emplace(new_len, &edge.getFinish());
                }
            }
        }
    }
    return std::move(res);
}

std::vector<DBGGraphPath>
FindPlausibleBulgeAlternatives(const DBGGraphPath &path, size_t max_diff, double min_cov) {
    size_t max_len = path.truncLen() + max_diff;
    std::unordered_map<dbg::Vertex *, size_t> reachable = findReachable(path.finish().rc(), min_cov, max_len);
    std::vector<DBGGraphPath> res;
    DBGGraphPath alternative(path.start());
    size_t iter_cnt = 0;
    size_t len = 0;
    bool forward = true;
    while(true) {
        iter_cnt += 1;
        if(iter_cnt > 10000)
            return {path};
        if(forward) {
            if(alternative.finish() == path.finish() && len + max_diff >= path.truncLen()) {
                res.emplace_back(alternative);
                if(res.size() > 30) {
                    return {path};
                }
            }
            forward = false;
            for(dbg::Edge &edge : alternative.finish()) {
                if((edge.getData().getCoverage() >= min_cov || edge.getData().is_reliable) && reachable.find(&edge.getFinish().rc()) != reachable.end() &&
                   reachable[&edge.getFinish().rc()] + edge.truncSize() + len <= max_len) {
                    len += edge.truncSize();
                    alternative+=Segment<dbg::Edge>(edge, 0, edge.truncSize());
                    forward = true;
                    break;
                }
            }
        } else {
            if (alternative.size() == 0)
                break;
            dbg::Edge &old_edge = alternative.back().contig();
            alternative.pop_back();
            len -= old_edge.truncSize();
            bool found = false;
            for(dbg::Edge &edge : alternative.finish()) {
                if((edge.getData().getCoverage() >= min_cov || edge.getData().is_reliable) &&
                   reachable.find(&edge.getFinish().rc()) != reachable.end() &&
                   reachable[&edge.getFinish().rc()] + edge.truncSize() + len <= max_len) {
                    if(found) {
                        len += edge.truncSize();
                        alternative+=Segment<dbg::Edge>(edge, 0, edge.truncSize());
                        forward = true;
                        break;
                    } else if (&edge == &old_edge) {
                        found = true;
                    }
                }
            }
        }
    }
    return std::move(res);
}

DBGGraphPath FindReliableExtension(dbg::Vertex &start, size_t len, double min_cov) {
    DBGGraphPath res(start);
    size_t clen = 0;
    while(clen < len) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.finish()) {
            if(edge.getData().is_reliable || edge.getData().getCoverage() >= min_cov) {
                if(next == nullptr)
                    next = &edge;
                else
                    return {};
            }
        }
        if(next == nullptr)
            return {};
        res += *next;
        clen += next->truncSize();
    }
    return std::move(res);
}

std::vector<DBGGraphPath>
FindPlausibleTipAlternatives(const DBGGraphPath &path, size_t max_diff, double min_cov) {
    size_t k = path.start().size();
    size_t max_len = path.truncLen() + max_diff;
    std::vector<DBGGraphPath> res;
    VERIFY(path.leftSkip() == 0);
    DBGGraphPath alternative(path.start());
    size_t iter_cnt = 0;
    size_t len = 0;
    size_t tip_len = path.truncLen();
    bool forward = true;
    while(true) {
        iter_cnt += 1;
        if(iter_cnt > 10000)
            return {path};
        if(forward) {
            forward = false;
            if(len >= tip_len + max_diff) {
                res.emplace_back(alternative);
                if(res.size() > 10) {
                    return {path};
                }
            } else {
                for (dbg::Edge &edge : alternative.finish()) {
                    if (edge.getData().getCoverage() >= min_cov || edge.getData().is_reliable) {
                        len += edge.truncSize();
                        alternative += Segment<dbg::Edge>(edge, 0, edge.truncSize());
                        forward = true;
                        break;
                    }
                }
            }
        } else {
            if (alternative.size() == 0)
                break;
            dbg::Edge &old_edge = alternative.back().contig();
            alternative.pop_back();
            len -= old_edge.truncSize();
            bool found = false;
            for(dbg::Edge &edge : alternative.finish()) {
                if(edge.getData().getCoverage() >= min_cov || edge.getData().is_reliable) {
                    if(found) {
                        len += edge.truncSize();
                        alternative += Segment<dbg::Edge>(edge, 0, edge.truncSize());
                        forward = true;
                        break;
                    } else if (&edge == &old_edge) {
                        found = true;
                    }
                }
            }
        }
    }
    return std::move(res);
}

DBGGraphPath FindLongestCoveredForwardExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov) {
    DBGGraphPath res(start);
    while(true) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.finish()) {
            if(edge == start) {
                return std::move(res);
            }
            if(edge.getData().getCoverage() >= min_rel_cov) {
                if(next == nullptr)
                    next = &edge;
                else {
                    return std::move(res);
                }
            } else if(edge.getData().getCoverage() > max_err_cov) {
                return std::move(res);
            }
        }
        for(dbg::Edge &edge : res.finish().rc()) {
            if(edge != res.back().contig().rc() && edge.getData().getCoverage() > max_err_cov) {
                return std::move(res);
            }
        }
        if(next == nullptr) {
            return res;
        }
        res += *next;
    }
}

DBGGraphPath FindLongestCoveredExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov) {
    DBGGraphPath res = FindLongestCoveredForwardExtension(start, min_rel_cov, max_err_cov);
    if(res.start() == res.finish())
        return std::move(res);
    return FindLongestCoveredForwardExtension(start.rc(), min_rel_cov, max_err_cov).subPath(1).RC() + res;
}
