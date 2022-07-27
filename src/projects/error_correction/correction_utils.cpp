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
                size_t new_len = next.first + edge.size();
                if((edge.getCoverage() >= min_cov || edge.is_reliable) && new_len <= max_dist) {
                    queue.emplace(new_len, edge.end());
                }
            }
        }
    }
    return std::move(res);
}

std::vector<dbg::GraphAlignment>
FindPlausibleBulgeAlternatives(const dbg::GraphAlignment &path, size_t max_diff, double min_cov) {
    size_t k = path.start().seq.size();
    size_t max_len = path.len() + max_diff;
    std::unordered_map<dbg::Vertex *, size_t> reachable = findReachable(path.finish().rc(), min_cov, max_len);
    std::vector<dbg::GraphAlignment> res;
    dbg::GraphAlignment alternative(path.start());
    size_t iter_cnt = 0;
    size_t len = 0;
    bool forward = true;
    while(true) {
        iter_cnt += 1;
        if(iter_cnt > 10000)
            return {path};
        if(forward) {
            if(alternative.finish() == path.finish() && len + max_diff >= path.len()) {
                res.emplace_back(alternative);
                if(res.size() > 30) {
                    return {path};
                }
            }
            forward = false;
            for(dbg::Edge &edge : alternative.finish()) {
                if((edge.getCoverage() >= min_cov || edge.is_reliable) && reachable.find(&edge.end()->rc()) != reachable.end() &&
                   reachable[&edge.end()->rc()] + edge.size() + len <= max_len) {
                    len += edge.size();
                    alternative.push_back(Segment<dbg::Edge>(edge, 0, edge.size()));
                    forward = true;
                    break;
                }
            }
        } else {
            if (alternative.size() == 0)
                break;
            dbg::Edge &old_edge = alternative.back().contig();
            alternative.pop_back();
            len -= old_edge.size();
            bool found = false;
            for(dbg::Edge &edge : alternative.finish()) {
                if((edge.getCoverage() >= min_cov || edge.is_reliable) &&
                   reachable.find(&edge.end()->rc()) != reachable.end() &&
                   reachable[&edge.end()->rc()] + edge.size() + len <= max_len) {
                    if(found) {
                        len += edge.size();
                        alternative.push_back(Segment<dbg::Edge>(edge, 0, edge.size()));
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

dbg::GraphAlignment FindReliableExtension(dbg::Vertex &start, size_t len, double min_cov) {
    dbg::GraphAlignment res(start);
    size_t clen = 0;
    while(clen < len) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.finish()) {
            if(edge.is_reliable || edge.getCoverage() >= min_cov) {
                if(next == nullptr)
                    next = &edge;
                else
                    return {};
            }
        }
        if(next == nullptr)
            return {};
        res += *next;
        clen += next->size();
    }
    return std::move(res);
}

std::vector<dbg::GraphAlignment>
FindPlausibleTipAlternatives(const dbg::GraphAlignment &path, size_t max_diff, double min_cov) {
    size_t k = path.start().seq.size();
    size_t max_len = path.len() + max_diff;
    std::vector<dbg::GraphAlignment> res;
    dbg::GraphAlignment alternative(path.start());
    size_t iter_cnt = 0;
    size_t len = 0;
    size_t tip_len = path.len();
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
                    if (edge.getCoverage() >= min_cov || edge.is_reliable) {
                        len += edge.size();
                        alternative.push_back(Segment<dbg::Edge>(edge, 0, edge.size()));
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
            len -= old_edge.size();
            bool found = false;
            for(dbg::Edge &edge : alternative.finish()) {
                if(edge.getCoverage() >= min_cov || edge.is_reliable) {
                    if(found) {
                        len += edge.size();
                        alternative.push_back(Segment<dbg::Edge>(edge, 0, edge.size()));
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

dbg::GraphAlignment FindLongestCoveredForwardExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov) {
    dbg::GraphAlignment res(start);
    while(true) {
        dbg::Edge *next = nullptr;
        for(dbg::Edge &edge : res.finish()) {
            if(edge == start) {
                return std::move(res);
            }
            if(edge.getCoverage() >= min_rel_cov) {
                if(next == nullptr)
                    next = &edge;
                else {
                    return std::move(res);
                }
            } else if(edge.getCoverage() > max_err_cov) {
                return std::move(res);
            }
        }
        for(dbg::Edge &edge : res.finish().rc()) {
            if(edge != res.back().contig().rc() && edge.getCoverage() > max_err_cov) {
                return std::move(res);
            }
        }
        if(next == nullptr) {
            return res;
        }
        res += *next;
    }
}

dbg::GraphAlignment FindLongestCoveredExtension(dbg::Edge &start, double min_rel_cov, double max_err_cov) {
    dbg::GraphAlignment res = FindLongestCoveredForwardExtension(start, min_rel_cov, max_err_cov);
    if(res.start() == res.finish())
        return std::move(res);
    return FindLongestCoveredForwardExtension(start.rc(), min_rel_cov, max_err_cov).subalignment(1).RC() + res;
}
