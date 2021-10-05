#include "correction_utils.hpp"

dbg::Edge *checkBorder(dbg::Vertex &v) {
    dbg::Edge * res = nullptr;
    size_t out_rel = 0;
    for(dbg::Edge &edge : v.rc()) {
        if(edge.is_reliable) {
            if (res == nullptr)
                res = &edge.rc();
            else
                return nullptr;
        }
    }
    for(dbg::Edge &edge : v) {
        if(edge.is_reliable)
            out_rel += 1;
    }
    if(out_rel == 0)
        return res;
    else
        return nullptr;
}

bool checkInner(dbg::Vertex &v) {
    size_t inc_rel = 0;
    size_t out_rel = 0;
    for(dbg::Edge &edge : v.rc()) {
        if(edge.is_reliable)
            inc_rel += 1;
    }
    for(dbg::Edge &edge : v) {
        if(edge.is_reliable)
            out_rel += 1;
    }
    return out_rel >= 1 && inc_rel >= 1;
}

void FillReliableWithConnections(logging::Logger &logger, dbg::SparseDBG &sdbg, double threshold) {
    logger.info() << "Marking reliable edges" << std::endl;
    for(auto &vit : sdbg) {
        for(dbg::Vertex * vp : {&vit.second, &vit.second.rc()}) {
            dbg::Vertex &v = *vp;
            for(dbg::Edge &edge : v) {
                edge.is_reliable = edge.getCoverage() >= threshold;
            }
        }
    }
    size_t cnt_paths = 0;
    std::vector<dbg::Edge *> new_reliable;
    for(auto &vit : sdbg) {
        for(dbg::Vertex * vp : {&vit.second, &vit.second.rc()}) {
            dbg::Vertex &v = *vp;
            dbg::Edge *last = checkBorder(v);
            if(last == nullptr)
                continue;
            typedef std::pair<double, dbg::Edge *> StoredValue;
            std::unordered_map<dbg::Vertex *, std::pair<double, dbg::Edge *>> res;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            queue.emplace(0, last);
            size_t cnt = 0;
            while(!queue.empty()) {
                cnt += 1;
                if(cnt > 10000) {
//                    logger << "Dijkstra too long" << std::endl;
                    break;
                }
                dbg::Edge *next = queue.top().second;
                double dist = queue.top().first;
                queue.pop();
                if(res.find(next->end()) != res.end() || checkInner(*next->end()))
                    continue;
                res.emplace(next->end(), std::make_pair(dist, next));
                if (checkBorder(next->end()->rc()) != nullptr) {
                    dbg::GraphAlignment al(next->end()->rc());
                    while(next != last) {
                        al.push_back(Segment<dbg::Edge>(next->rc(), 0, next->size()));
                        new_reliable.emplace_back(next);
//                        logger << "New edge marked as reliable " << next->size() << "(" << next->getCoverage() << ")" << std::endl;
                        next = res[next->start()].second;
                    }
//                    logger << "Path of size " << al.size() << " and length " << al.len() << " marked as reliable." << std::endl;
                    cnt_paths += 1;
                    break;
                } else {
                    for(dbg::Edge &edge : *next->end()) {
                        double score = edge.size() / std::max<double>(std::min(edge.getCoverage(), threshold), 1);
                        queue.emplace(dist + score, &edge);
                    }
                }
            }
        }
    }
    logger.info() << "Marked " << new_reliable.size() << " edges in " << cnt_paths << " paths as reliable" << std::endl;
    for(dbg::Edge *edge : new_reliable) {
        edge->is_reliable = true;
        edge->rc().is_reliable = true;
    }
}

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
