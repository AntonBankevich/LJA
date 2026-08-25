#include "correction_utils.hpp"
using namespace ag;
namespace dbg {
    std::unordered_map<Vertex *, size_t> findReachable(Vertex &start, const std::function<bool(const Edge &)> &isReliable, size_t max_dist) {
        typedef std::pair<size_t, Vertex *> StoredValue;
        std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
        std::unordered_map<Vertex *, size_t> res;
        queue.emplace(0, &start);
        while (!queue.empty()) {
            StoredValue next = queue.top();
            queue.pop();
            if (res.find(next.second) == res.end()) {
                res[next.second] = next.first;
                for (Edge &edge: *next.second) {
                    size_t new_len = next.first + edge.truncSize();
                    if (isReliable(edge) && new_len <= max_dist) {
                        queue.emplace(new_len, &edge.getFinish());
                    }
                }
            }
        }
        return std::move(res);
    }

    std::vector<GraphPath>
    FindPlausibleBulgeAlternatives(const GraphPath &path, size_t max_diff, const std::function<bool(const Edge &)> &isReliable) {
        size_t max_len = path.truncLen() + max_diff;
        std::unordered_map<Vertex *, size_t> reachable = findReachable(path.getFinish().rc(), isReliable, max_len);
        std::vector<GraphPath> res;
        GraphPath alternative(path.getStart());
        size_t iter_cnt = 0;
        size_t len = 0;
        bool forward = true;
        while (true) {
            iter_cnt += 1;
            if (iter_cnt > 10000)
                return {path};
            if (forward) {
                if (alternative.getFinish() == path.getFinish() && len + max_diff >= path.truncLen()) {
                    res.emplace_back(alternative);
                    if (res.size() > 30) {
                        return {path};
                    }
                }
                forward = false;
                for (Edge &edge: alternative.getFinish()) {
                    if (isReliable(edge) &&
                        reachable.find(&edge.getFinish().rc()) != reachable.end() &&
                        reachable[&edge.getFinish().rc()] + edge.truncSize() + len <= max_len) {
                        len += edge.truncSize();
                        alternative += edge;
                        forward = true;
                        break;
                    }
                }
            } else {
                if (alternative.empty())
                    break;
                Edge &old_edge = alternative.back().contig();
                alternative.pop_back();
                len -= old_edge.truncSize();
                bool found = false;
                for (Edge &edge: alternative.getFinish()) {
                    if (isReliable(edge) &&
                        reachable.find(&edge.getFinish().rc()) != reachable.end() &&
                        reachable[&edge.getFinish().rc()] + edge.truncSize() + len <= max_len) {
                        if (found) {
                            len += edge.truncSize();
                            alternative += edge;
                            forward = true;
                            break;
                        } else if (&edge == &old_edge) {
                            found = true;
                        }
                    }
                }
            }
        }
        return oneline::removeValue(res.begin(), res.end(), path);
    }

    std::vector<GraphPath>
    FindPlausibleBulgeAlternatives(const GraphPath &path, size_t max_diff, double min_cov) {
        return FindPlausibleBulgeAlternatives(path, max_diff, [min_cov](const Edge &e) {
            return e.getCoverage() >= min_cov || e.is_reliable;
        });
    }

    GraphPath FindReliableExtension(Vertex &start, size_t len, double min_cov) {
        GraphPath res(start);
        size_t clen = 0;
        while (clen < len) {
            Edge *next = nullptr;
            for (Edge &edge: res.getFinish()) {
                if (edge.is_reliable || edge.getCoverage() >= min_cov) {
                    if (next == nullptr)
                        next = &edge;
                    else
                        return {};
                }
            }
            if (next == nullptr)
                return {};
            res += *next;
            clen += next->truncSize();
        }
        return std::move(res);
    }

    std::vector<GraphPath>
    FindPlausibleTipAlternatives(const GraphPath &path, size_t max_diff, double min_cov) {
        size_t k = path.getStart().size();
        size_t max_len = path.truncLen() + max_diff;
        std::vector<GraphPath> res;
        VERIFY(path.leftCut() == 0);
        GraphPath alternative(path.getStart());
        size_t iter_cnt = 0;
        size_t len = 0;
        size_t tip_len = path.truncLen();
        bool forward = true;
        while (true) {
            iter_cnt += 1;
            if (iter_cnt > 10000)
                return {path};
            if (forward) {
                forward = false;
                if (len >= tip_len + max_diff) {
                    res.emplace_back(alternative);
                    if (res.size() > 10) {
                        return {path};
                    }
                } else {
                    for (Edge &edge: alternative.getFinish()) {
                        if (edge.getCoverage() >= min_cov || edge.is_reliable) {
                            len += edge.truncSize();
                            alternative += edge;
                            forward = true;
                            break;
                        }
                    }
                }
            } else {
                if (alternative.empty())
                    break;
                Edge &old_edge = alternative.back().contig();
                alternative.pop_back();
                len -= old_edge.truncSize();
                bool found = false;
                for (Edge &edge: alternative.getFinish()) {
                    if (edge.getCoverage() >= min_cov || edge.is_reliable) {
                        if (found) {
                            len += edge.truncSize();
                            alternative += edge;
                            forward = true;
                            break;
                        } else if (&edge == &old_edge) {
                            found = true;
                        }
                    }
                }
            }
        }
        return oneline::filter<GraphPath, std::vector<GraphPath>::iterator>(res.begin(), res.end(),
                                            [&path](const GraphPath &other) -> bool {return !other.startsWith(path);});
    }

    GraphPath FindLongestCoveredForwardExtension(Edge &start, size_t max_size, double min_rel_cov, double max_err_cov) {
        GraphPath res(start);
        size_t sz = 0;
        while (sz < max_size) {
            Edge *next = nullptr;
            for (Edge &edge: res.getFinish()) {
                if (edge == start) {
                    return std::move(res);
                }
                if (edge.getCoverage() >= min_rel_cov) {
                    if (next == nullptr)
                        next = &edge;
                    else {
                        return std::move(res);
                    }
                } else if (edge.getCoverage() > max_err_cov) {
                    return std::move(res);
                }
            }
            for (Edge &edge: res.getFinish().rc()) {
                if (edge != res.back().contig().rc() && edge.getCoverage() > max_err_cov) {
                    return std::move(res);
                }
            }
            if (next == nullptr) {
                return res;
            }
            res += *next;
            sz++;
        }
        return std::move(res);
    }

    GraphPath FindLongestCoveredExtension(Edge &start, size_t max_size, double min_rel_cov, double max_err_cov) {
        GraphPath res = FindLongestCoveredForwardExtension(start, max_size, min_rel_cov, max_err_cov);
        if (res.getStart() == res.getFinish())
            return std::move(res);
        ag::GraphPath tmp = FindLongestCoveredForwardExtension(start.rc(), max_size, min_rel_cov, max_err_cov);
        tmp.pop_front();
        return tmp.RC() + res;
    }

    std::vector<GraphPath>
    SuffixSupportedBulgeAlternatives(const ag::SuffixTracker &tracker, const GraphPath &bulge, size_t threshold) {
        Vertex &start = bulge.getStart();
        Vertex &end = bulge.getFinish();
        std::vector<std::pair<Sequence, int>> candidates;
        for(Edge &first_edge : start) {
            const ag::SuffixRecord &record = tracker.getSuffixRecord(first_edge);
            for (const auto &extension: record.getSuffixes()) {
                if (extension.second == 0)
                    continue;
                GraphPath unpacked(first_edge.getFinish(), extension.first);
                GraphPath prefix(first_edge);
                for (Edge &edge: unpacked.edges()) {
                    if (prefix.getFinish() == end) {
                        candidates.emplace_back(Sequence(prefix.getFSplits().begin(), prefix.getFSplits().end()),
                                                extension.second);
                    }
                    prefix += edge;
                }
                if(prefix.getFinish() == end)
                    candidates.emplace_back(Sequence(prefix.getFSplits().begin(), prefix.getFSplits().end()),
                                            extension.second);
            }
        }
//        unlock();
        if (candidates.empty())
            return {};
        std::sort(candidates.begin(), candidates.end());
        candidates.emplace_back(Sequence(), 0);
        std::vector<GraphPath> res;
        size_t cnt = 0;
        for (size_t i = 0; i < candidates.size(); i++) {
            if (i > 0 && candidates[i - 1].first != candidates[i].first) {
                if (cnt >= threshold)
                    res.emplace_back(start, candidates[i - 1].first);
                cnt = 0;
            }
            cnt += candidates[i].second;
        }
        return oneline::removeValue(res.begin(), res.end(), bulge);
    }

    EdgeId SuffixSupportedExtension(const ag::SuffixRecord &record, const GraphPath &start, size_t min_good,
                             size_t max_bad) {
        size_t bad = 0;
        size_t good = 0;
        EdgeId res;
        for(Edge &next : start.getFinish()) {
            size_t cnt = record.countStartsWith(start + next);
            if(cnt <= max_bad) bad++;
            if(cnt >= min_good) {
                good++;
                res = next.getId();
            }
        }
        if (bad + 1 != start.getFinish().outDeg() || good != 1)
            return {};
        return res;
    }

    GraphPath
    FullSuffixSupportedExtension(const ag::SuffixRecord &record, GraphPath start, size_t min_good_cov,
                                 size_t max_bad_cov, size_t max_size) {
        if(!start.valid())
            start = {record.getEdge().getFinish()};
        ag::PathPosition pos = start.lastPosition();
        for(size_t i = 0; i < max_size; i++) {
            EdgeId next = SuffixSupportedExtension(record, start, min_good_cov, max_bad_cov);
            if (!next.valid())
                break;
            start += *next;
        }
        return start.subPath(pos);
    }

    std::vector<GraphPath>
    SuffixSupportedTipAlternatives(const ag::SuffixTracker &tracker, const GraphPath &tip, double threshold) {
        Vertex &start = tip.getStart();
        size_t len = tip.truncLen();
        len += std::max<size_t>(30, len / 20);
        std::vector<std::pair<Sequence, int>> candidates;
        for(Edge &first_edge : start) {
            const ag::SuffixRecord &record = tracker.getSuffixRecord(first_edge);
            if(first_edge.truncSize() >= len) {
                candidates.emplace_back(first_edge.getCode(), record.countStartsWith(GraphPath()));
                continue;
            }
            for (const auto &extension: record.getSuffixes()) {
                if (extension.second == 0)
                    continue;
                GraphPath unpacked(first_edge.getFinish(), extension.first);
                unpacked.push_front(first_edge);
                if (unpacked.truncLen() >= len) {
                    unpacked.cutBack(unpacked.truncLen() - len);
                    candidates.emplace_back(unpacked.getFSplits(), extension.second);
                }
            }
        }
        if (candidates.empty())
            return {};
        std::sort(candidates.begin(), candidates.end());
        candidates.emplace_back(Sequence(), 0);
        std::vector<GraphPath> res;
        size_t cnt = 0;
        for (size_t i = 0; i<candidates.size(); i++) {
            if (i > 0 && candidates[i - 1].first != candidates[i].first) {
                if (cnt > threshold) {
                    GraphPath cp(start, candidates[i - 1].first);
                    cp.cutBack(cp.truncLen() - len);
                    res.emplace_back(cp);
                }
                cnt = 0;
            }
            cnt += candidates[i].second;
        }
        return oneline::filter<GraphPath, std::vector<GraphPath>::iterator>(res.begin(), res.end(),
                                            [&tip](const GraphPath &path) -> bool {return !path.startsWith(tip);});
    }
}