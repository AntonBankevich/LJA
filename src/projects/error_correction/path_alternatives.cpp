#include "path_alternatives.hpp"
#include "common/oneline_utils.hpp"
#include "sequences/edit_distance.hpp"
#include <queue>

namespace ag {
    std::unordered_map<Vertex *, size_t> findReachable(Vertex &start,
                                                        const std::function<bool(const Edge &)> &isTraversable,
                                                        size_t max_dist) {
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
                    size_t new_len = next.first + edge.rc().truncSize();
                    if (isTraversable(edge) && new_len <= max_dist) {
                        queue.emplace(new_len, &edge.getFinish());
                    }
                }
            }
        }
        return std::move(res);
    }

    std::vector<GraphPath> FindAlternativeSegments(const GraphPath &path, size_t max_diff,
                                                    const std::function<bool(const Edge &)> &isTraversable) {
        size_t max_flen = path.truncLen() + max_diff;
        size_t flen = path.truncLen();
        std::unordered_map<Vertex *, size_t> reachable = findReachable(path.getFinish().rc(), isTraversable, max_flen);
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
                if (alternative.getFinish() == path.getFinish() && len + max_diff >= flen) {
                    res.emplace_back(alternative);
                    if (res.size() > 30) {
                        return {path};
                    }
                }
                forward = false;
                for (Edge &edge: alternative.getFinish()) {
                    if (isTraversable(edge) &&
                        reachable.find(&edge.getFinish().rc()) != reachable.end() &&
                        reachable[&edge.getFinish().rc()] + edge.truncSize() + len <= max_flen) {
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
                    if (isTraversable(edge) &&
                        reachable.find(&edge.getFinish().rc()) != reachable.end() &&
                        reachable[&edge.getFinish().rc()] + edge.truncSize() + len <= max_flen) {
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

    std::vector<GraphPath> FindAlternativeTips(const GraphPath &path, size_t max_diff,
                                                const std::function<bool(const Edge &)> &isTraversable) {
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
                        if (isTraversable(edge)) {
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
                    if (isTraversable(edge)) {
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

    size_t tournament(const Sequence &original, const std::vector<Sequence> &candidates) {
        size_t winner = 0;
        std::vector<size_t> dists;
        size_t max_dist = std::max<size_t>(20, original.size() / 100);
        for (size_t i = 0; i < candidates.size(); i++) {
            dists.push_back(edit_distance(original, candidates[i], max_dist));
            if (dists.back() < dists[winner])
                winner = i;
        }
        if (dists[winner] >= max_dist)
            return -1;
        for (size_t i = 0; i < candidates.size(); i++) {
            if (i != winner && dists[i] < max_dist) {
                size_t diff = edit_distance(candidates[winner], candidates[i], max_dist);
                VERIFY(dists[winner] <= dists[i] + diff);
                VERIFY(dists[i] <= dists[winner] + diff);
                if (dists[i] < max_dist && dists[i] != dists[winner] + diff)
                    return -1;
            }
        }
        return winner;
    }

    std::pair<GraphPath, size_t> bestAlignmentPrefix(const GraphPath &al, const Sequence &seq, size_t max_diff) {
        Sequence candSeq = al.truncSeq();
        VERIFY(!candSeq.empty());
        std::pair<size_t, size_t> bp = bestPrefix(seq, candSeq, max_diff);
        size_t len = bp.first;
        VERIFY(len > 0);
        Sequence prefix = candSeq.Subseq(0, len);
        GraphPath res(al.getStart());
        res.extend(prefix);
        return {res, bp.second};
    }
}
