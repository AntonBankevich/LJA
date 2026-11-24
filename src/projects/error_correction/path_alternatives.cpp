#include "path_alternatives.hpp"
#include "common/oneline_utils.hpp"
#include "sequences/edit_distance.hpp"
#include <algorithm>
#include <map>
#include <queue>
#include <set>

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

    namespace {
        typedef std::function<bool(const Edge &)> EdgePredicate;
        const size_t NO_COST = size_t(-1);

//        Costs of aligning a path's sequence to the query prefixes [0, j), for j in [lo, lo + cost.size()).
        struct AlignmentBand {
            size_t lo = 0;
            std::vector<size_t> cost;

            bool empty() const { return cost.empty(); }
            size_t min() const { return empty() ? NO_COST : *std::min_element(cost.begin(), cost.end()); }
//            The first query position with the minimal cost.
            size_t argmin() const { return lo + (std::min_element(cost.begin(), cost.end()) - cost.begin()); }
            size_t hi() const { return lo + cost.size(); }
            size_t at(size_t j) const { return j >= lo && j < hi() ? cost[j - lo] : NO_COST; }

//            Drops the costs at either end that are at least min() + slack or above max_cost; the band becomes empty
//            if all costs are above max_cost.
            void trim(size_t slack, size_t max_cost) {
                const size_t limit = std::min(min() + slack, max_cost + 1);
                size_t first = 0, last = cost.size();
                while (first < last && cost[first] >= limit)
                    first++;
                while (last > first && cost[last - 1] >= limit)
                    last--;
                cost = std::vector<size_t>(cost.begin() + first, cost.begin() + last);
                lo += first;
            }
        };

//        The band after one more path letter c: the standard edit-distance DP row, trimmed to costs below row minimum
//        + slack and at most max_cost. Beyond the previous band only query letters can be added (horizontal moves),
//        so the row stops once those reach either limit.
        AlignmentBand nextRow(const AlignmentBand &prev, unsigned char c, const Sequence &query, size_t slack,
                              size_t max_cost) {
            AlignmentBand row{prev.lo, {}};
            size_t row_min = NO_COST;
            for (size_t j = prev.lo; j <= query.size(); j++) {
                size_t best = NO_COST;
                if (j < prev.hi())
                    best = prev.at(j) + 1;  // c not in the query
                if (j > prev.lo && j - 1 < prev.hi())
                    best = std::min(best, prev.at(j - 1) + (query[j - 1] != c ? 1 : 0));
                if (!row.cost.empty())
                    best = std::min(best, row.cost.back() + 1);  // query letter not on the path
                if (j >= prev.hi() && (best >= row_min + slack || best > max_cost))
                    break;
                row.cost.push_back(best);
                row_min = std::min(row_min, best);
            }
            row.trim(slack, max_cost);
            return row;
        }

//        Where alternatives to a bulge end. F is the bulge's finish vertex, L(v) the length from the end of v to the end
//        of F along a route; v overlaps F if L(v) < |F|. Past a vertex overlapping F, a route to F only spells F.
        struct BulgeEnds {
//            How a route ends through a final edge u -> x (u does not overlap F, x does): F's sequence starts after
//            `letters` letters of the edge, and tail (x -> F) completes the route.
            struct Ending {
                size_t letters;
                std::vector<Edge *> tail;
            };
            std::unordered_map<const Edge *, std::vector<Ending>> endings;  // final edge -> one tail per letters value
//            Shortest length from the end of a vertex to the start of F's sequence.
            std::unordered_map<Vertex *, size_t> to_finish_start;
//            Routes start -> F along which the start overlaps F, with L(start): they spell the last L(start) letters of F.
            std::vector<std::pair<size_t, std::vector<Edge *>>> overlapping_routes;

            const std::vector<Ending> *endingsOf(const Edge &edge) const {
                auto it = endings.find(&edge);
                return it == endings.end() ? nullptr : &it->second;
            }

            size_t toFinishStart(Vertex &v) const {
                auto it = to_finish_start.find(&v);
                return it == to_finish_start.end() ? NO_COST : it->second;
            }
        };

//        Multi-source Dijkstra backwards: lowers dist[v] to dist[x] + |e| over traversable edges e = v -> x, keeping
//        only distances up to max_dist.
        void relaxBackwards(std::unordered_map<Vertex *, size_t> &dist, const EdgePredicate &isTraversable,
                            size_t max_dist) {
            typedef std::pair<size_t, Vertex *> Entry;
            std::priority_queue<Entry, std::vector<Entry>, std::greater<>> queue;
            for (const auto &[v, d]: dist)
                queue.emplace(d, v);
            while (!queue.empty()) {
                auto [d, x] = queue.top();
                queue.pop();
                if (d > dist[x])
                    continue;
                for (Edge &in: x->incoming()) {
                    const size_t len = d + in.truncSize();
                    Vertex *v = &in.getStart();
                    if (!isTraversable(in) || len > max_dist)
                        continue;
                    auto it = dist.find(v);
                    if (it == dist.end() || len < it->second) {
                        dist[v] = len;
                        queue.emplace(len, v);
                    }
                }
            }
        }

//        Walks back from F over the vertices overlapping F (all their routes, as (vertex, L) states), collecting the
//        final edges and the routes from start; then extends the lengths to the start of F's sequence backwards from the
//        final edges, up to max_len.
        BulgeEnds findBulgeEnds(Vertex &finish, Vertex &start, const EdgePredicate &isTraversable, size_t max_len) {
            BulgeEnds res;
            const size_t f = finish.size();
            struct Overlapping {
                Vertex *v;
                size_t len;               // L(v)
                std::vector<Edge *> tail; // v -> F
            };
            std::vector<Overlapping> overlapping = {{&finish, 0, {}}};
            std::set<std::pair<Vertex *, size_t>> seen = {{&finish, 0}};
            for (size_t i = 0; i < overlapping.size(); i++) {
                const Overlapping x = overlapping[i];  // a copy: the vector grows below
                for (Edge &in: x.v->incoming()) {
                    if (!isTraversable(in))
                        continue;
                    Vertex &u = in.getStart();
                    const size_t len = x.len + in.truncSize();
                    if (len >= f) {  // u does not overlap F: in is a final edge
                        const size_t letters = len - f;
                        std::vector<BulgeEnds::Ending> &ends = res.endings[&in];
                        if (std::none_of(ends.begin(), ends.end(),
                                         [&](const BulgeEnds::Ending &e) { return e.letters == letters; }))
                            ends.push_back({letters, x.tail});
                        if (letters < res.toFinishStart(u))
                            res.to_finish_start[&u] = letters;
                    } else if (seen.emplace(&u, len).second) {
                        std::vector<Edge *> tail = {&in};
                        tail.insert(tail.end(), x.tail.begin(), x.tail.end());
                        if (&u == &start)
                            res.overlapping_routes.emplace_back(len, tail);
                        overlapping.push_back({&u, len, std::move(tail)});
                    }
                }
            }
            relaxBackwards(res.to_finish_start, isTraversable, max_len);
            return res;
        }

//        Best-first search over the tree of paths from a start vertex, each path aligned to a query (the core of
//        FindAlignedAlternatives). A path's priority is its band minimum; all costs are alignment cost + base_cost.
//        Without bulge_ends (tip) an alignment of the whole query may end after any letter, keeping the best one per
//        tree node; with them (bulge) only at the endings of final edges, and paths only continue into vertices that
//        can still reach the start of the finish's sequence.
        class PathAligner {
        public:
//            An alignment of the whole query ending after `letters` letters of the last edge of tree node `node`, then
//            completed by tail (bulge only).
            struct End {
                size_t node;
                size_t letters;
                const std::vector<Edge *> *tail;
                size_t cost;
            };

            PathAligner(Vertex &start, Sequence query, size_t base_cost, size_t max_cost,
                        const AlignedSearchParams &params, const EdgePredicate &isTraversable,
                        const BulgeEnds *bulge_ends) :
                    start(start), query(std::move(query)), base_cost(base_cost), max_cost(max_cost),
                    slack(params.band_slack), max_choices(params.max_choices), max_steps(params.max_steps),
                    divergence(params.divergence), isTraversable(isTraversable), bulge_ends(bulge_ends) {
                VERIFY(slack > 0);
                AlignmentBand empty_path;
                for (size_t j = 0; j <= this->query.size() && j < slack; j++)
                    empty_path.cost.push_back(j);
                nodes.push_back({size_t(-1), nullptr, std::move(empty_path), 0});
            }

//            Runs the search. best_known is the cost of an alternative found otherwise (NO_COST if none). Stops once
//            the next path costs more than max_cost or best + slack; returns false if it gives up after max_steps.
            bool run(size_t best_known) {
                best = best_known;
                queue.emplace(cost(nodes[0].band.min()), 0);
                for (size_t steps = 0; !queue.empty(); steps++) {
                    auto [score, id] = queue.top();
                    queue.pop();
                    if (score > max_cost || (best != NO_COST && score > best + slack))
                        return true;
                    if (steps == max_steps)
                        return false;
                    if (!queryUsedUp(nodes[id].band))
                        for (Edge &edge: id == 0 ? start : nodes[id].edge->getFinish())
                            if (isTraversable(edge))
                                addChild(id, edge);
                }
                return true;
            }

            size_t bestCost() const { return best; }
            const std::vector<End> &ends() const { return found; }

//            The ends that are the cheapest on their route (a path and its extensions are one route; earliest on ties).
            std::vector<End> cheapestPerRoute() const {
                std::vector<size_t> node_cost(nodes.size(), NO_COST);
                for (const End &e: found)
                    node_cost[e.node] = e.cost;
                std::vector<bool> beaten(nodes.size(), false);
                for (const End &e: found) {
                    for (size_t x = nodes[e.node].parent; x != 0; x = nodes[x].parent) {
                        if (node_cost[x] <= e.cost)
                            beaten[e.node] = true;
                        else if (node_cost[x] != NO_COST)
                            beaten[x] = true;
                    }
                }
                std::vector<End> res;
                for (const End &e: found)
                    if (!beaten[e.node])
                        res.push_back(e);
                return res;
            }

//            The edges of the path of a tree node.
            std::vector<Edge *> route(size_t node) const {
                std::vector<Edge *> res;
                for (size_t x = node; x != 0; x = nodes[x].parent)
                    res.push_back(nodes[x].edge);
                std::reverse(res.begin(), res.end());
                return res;
            }

        private:
            struct Node {
                size_t parent;
                Edge *edge;
                AlignmentBand band;  // at the end of edge
                size_t len;          // length of the path's truncSeq
            };

            Vertex &start;
            const Sequence query;
            const size_t base_cost, max_cost, slack, max_choices, max_steps;
            const DivergenceLimit divergence;
            const EdgePredicate &isTraversable;
            const BulgeEnds *bulge_ends;
            std::vector<Node> nodes;
            std::vector<End> found;
            size_t best = NO_COST;
            typedef std::pair<size_t, size_t> QueueEntry;  // (cost, node)
            std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<>> queue;
//            (end vertex, query position of the band minimum) -> the paths kept for extension.
            std::map<std::pair<const Vertex *, size_t>, std::vector<size_t>> choices;

//            Whether the new path `id` ending at v with this band is kept for extension: if it is among the
//            max_choices kept paths with the same key with the lowest band minimum, or if at some query position it is
//            better than all of them.
            bool keepChoice(const Vertex &v, const AlignmentBand &band, size_t id) {
                std::vector<size_t> &kept = choices[{&v, band.argmin()}];
                const size_t band_min = band.min();
                size_t at_least_as_good = 0;
                for (size_t k: kept)
                    if (nodes[k].band.min() <= band_min)
                        at_least_as_good++;
                bool keep = at_least_as_good < max_choices;
                for (size_t j = band.lo; !keep && j < band.hi(); j++) {
                    keep = true;
                    for (size_t k: kept)
                        if (nodes[k].band.at(j) <= band.at(j))
                            keep = false;
                }
                if (keep)
                    kept.push_back(id);
                return keep;
            }

            size_t cost(size_t alignment_cost) const {
                return alignment_cost == NO_COST ? NO_COST : alignment_cost + base_cost;
            }

//            Tip: once only the whole-query cost is left, extending the path only adds to it.
            bool queryUsedUp(const AlignmentBand &band) const {
                return bulge_ends == nullptr && band.lo == query.size();
            }

//            Cost 0 means the original sequence itself: kept (it is dropped at the end), but never the best cost the
//            other candidates are measured against.
            void addEnd(const End &end) {
                if (end.cost > max_cost)
                    return;
                found.push_back(end);
                if (end.cost > 0)
                    best = std::min(best, end.cost);
            }

//            Adds the tree node id + edge: runs the DP along the edge, records the alignment ends on it, and queues the
//            new node if paths continue from it.
            void addChild(size_t id, Edge &edge) {
                const size_t len = nodes[id].len + edge.truncSize();
                const std::vector<BulgeEnds::Ending> *endings = nullptr;
                bool continues = true;
                if (bulge_ends != nullptr) {
                    endings = bulge_ends->endingsOf(edge);
                    const size_t rest = bulge_ends->toFinishStart(edge.getFinish());
                    continues = rest != NO_COST && len + rest <= query.size() + max_cost - base_cost;
                    if (!continues && endings == nullptr)
                        return;
                }
//                Letters of the edge to align: all of them, or, if paths stop here, up to the last ending.
                size_t letters = edge.truncSize();
                if (!continues) {
                    letters = 0;
                    for (const BulgeEnds::Ending &e: *endings)
                        letters = std::max(letters, e.letters);
                }
                const size_t child = nodes.size();
                const size_t n = query.size();
                AlignmentBand band = nodes[id].band;
                size_t tip_cost = NO_COST, tip_letters = 0;  // tip: the best end on this edge
                for (size_t k = 0;; k++) {
                    if (endings != nullptr)
                        for (const BulgeEnds::Ending &e: *endings)
                            if (e.letters == k)
                                addEnd({child, k, &e.tail, cost(band.at(n))});
                    if (k == letters || cost(band.min()) > max_cost || queryUsedUp(band))
                        break;
                    band = nextRow(band, edge.truncSeq()[k], query, slack, divergence.at(nodes[id].len + k + 1));
                    if (bulge_ends == nullptr && band.at(n) < tip_cost) {
                        tip_cost = band.at(n);
                        tip_letters = k + 1;
                    }
                }
                if (tip_cost != NO_COST)
                    addEnd({child, tip_letters, nullptr, cost(tip_cost)});
                const size_t score = cost(band.min());
                if (continues && score <= max_cost && keepChoice(edge.getFinish(), band, child))
                    queue.emplace(score, child);
                nodes.push_back({id, &edge, std::move(band), len});
            }
        };

        GraphPath makePath(Vertex &start, const std::vector<Edge *> &route, size_t cut_right = 0) {
            GraphPath res(start);
            for (Edge *edge: route)
                res += *edge;
            res.setCutRight(cut_right);
            return res;
        }

//        Drops the candidates with cost 0 (the original itself, or a path spelling exactly its sequence) and those
//        costing more than best + slack.
        std::vector<AlignedAlternative> finalCandidates(std::vector<AlignedAlternative> candidates, size_t best,
                                                        size_t slack) {
            candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [&](const AlignedAlternative &a) {
                return a.cost == 0 || a.cost > best + slack;
            }), candidates.end());
            return candidates;
        }

    }

    std::vector<AlignedAlternative> FindAlignedTipAlternatives(const GraphPath &tip, size_t max_cost,
                                                               const EdgePredicate &isTraversable,
                                                               const AlignedSearchParams &params) {
        VERIFY(tip.leftCut() == 0);
        PathAligner aligner(tip.getStart(), tip.truncSeq(), 0, max_cost, params, isTraversable, nullptr);
        if (!aligner.run(NO_COST))
            return {};
        std::vector<AlignedAlternative> res;
        for (const PathAligner::End &e: aligner.cheapestPerRoute()) {
            std::vector<Edge *> route = aligner.route(e.node);
            res.push_back({makePath(tip.getStart(), route, route.back()->truncSize() - e.letters), e.cost});
        }
        return finalCandidates(std::move(res), aligner.bestCost(), params.band_slack);
    }

//    A route from the bulge's start to its finish F spells inner part + F, or, if the start overlaps F along it, a
//    suffix of F. Only inner parts are aligned. If one of two such sequences is a suffix of F, it is a suffix of the
//    other one too, so their edit distance is the length difference; for an overlapping bulge vs an alternative with
//    inner part I that is |I| + base_cost, the cost of aligning I to an empty query.
    std::vector<AlignedAlternative> FindAlignedBulgeAlternatives(const GraphPath &bulge, size_t max_cost,
                                                                 const EdgePredicate &isTraversable,
                                                                 const AlignedSearchParams &params) {
        VERIFY(bulge.leftCut() == 0 && bulge.rightCut() == 0);
        Vertex &start = bulge.getStart();
        Vertex &finish = bulge.getFinish();
        const Sequence seq = bulge.truncSeq();
        const size_t f = finish.size();
        const bool overlap = seq.size() < f;
        const size_t base_cost = overlap ? f - seq.size() : 0;
        if (base_cost > max_cost)
            return {};
        const Sequence inner = overlap ? Sequence() : seq.Subseq(0, seq.size() - f);
        const BulgeEnds ends = findBulgeEnds(finish, start, isTraversable, inner.size() + max_cost - base_cost);
        std::vector<AlignedAlternative> res;
        size_t best = NO_COST;
        for (const auto &[len, route]: ends.overlapping_routes) {
            const size_t cost = std::max(len, seq.size()) - std::min(len, seq.size());
            if (cost > max_cost)
                continue;
            res.push_back({makePath(start, route), cost});
            if (cost > 0)
                best = std::min(best, cost);
        }
        PathAligner aligner(start, inner, base_cost, max_cost, params, isTraversable, &ends);
        if (!aligner.run(best))
            return {};
        for (const PathAligner::End &e: aligner.ends()) {
            std::vector<Edge *> route = aligner.route(e.node);
            route.insert(route.end(), e.tail->begin(), e.tail->end());
            res.push_back({makePath(start, route), e.cost});
        }
        return finalCandidates(std::move(res), aligner.bestCost(), params.band_slack);
    }

    Tournament::Tournament(double divergence_factor, double min_diff_factor, DivergenceLimit limit) :
            divergence_factor(divergence_factor), min_diff_factor(min_diff_factor), limit(limit) {}

    size_t Tournament::maxDist(size_t original_size) const {
        return limit.at(original_size);
    }

    size_t Tournament::tournamentOnDistances(size_t original_size, const std::vector<size_t> &dists,
                                             const std::vector<Sequence> &candidates) const {
        const size_t max_dist = maxDist(original_size);
        const size_t winner = std::min_element(dists.begin(), dists.end()) - dists.begin();
        if (dists[winner] >= max_dist)
            return -1;
        const double error_rate = double(dists[winner]) / double(std::max<size_t>(original_size, 1));
        const double factor = std::max(min_diff_factor, 1 - divergence_factor * error_rate);
        for (size_t i = 0; i < candidates.size(); i++) {
            if (i != winner && dists[i] < max_dist) {
                size_t diff = edit_distance(candidates[winner], candidates[i], max_dist);
                VERIFY(dists[winner] <= dists[i] + diff);
                VERIFY(dists[i] <= dists[winner] + diff);
                if (double(dists[winner]) > double(dists[i]) - double(diff) * factor + 1e-9)
                    return -1;
            }
        }
        return winner;
    }

    size_t Tournament::tournament(const Sequence &original, const std::vector<Sequence> &candidates) const {
        const size_t max_dist = maxDist(original.size());
        std::vector<size_t> dists;
        for (const Sequence &candidate: candidates)
            dists.push_back(edit_distance(original, candidate, max_dist));
        return tournamentOnDistances(original.size(), dists, candidates);
    }

    size_t Tournament::tournament(size_t original_size, const std::vector<AlignedAlternative> &alternatives) const {
        const size_t max_dist = maxDist(original_size);
        std::vector<size_t> dists;
        std::vector<Sequence> candidates;
        for (const AlignedAlternative &a: alternatives) {
            dists.push_back(std::min(a.cost, max_dist));
            candidates.push_back(a.path.truncSeq());
        }
        return tournamentOnDistances(original_size, dists, candidates);
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
