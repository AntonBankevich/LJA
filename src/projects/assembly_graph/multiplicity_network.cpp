#include "multiplicity_network.hpp"
#include "common/verify.hpp"
#include <lemon/smart_graph.h>
#include <lemon/network_simplex.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>

using namespace ag;

namespace {
//    Values of the lower convex hull of points (i, values[i]) at every i.
    std::vector<double> LowerConvexHull(const std::vector<double> &values) {
        std::vector<size_t> hull;
        for (size_t i = 0; i < values.size(); i++) {
            while (hull.size() >= 2) {
                size_t a = hull[hull.size() - 2];
                size_t b = hull.back();
//                b is above segment a-i (or on it): drop it
                double cross = (values[b] - values[a]) * double(i - a) - (values[i] - values[a]) * double(b - a);
                if (cross < 0)
                    break;
                hull.pop_back();
            }
            hull.push_back(i);
        }
        std::vector<double> res(values.size());
        for (size_t j = 0; j + 1 < hull.size(); j++) {
            size_t a = hull[j];
            size_t b = hull[j + 1];
            for (size_t i = a; i <= b; i++)
                res[i] = values[a] + (values[b] - values[a]) * double(i - a) / double(b - a);
        }
        res[hull.back()] = values[hull.back()];
        return res;
    }
}

size_t MultiplicityNetwork::addNode() {
    VERIFY(!solved);
    return node_count++;
}

int64_t MultiplicityNetwork::scaled(double cost) const {
    VERIFY(std::isfinite(cost));
    return int64_t(std::llround(cost * cost_scale));
}

size_t MultiplicityNetwork::addArc(size_t from, size_t to, int64_t lower, int64_t upper, int64_t cost, size_t owner) {
    VERIFY(!solved);
    VERIFY(from < node_count && to < node_count);
    arcs.push_back({from, to, lower, upper, cost, owner});
    return arcs.size() - 1;
}

size_t MultiplicityNetwork::addLongArc(size_t from, size_t to, const std::vector<double> &costs) {
    VERIFY(!costs.empty());
    std::vector<double> hull = LowerConvexHull(costs);
    size_t id = long_arcs.size();
    LongArc long_arc{from, to, 0, {}};
    std::vector<int64_t> rounded;
    for (double c : hull)
        rounded.push_back(scaled(c));
    long_arc.base = addArc(from, to, 1, 1, rounded[0], id);
    for (size_t i = 1; i < rounded.size(); i++)
        long_arc.increments.push_back(addArc(from, to, 0, 1, rounded[i] - rounded[i - 1], id));
    long_arcs.push_back(std::move(long_arc));
    return id;
}

size_t MultiplicityNetwork::addFreeArc(size_t from, size_t to) {
    return addArc(from, to, 0, -1, 0, NO_OWNER);
}

size_t MultiplicityNetwork::addStart(size_t node, double cost) {
    return addArc(S, node, 0, -1, scaled(cost), NO_OWNER);
}

size_t MultiplicityNetwork::addEnd(size_t node, double cost) {
    return addArc(node, T, 0, -1, scaled(cost), NO_OWNER);
}

void MultiplicityNetwork::addFixedStart(size_t node, size_t amount) {
    addArc(S, node, int64_t(amount), int64_t(amount), 0, NO_OWNER);
}

void MultiplicityNetwork::addFixedEnd(size_t node, size_t amount) {
    addArc(node, T, int64_t(amount), int64_t(amount), 0, NO_OWNER);
}

bool MultiplicityNetwork::solve() {
    VERIFY(!solved);
//    Every cycle with negative cost passes through a bounded unit arc, so no optimal circulation needs more flow on an
//    arc than all bounded capacities together.
    int64_t big = 1;
    for (const Arc &arc : arcs)
        if (arc.upper >= 0)
            big += arc.upper;
    addArc(T, S, 0, big, 0, NO_OWNER);
    for (Arc &arc : arcs)
        if (arc.upper < 0)
            arc.upper = big;
    solved = true;

    lemon::SmartDigraph g;
    g.reserveNode(int(node_count));
    g.reserveArc(int(arcs.size()));
    std::vector<lemon::SmartDigraph::Node> nodes;
    for (size_t i = 0; i < node_count; i++)
        nodes.push_back(g.addNode());
    std::vector<lemon::SmartDigraph::Arc> lemon_arcs;
    for (const Arc &arc : arcs)
        lemon_arcs.push_back(g.addArc(nodes[arc.from], nodes[arc.to]));
    lemon::SmartDigraph::ArcMap<int64_t> lower(g), upper(g), cost(g);
    for (size_t i = 0; i < arcs.size(); i++) {
        lower[lemon_arcs[i]] = arcs[i].lower;
        upper[lemon_arcs[i]] = arcs[i].upper;
        cost[lemon_arcs[i]] = arcs[i].cost;
    }
    lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, int64_t> simplex(g);
    simplex.lowerMap(lower).upperMap(upper).costMap(cost);
    auto status = simplex.run();
    if (status != lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, int64_t>::OPTIMAL)
        return false;
    flows.clear();
    for (const lemon::SmartDigraph::Arc &arc : lemon_arcs)
        flows.push_back(simplex.flow(arc));
    potentials.clear();
    for (const lemon::SmartDigraph::Node &node : nodes)
        potentials.push_back(simplex.potential(node));
    residual.assign(node_count, {});
    for (size_t i = 0; i < arcs.size(); i++) {
        if (flows[i] < arcs[i].upper)
            residual[arcs[i].from].emplace_back(i, true);
        if (flows[i] > arcs[i].lower)
            residual[arcs[i].to].emplace_back(i, false);
    }
    return true;
}

size_t MultiplicityNetwork::multiplicity(size_t long_arc) const {
    VERIFY(solved && !flows.empty());
    const LongArc &la = long_arcs[long_arc];
    int64_t res = flows[la.base];
    for (size_t arc : la.increments)
        res += flows[arc];
    return size_t(res);
}

size_t MultiplicityNetwork::flow(size_t arc) const {
    VERIFY(solved && !flows.empty());
    return size_t(flows[arc]);
}

double MultiplicityNetwork::totalCost() const {
    VERIFY(solved && !flows.empty());
    int64_t res = 0;
    for (size_t i = 0; i < arcs.size(); i++)
        res += arcs[i].cost * flows[i];
    return double(res) / cost_scale;
}

int64_t MultiplicityNetwork::reducedCost(const Arc &arc) const {
    return arc.cost + potentials[arc.from] - potentials[arc.to];
}

int64_t MultiplicityNetwork::residualDistance(size_t from, size_t to, size_t excluded_owner, int64_t limit) const {
    if (limit < 0)
        return -1;
    if (from == to)
        return 0;
    std::vector<int64_t> dist(node_count, -1);
    typedef std::pair<int64_t, size_t> Entry;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    dist[from] = 0;
    queue.emplace(0, from);
    while (!queue.empty()) {
        auto [d, node] = queue.top();
        queue.pop();
        if (d != dist[node])
            continue;
        if (node == to)
            return d;
        for (const std::pair<size_t, bool> &step : residual[node]) {
            const Arc &arc = arcs[step.first];
            if (arc.owner == excluded_owner)
                continue;
            int64_t reduced = step.second ? reducedCost(arc) : -reducedCost(arc);
            VERIFY(reduced >= 0);
            size_t next = step.second ? arc.to : arc.from;
            int64_t nd = d + reduced;
            if (nd > limit)
                continue;
            if (dist[next] < 0 || nd < dist[next]) {
                dist[next] = nd;
                queue.emplace(nd, next);
            }
        }
    }
    return -1;
}

double MultiplicityNetwork::increaseCost(size_t long_arc, double max_cost) const {
    VERIFY(solved && !flows.empty());
    const LongArc &la = long_arcs[long_arc];
    int64_t best_unit = -1;
    for (size_t arc : la.increments) {
        if (flows[arc] == 0) {
            int64_t reduced = reducedCost(arcs[arc]);
            if (best_unit < 0 || reduced < best_unit)
                best_unit = reduced;
        }
    }
    if (best_unit < 0)
        return INF;
    int64_t limit = scaled(max_cost);
    int64_t rest = residualDistance(la.to, la.from, long_arc, limit - best_unit);
    if (rest < 0)
        return INF;
    return double(best_unit + rest) / cost_scale;
}

double MultiplicityNetwork::decreaseCost(size_t long_arc, double max_cost) const {
    VERIFY(solved && !flows.empty());
    const LongArc &la = long_arcs[long_arc];
    int64_t best_unit = -1;
    for (size_t arc : la.increments) {
        if (flows[arc] == 1) {
            int64_t reduced = -reducedCost(arcs[arc]);
            if (best_unit < 0 || reduced < best_unit)
                best_unit = reduced;
        }
    }
    if (best_unit < 0)
        return INF;
    int64_t limit = scaled(max_cost);
    int64_t rest = residualDistance(la.from, la.to, long_arc, limit - best_unit);
    if (rest < 0)
        return INF;
    return double(best_unit + rest) / cost_scale;
}
