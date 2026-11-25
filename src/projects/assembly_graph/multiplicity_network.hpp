#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace ag {

//    Min-cost circulation that assigns multiplicities to long vertices. Graph independent: callers add nodes and arcs.
//
//    Every long vertex is one "long arc" from its in-node to its out-node. Its flow is the multiplicity, at least 1, and
//    costs cost(m) given by the caller (convexified here). Free arcs (cost 0, unbounded) connect long arcs through the
//    short part of the graph. The genome may start or end at any node through start arcs (S -> node) and end arcs
//    (node -> T) with a given cost per unit; T -> S closes the circulation. Fixed arcs carry an exact amount of flow
//    from S or into T: they inject flow of boundary vertices whose multiplicity is already known.
//
//    After solve(), increaseCost/decreaseCost give the minimal extra total cost of a solution where one long arc carries
//    one more / one less unit: the cheapest residual cycle through that arc, found by Dijkstra on reduced costs.
//    All costs are in nats; internally they are rounded to integers in units of 1 / cost_scale.
    class MultiplicityNetwork {
    public:
        static constexpr double INF = std::numeric_limits<double>::infinity();

        explicit MultiplicityNetwork(double cost_scale = 1000.0) : cost_scale(cost_scale) {}

        size_t addNode();
        size_t nodeCount() const {return node_count;}
//        costs[i] is the cost of multiplicity i + 1; the arc can not carry more than costs.size() units. costs must not be
//        empty. Costs are replaced by their lower convex hull so that unit arcs with non-decreasing costs represent them.
        size_t addLongArc(size_t from, size_t to, const std::vector<double> &costs);
//        Cost 0, unbounded capacity. Returns an arc id usable with flow().
        size_t addFreeArc(size_t from, size_t to);
        size_t addStart(size_t node, double cost);
        size_t addEnd(size_t node, double cost);
//        Exactly amount units enter the network at node (from S) or leave it at node (to T).
        void addFixedStart(size_t node, size_t amount);
        void addFixedEnd(size_t node, size_t amount);

//        Returns false if there is no feasible circulation (can not happen when every node with a fixed arc or a long arc
//        endpoint also has start and end arcs).
        bool solve();

        size_t multiplicity(size_t long_arc) const;
        size_t maxMultiplicity(size_t long_arc) const {return long_arcs[long_arc].increments.size() + 1;}
//        Flow on a free, start or end arc (ids returned by the add methods).
        size_t flow(size_t arc) const;
        double totalCost() const;
//        Minimal increase of total cost if the long arc carries one unit more / less. INF if impossible (multiplicity
//        cap reached, multiplicity 1, or no way to reroute) or if the cost exceeds max_cost.
        double increaseCost(size_t long_arc, double max_cost) const;
        double decreaseCost(size_t long_arc, double max_cost) const;

    private:
        struct Arc {
            size_t from;
            size_t to;
            int64_t lower;
            int64_t upper;
            int64_t cost;
            size_t owner;  // long arc id for unit arcs of long arcs, NO_OWNER otherwise
        };
        struct LongArc {
            size_t from;
            size_t to;
            size_t base;  // arc with exactly one unit: multiplicity is at least 1
            std::vector<size_t> increments;  // unit arcs for multiplicities 2, 3, ...; non-decreasing costs
        };
        static constexpr size_t NO_OWNER = size_t(-1);

        size_t addArc(size_t from, size_t to, int64_t lower, int64_t upper, int64_t cost, size_t owner);
        int64_t scaled(double cost) const;
        int64_t reducedCost(const Arc &arc) const;
//        Cheapest path from -> to in the residual graph by reduced costs, not using arcs owned by excluded_owner.
//        Returns -1 if there is none within limit.
        int64_t residualDistance(size_t from, size_t to, size_t excluded_owner, int64_t limit) const;

        double cost_scale;
        size_t node_count = 2;  // 0 = S, 1 = T
        static constexpr size_t S = 0;
        static constexpr size_t T = 1;
        std::vector<Arc> arcs;
        std::vector<LongArc> long_arcs;
        std::vector<int64_t> flows;
        std::vector<int64_t> potentials;
//        node -> (arc, forward) for every arc with residual capacity leaving node; built by solve()
        std::vector<std::vector<std::pair<size_t, bool>>> residual;
        bool solved = false;
    };
}
