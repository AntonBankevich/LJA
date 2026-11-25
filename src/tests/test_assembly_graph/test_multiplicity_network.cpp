// MultiplicityNetwork: multiplicities from convex per-arc costs, genome ends, and the +-1 cost queries.

#include "gtest/gtest.h"
#include "assembly_graph/multiplicity_network.hpp"
#include <cmath>

using ag::MultiplicityNetwork;

namespace {
//    Quadratic cost around an expected multiplicity, like a coverage likelihood.
    std::vector<double> Costs(double expected, size_t max_mult = 6, double weight = 10) {
        std::vector<double> res;
        for (size_t m = 1; m <= max_mult; m++)
            res.push_back(weight * (double(m) - expected) * (double(m) - expected));
        return res;
    }

    struct LongArc {
        size_t in;
        size_t out;
        size_t arc;
    };

    LongArc AddLong(MultiplicityNetwork &net, double expected, double end_cost) {
        size_t in = net.addNode();
        size_t out = net.addNode();
        size_t arc = net.addLongArc(in, out, Costs(expected));
        net.addStart(in, end_cost);
        net.addEnd(out, end_cost);
        return {in, out, arc};
    }
}

//    A, B -> X -> D, E: X is a repeat of multiplicity 2 between four unique vertices.
TEST(MultiplicityNetwork, RepeatBetweenUniqueVertices) {
    MultiplicityNetwork net;
    double end = 5, brk = 15;
    LongArc a = AddLong(net, 1, end), b = AddLong(net, 1, end), x = AddLong(net, 2, brk);
    LongArc d = AddLong(net, 1.1, end), e = AddLong(net, 0.9, end);
    size_t hub_in = net.addNode(), hub_out = net.addNode();
    net.addFreeArc(a.out, hub_in);
    net.addFreeArc(b.out, hub_in);
    net.addFreeArc(hub_in, x.in);
    net.addFreeArc(x.out, hub_out);
    net.addFreeArc(hub_out, d.in);
    net.addFreeArc(hub_out, e.in);
    ASSERT_TRUE(net.solve());
    EXPECT_EQ(net.multiplicity(a.arc), 1u);
    EXPECT_EQ(net.multiplicity(b.arc), 1u);
    EXPECT_EQ(net.multiplicity(x.arc), 2u);
    EXPECT_EQ(net.multiplicity(d.arc), 1u);
    EXPECT_EQ(net.multiplicity(e.arc), 1u);
//    X can not go to 3 without another vertex changing (cost >= 10) or new ends; going down costs at least its own 10.
    EXPECT_GT(net.increaseCost(x.arc, 1000), 9.9);
    EXPECT_GT(net.decreaseCost(x.arc, 1000), 9.9);
//    Multiplicity 1 is the lower bound.
    EXPECT_TRUE(std::isinf(net.decreaseCost(a.arc, 1000)));
}

//    Two parallel vertices with coverage 1.5 copies each between flanks of multiplicity 3, on a circle (so no genome
//    ends are needed): (1, 2) and (2, 1) cost the same, so the change is free.
TEST(MultiplicityNetwork, AmbiguousSplitIsFree) {
    MultiplicityNetwork net;
    double end = 50;
    LongArc left = AddLong(net, 3, end), right = AddLong(net, 3, end);
    LongArc p = AddLong(net, 1.5, end), q = AddLong(net, 1.5, end);
    net.addFreeArc(left.out, p.in);
    net.addFreeArc(left.out, q.in);
    net.addFreeArc(p.out, right.in);
    net.addFreeArc(q.out, right.in);
    net.addFreeArc(right.out, left.in);
    ASSERT_TRUE(net.solve());
    EXPECT_EQ(net.multiplicity(left.arc), 3u);
    EXPECT_EQ(net.multiplicity(right.arc), 3u);
    EXPECT_EQ(net.multiplicity(p.arc) + net.multiplicity(q.arc), 3u);
    size_t up = net.multiplicity(p.arc) == 1 ? p.arc : q.arc;
    size_t down = up == p.arc ? q.arc : p.arc;
    EXPECT_NEAR(net.increaseCost(up, 1000), 0, 1e-6);
    EXPECT_NEAR(net.decreaseCost(down, 1000), 0, 1e-6);
}

//    A vertex alone with coverage of 2 copies: it is a circle-free piece, so multiplicity 2 needs two extra ends.
//    Cheap ends: 2 wins. Expensive ends: 1 wins, and the cap hides the increase cost.
TEST(MultiplicityNetwork, EndsLimitMultiplicity) {
    {
        MultiplicityNetwork net;
        LongArc a = AddLong(net, 2, 1);
        ASSERT_TRUE(net.solve());
        EXPECT_EQ(net.multiplicity(a.arc), 2u);
        EXPECT_NEAR(net.decreaseCost(a.arc, 1000), 10 - 2, 1e-6);
    }
    {
        MultiplicityNetwork net;
        LongArc a = AddLong(net, 2, 100);
        ASSERT_TRUE(net.solve());
        EXPECT_EQ(net.multiplicity(a.arc), 1u);
        EXPECT_NEAR(net.increaseCost(a.arc, 1000), 200 - 10, 1e-6);
        EXPECT_TRUE(std::isinf(net.increaseCost(a.arc, 100)));
    }
}

//    Fixed boundary flow forces the multiplicity of the vertex between two fixed ports.
TEST(MultiplicityNetwork, FixedBoundaryFlow) {
    MultiplicityNetwork net;
    size_t source = net.addNode(), sink = net.addNode();
    net.addFixedStart(source, 4);
    net.addFixedEnd(sink, 4);
    net.addEnd(source, 30);
    net.addStart(sink, 30);
    LongArc a = AddLong(net, 3.4, 30);
    net.addFreeArc(source, a.in);
    net.addFreeArc(a.out, sink);
    ASSERT_TRUE(net.solve());
    EXPECT_EQ(net.multiplicity(a.arc), 4u);
}

//    Non-convex costs are convexified: a cost dip at 3 between higher costs still yields a consistent answer.
TEST(MultiplicityNetwork, NonConvexCostsAreConvexified) {
    MultiplicityNetwork net;
    size_t in = net.addNode(), out = net.addNode();
    size_t arc = net.addLongArc(in, out, {0, 10, 1, 20});
    net.addStart(in, 0);
    net.addEnd(out, 0);
    ASSERT_TRUE(net.solve());
    size_t m = net.multiplicity(arc);
    EXPECT_TRUE(m == 1 || m == 3);
    EXPECT_NEAR(net.totalCost(), 0, 1e-6);
}
