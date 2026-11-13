#pragma once

#include "supregraph.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"
#include "common/logging.hpp"
#include <optional>

namespace spg {
//    A simple loop: a fork core (in- and out-degree 2) with a loop core -l_out-> body -l_in-> core through a body
//    with in- and out-degree 1, entered from outside by e_in and left by e_out. The walk through it is
//    e_in (l_out l_in)^n e_out for some number n >= 1 of body copies.
    struct SimpleLoop {
        ag::VertexId core;
        ag::VertexId body;
        ag::EdgeId e_in;
        ag::EdgeId l_out;
        ag::EdgeId l_in;
        ag::EdgeId e_out;
    };

//    The simple loop at core, if core has one. No vertex of the loop, nor its entrance start or exit finish, may coincide
//    with another or with another's rc, and the entrance and exit must lead to different vertices.
    std::optional<SimpleLoop> DetectSimpleLoop(ag::Vertex &core);

//    Number of body copies the coverage suggests, as a real number: the entrance start and exit finish give the
//    coverage of one copy; the body is traversed n times and the core n + 1 times, so both give an estimate of n, and
//    they are averaged with their lengths as weights. Negative if the neighbours have no coverage.
    double EstimateLoopCopies(const SimpleLoop &loop);

//    Unrolls every simple loop of the graph into the number of body copies estimated from coverage (loops estimated at
//    less than one copy or more than max_copies are left alone). Each copy but the last is one resolveVertex that keeps
//    the loop pair (in a fresh core and body), the last one resolves the loop pair away. Plans never break a read: a
//    loop is left (possibly partly unrolled) when some read bypasses it, when the read support of the core is unknown,
//    or when the loop does not return to the simple shape after a step; reads going around more often than the
//    coverage says get the extra copies. Run on a merged graph (MergeAllSPG), single-threaded; leaves unmerged paths
//    behind. Returns the number of loops fully unrolled.
    size_t ResolveLoops(logging::Logger &logger, Supregraph &spg, ag::SuffixTracker &suffixes, size_t max_copies);
}
