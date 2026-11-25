#pragma once

#include "assembly_graph.hpp"
#include "dll_path_storage.hpp"
#include "visualization.hpp"
#include <experimental/filesystem>
#include <functional>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ag {

//    Coverage of a long vertex given its multiplicity m (number of times the genome passes through it).
//
//    The observed value is the mean coverage x of a vertex of length l. It is treated as n = x * l / R independent read
//    samples, where R is the effective sample length (about the read length: neighbouring positions share reads). n
//    follows a negative binomial distribution with mean m * c * l / R and dispersion 1 / sigma^2, so in coverage units
//        E[x] = m * c,   Var[x] = m * c * R / l + (sigma * m * c)^2.
//    The first term is counting noise and vanishes for long vertices; the second is multiplicative bias (GC content,
//    mappability, library) that does not average out with length. prior[m - 1] is the fraction of long-vertex length
//    with multiplicity m.
    class CoverageModel {
    public:
        double unit_coverage = 1;   // c: coverage of one copy
        double sample_length = 1;   // R
        double dispersion = 0.1;    // sigma
        std::vector<double> prior;  // prior[m - 1] for m = 1..prior.size()

        double logDensity(double coverage, double length, size_t multiplicity) const;
        double logPrior(size_t multiplicity) const;
//        -log(density * prior): cost of explaining the vertex with this multiplicity.
        double cost(double coverage, double length, size_t multiplicity) const;
//        log of the mixture density over all multiplicities in the prior.
        double logLikelihood(double coverage, double length) const;
        size_t mostLikely(double coverage, double length) const;
        std::string str() const;

//        Fits the model to (coverage, length) of long vertices (one per rc pair) by EM, length weighted.
//        Assumes that most of the genome is present in one copy: c starts at the mode of the read-mass weighted
//        coverage histogram and EM only refines it (multiplying all multiplicities by j and dividing c by j would
//        explain the data as well, so the data alone can not fix c). Inbred genomes, mostly one copy per haplotype
//        pair, need a separate model. unit_hint > 0 fixes c instead.
        static CoverageModel Fit(logging::Logger &logger, const std::vector<std::pair<double, double>> &data,
                                 size_t max_multiplicity, double unit_hint);
    };

    struct MultiplicityParams {
//        Vertices at least this long get a multiplicity. Must be at least 3-4 typical read lengths: shorter vertices
//        include single uncorrected reads, and their coverage is too noisy to estimate multiplicity.
        size_t min_long = 40000;
//        Each flow problem gets at most about this many long vertices (longest first). Larger regions are solved for
//        their longest vertices first, then split at them and the parts are solved separately (see Estimate).
        size_t max_long_per_problem = 2000;
        size_t max_multiplicity = 100;
//        Cost of one genome end (start) at the side of a long vertex of length l that connects to no other long vertex:
//        -log(rate * l), where rate is the frequency per bp of obvious breaks (dead out-sides of long vertices with
//        about one copy of coverage, per bp of their total length; P capped at 1/2). Anywhere else (a coverage break
//        inside a component) the end costs break_cost_factor times as much.
        double break_cost_factor = 3;
//        Increase/decrease costs above ambiguity_cap_factor * break cost are reported as infinite.
        double ambiguity_cap_factor = 20;
//        Unit (one copy) coverage; 0 to estimate.
        double unit_coverage = 0;
//        Short dead ends of the graph with at least this coverage may be real sequence ending in a coverage break: the
//        genome may end (start) there, at the end cost - log P(real) (see DeadEndEstimate). Lower ones are errors.
        double min_dead_end_coverage = 2;
//        Prior probability that such a dead end is real.
        double dead_end_prior = 0.5;
    };

//    A short dead end of the graph (a vertex without outgoing edges, or without incoming ones for a start) where the
//    genome may end (start). Judged by the unbranching stretch of short vertices leading to it: its spelled length and
//    its length-weighted coverage, compared under the one-copy coverage model (real sequence whose reads stop at a
//    coverage break) and under the error coverage model (fitted to vertices below a quarter of a copy).
    struct DeadEndEstimate {
        bool forward;           // the genome ends here; otherwise it starts here
        size_t length;          // of the stretch
        double coverage;        // of the stretch
        double p_real;
        double cost;            // of one genome end (start) here
        size_t flow;            // genome ends (starts) placed here by the solution
        size_t level;
    };

    struct MultiplicityEstimate {
        size_t multiplicity = 0;
//        Minimal increase of total cost (nats) of a solution where this vertex has one copy more / less, with everything
//        else free to adjust. Infinite if impossible or above the cap.
        double increase_cost = 0;
        double decrease_cost = 0;
//        coverage / unit coverage; -1 if the vertex has no coverage information
        double expected = 0;
//        Recursion level of the flow problem that fixed it and that problem's length threshold.
        size_t level = 0;
        size_t threshold = 0;
    };

//    Multiplicities of long vertices by min-cost circulation.
//
//    The graph is cut at long vertices (each becomes an arc from its in-half to its out-half carrying its multiplicity,
//    with cost from CoverageModel). The short vertices between them are condensed into strongly connected components
//    and contracted, which keeps exactly the information which entry can reach which exit. Genome ends cost the end cost
//    where the contracted graph has a dead end and the larger break cost anywhere else (see MultiplicityParams). Short
//    dead ends with enough coverage stay in the network: the genome may end (start) there too, at a cost growing as
//    they look less like real sequence (see DeadEndEstimate).
//
//    Regions with more than max_long_per_problem long vertices are processed iteratively: the longest ones are solved
//    first, then each component of the remaining shorter vertices is solved separately (in parallel), with the flow
//    entering and leaving it through the already solved vertices fixed by the parent solution.
//    Weakly connected components that are not their own reverse complement are solved for one strand and mirrored.
    class LongVertexMultiplicities {
    public:
        static LongVertexMultiplicities Estimate(logging::Logger &logger, size_t threads, const AssemblyGraph &graph,
                                                 const MultiplicityParams &params);

        const MultiplicityEstimate *get(const Vertex &vertex) const;
        const CoverageModel &getModel() const {return model;}
//        End and break costs after a long vertex of mean length (see MultiplicityParams::break_cost_factor).
        double getEndCost() const {return end_cost;}
        double getBreakCost() const {return break_cost;}
        double getBreakRate() const {return break_rate;}
//        How a reference passes through a vertex (see DLLAlignmentStorage::vertexTraversals): traversals away from the
//        ends of reference paths (the true multiplicity), traversals touching an end (the copy may be cut there), and
//        the smallest fraction of the vertex spelled by one traversal (graph and reference may differ).
        struct ReferenceCount {
            size_t traversals = 0;
            size_t at_ends = 0;
            double min_covered = 1;
        };
        typedef std::function<ReferenceCount(const Vertex &)> ReferenceCounter;
//        Counts from the contig paths tracked by paths, collected once by this call.
        static ReferenceCounter FromPaths(DLLAlignmentStorage &paths);

//        One tab-separated line per vertex with an estimate (both strands); with reference columns if one is given.
        void printReport(const std::experimental::filesystem::path &path, const ReferenceCounter &reference = {}) const;
//        Logs estimated against true multiplicity of canonical vertices: counts and total length of correct, too high
//        and too low estimates, how many wrong ones were flagged (the change towards the truth costs less than one
//        genome end), and the table of (estimated, true) pairs. Vertices touching a reference end are listed apart.
        void compareWithReference(logging::Logger &logger, const ReferenceCounter &reference) const;
//        One tab-separated line per dead end where the genome may end or start (see DeadEndEstimate), as judged by the
//        last flow problem that contained it; with the reference traversals of the dead-end vertex if a reference is
//        given.
        void printDeadEnds(const std::experimental::filesystem::path &path, const ReferenceCounter &reference = {}) const;
//        "Mult:m (+increase/-decrease)" for vertices with an estimate; "End/Start p=P(real) used=flow" for dead ends.
        VertexInfo vertexInfo() const;

    private:
        CoverageModel model;
//        Coverage of erroneous vertices, for judging dead ends.
        CoverageModel error_model;
        double break_rate = 0;
        double end_cost = 0;
        double break_cost = 0;
        std::unordered_map<const Vertex *, MultiplicityEstimate> estimates;
//        (vertex, forward) -> estimate
        std::map<std::pair<const Vertex *, bool>, DeadEndEstimate> dead_ends;
    };
}
