#include "long_vertex_multiplicity.hpp"
#include "multiplicity_network.hpp"
#include "common/verify.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_set>

using namespace ag;

namespace {
    constexpr double PRIOR_FLOOR = 1e-4;
    constexpr double MIN_DISPERSION = 0.01;
    constexpr double MAX_SAMPLE_LENGTH = 1e6;
    constexpr size_t NO_TAG = size_t(-1);

    double LogSumExp(const std::vector<double> &values) {
        double mx = *std::max_element(values.begin(), values.end());
        double sum = 0;
        for (double v : values)
            sum += std::exp(v - mx);
        return mx + std::log(sum);
    }

    size_t ModelSize(double max_coverage, double unit, size_t max_multiplicity) {
        double size = std::ceil(max_coverage / unit) + 2;
        return std::max<size_t>(1, std::min<size_t>(max_multiplicity, size_t(size)));
    }

//    Mode of the coverage histogram in log scale, weighted by read mass (length * coverage): long vertices made of single
//    uncorrected reads can outweigh real ones by length, but not by the reads they contain.
    double CoverageMode(const std::vector<std::pair<double, double>> &data) {
        const double bin = 0.02;
        const int smoothing = 10;
        double lo = std::numeric_limits<double>::infinity();
        double hi = -lo;
        for (const auto &[cov, len] : data) {
            lo = std::min(lo, std::log(cov));
            hi = std::max(hi, std::log(cov));
        }
        size_t bins = size_t((hi - lo) / bin) + 1;
        std::vector<double> hist(bins, 0);
        for (const auto &[cov, len] : data)
            hist[size_t((std::log(cov) - lo) / bin)] += len * cov;
        size_t best = 0;
        double best_value = -1;
        for (size_t i = 0; i < bins; i++) {
            double value = 0;
            for (int d = -smoothing; d <= smoothing; d++) {
                long long j = (long long)i + d;
                if (j >= 0 && j < (long long)bins)
                    value += hist[j] * double(smoothing + 1 - std::abs(d));
            }
            if (value > best_value) {
                best_value = value;
                best = i;
            }
        }
        return std::exp(lo + (double(best) + 0.5) * bin);
    }

//    EM from the initial unit coverage; returns the model and its length-weighted mean log-likelihood.
    std::pair<CoverageModel, double> FitEM(const std::vector<std::pair<double, double>> &data, double unit,
                                           bool fix_unit, size_t max_multiplicity) {
        double max_coverage = 0;
        double total_length = 0;
        for (const auto &[cov, len] : data) {
            max_coverage = std::max(max_coverage, cov);
            total_length += len;
        }
        CoverageModel model;
        model.unit_coverage = unit;
        model.sample_length = 10000;
        model.dispersion = 0.1;
        double score = -std::numeric_limits<double>::infinity();
        for (size_t iteration = 0; iteration < 200; iteration++) {
            size_t size = ModelSize(max_coverage, model.unit_coverage, max_multiplicity);
            if (model.prior.size() != size)
                model.prior.assign(size, 1.0 / double(size));
//            E-step
            std::vector<std::vector<double>> weights;
            double new_score = 0;
            for (const auto &[cov, len] : data) {
                std::vector<double> logw;
                for (size_t m = 1; m <= size; m++)
                    logw.push_back(model.logPrior(m) + model.logDensity(cov, len, m));
                double norm = LogSumExp(logw);
                new_score += len * norm;
                for (double &w : logw)
                    w = std::exp(w - norm);
                weights.push_back(std::move(logw));
            }
            new_score /= total_length;
//            M-step: prior and unit coverage
            std::vector<double> prior(size, 0);
            double cov_sum = 0;
            double mult_sum = 0;
            for (size_t i = 0; i < data.size(); i++) {
                const auto &[cov, len] = data[i];
                for (size_t m = 1; m <= size; m++) {
                    prior[m - 1] += len * weights[i][m - 1];
                    cov_sum += len * weights[i][m - 1] * cov;
                    mult_sum += len * weights[i][m - 1] * double(m);
                }
            }
            for (double &p : prior)
                p = std::max(p / total_length, PRIOR_FLOOR);
            double prior_sum = std::accumulate(prior.begin(), prior.end(), 0.0);
            for (double &p : prior)
                p /= prior_sum;
            model.prior = prior;
            if (!fix_unit)
                model.unit_coverage = cov_sum / mult_sum;
//            M-step: variance terms by weighted least squares of squared residuals on (m c / l, (m c)^2),
//            reweighted by the previous predicted variance (squared residuals have variance ~ prediction^2).
            double c = model.unit_coverage;
            double s_aa = 0, s_ab = 0, s_bb = 0, s_ay = 0, s_by = 0;
            for (size_t i = 0; i < data.size(); i++) {
                const auto &[cov, len] = data[i];
                for (size_t m = 1; m <= size; m++) {
                    double w = weights[i][m - 1];
                    if (w < 1e-6)
                        continue;
                    double mean = double(m) * c;
                    double a = mean / len;
                    double b = mean * mean;
                    double y = (cov - mean) * (cov - mean);
                    double pred = model.sample_length * a + model.dispersion * model.dispersion * b;
                    double u = len * w / (pred * pred);
                    s_aa += u * a * a;
                    s_ab += u * a * b;
                    s_bb += u * b * b;
                    s_ay += u * a * y;
                    s_by += u * b * y;
                }
            }
            double det = s_aa * s_bb - s_ab * s_ab;
            double r = det > 0 ? (s_ay * s_bb - s_by * s_ab) / det : -1;
            double s = det > 0 ? (s_by * s_aa - s_ay * s_ab) / det : -1;
            if (r < 1 || s < MIN_DISPERSION * MIN_DISPERSION) {
//                One of the terms is not supported by the data: clamp it and refit the other one alone.
                double r_only = s_aa > 0 ? (s_ay - MIN_DISPERSION * MIN_DISPERSION * s_ab) / s_aa : 1;
                double s_only = s_bb > 0 ? (s_by - s_ab) / s_bb : MIN_DISPERSION * MIN_DISPERSION;
                if (s_only >= MIN_DISPERSION * MIN_DISPERSION) {
                    r = 1;
                    s = s_only;
                } else {
                    r = std::max(1.0, r_only);
                    s = MIN_DISPERSION * MIN_DISPERSION;
                }
            }
            model.sample_length = std::min(r, MAX_SAMPLE_LENGTH);
            model.dispersion = std::sqrt(s);
            bool converged = std::abs(new_score - score) < 1e-7;
            score = new_score;
            if (converged)
                break;
        }
        return {model, score};
    }

    struct Port {
        const Vertex *vertex;
        bool out_side;  // flow leaves vertex through its end into the region; otherwise it enters vertex from the region
        size_t flow;
    };

    struct Region {
        std::vector<const Vertex *> vertices;
        std::vector<Port> ports;
        size_t level = 0;
        bool mirrored = false;
    };

    struct RegionResult {
        std::vector<std::pair<const Vertex *, MultiplicityEstimate>> estimates;
        std::vector<Region> children;
        size_t threshold = 0;
        size_t long_count = 0;
        size_t nodes = 0;
        size_t end_flow = 0;
        size_t break_flow = 0;
        std::vector<std::pair<std::pair<const Vertex *, bool>, DeadEndEstimate>> dead_ends;  // ((vertex, forward), ...)
        size_t dead_end_flow = 0;
    };

//    Sources (flow enters the short part here) are out-halves of long vertices and port_out; sinks are in-halves and
//    port_in. scc nodes are condensed strongly connected components of short vertices. dead_end (dead_start) is where
//    the genome may end (start) at a short dead end of the graph that may be real sequence with a coverage break.
    enum class NodeKind { long_in, long_out, port_in, port_out, scc, dead_end, dead_start };

    bool IsPortSource(NodeKind kind) {return kind == NodeKind::long_out || kind == NodeKind::port_out;}
    bool IsPortSink(NodeKind kind) {return kind == NodeKind::long_in || kind == NodeKind::port_in;}
    bool IsSource(NodeKind kind) {return IsPortSource(kind) || kind == NodeKind::dead_start;}
    bool IsSink(NodeKind kind) {return IsPortSink(kind) || kind == NodeKind::dead_end;}

//    Network of the short part between ports. Every arc carries the tag of the component of short vertices it passes
//    through (NO_TAG for an edge between two ports). trim() and contract() keep exactly which sources reach which sinks.
    class Gadget {
    public:
        size_t addNode(NodeKind kind, size_t tag) {
            kinds.push_back(kind);
            tags.push_back(tag);
            alive.push_back(true);
            out.emplace_back();
            in.emplace_back();
            return kinds.size() - 1;
        }

        void addArc(size_t from, size_t to, size_t tag) {
            if (from == to)
                return;
            out[from].emplace(to, tag);
            in[to].emplace(from, tag);
        }

//        Removes scc nodes that no source reaches or that reach no sink.
        void trim() {
            std::vector<bool> forward = reach(true);
            std::vector<bool> backward = reach(false);
            for (size_t x = 0; x < kinds.size(); x++)
                if (kinds[x] == NodeKind::scc && !(forward[x] && backward[x]))
                    remove(x);
        }

//        Repeatedly merges an scc node with a single in-neighbour (or out-neighbour) into it.
        void contract() {
            std::vector<size_t> work;
            for (size_t x = 0; x < kinds.size(); x++)
                if (alive[x] && kinds[x] == NodeKind::scc)
                    work.push_back(x);
            while (!work.empty()) {
                size_t x = work.back();
                work.pop_back();
                if (!alive[x])
                    continue;
                if (in[x].size() == 1) {
                    size_t p = in[x].begin()->first;
                    std::vector<std::pair<size_t, size_t>> next(out[x].begin(), out[x].end());
                    remove(x);
                    for (const auto &[y, tag] : next) {
                        addArc(p, y, tags[x]);
                        push(work, y);
                    }
                    push(work, p);
                } else if (out[x].size() == 1) {
                    size_t q = out[x].begin()->first;
                    std::vector<std::pair<size_t, size_t>> prev(in[x].begin(), in[x].end());
                    remove(x);
                    for (const auto &[y, tag] : prev) {
                        addArc(y, q, tags[x]);
                        push(work, y);
                    }
                    push(work, q);
                }
            }
        }

        std::vector<NodeKind> kinds;
        std::vector<size_t> tags;
        std::vector<bool> alive;
        std::vector<std::set<std::pair<size_t, size_t>>> out;  // (neighbour, tag)
        std::vector<std::set<std::pair<size_t, size_t>>> in;

//        Nodes reachable from a source (forward) or reaching a sink (backward); with ports_only, dead_end and dead_start
//        nodes do not count as sinks and sources.
        std::vector<bool> reach(bool forward, bool ports_only = false) const {
            std::vector<bool> res(kinds.size(), false);
            std::vector<size_t> stack;
            for (size_t x = 0; x < kinds.size(); x++) {
                bool seed = forward ? (ports_only ? IsPortSource(kinds[x]) : IsSource(kinds[x]))
                                    : (ports_only ? IsPortSink(kinds[x]) : IsSink(kinds[x]));
                if (alive[x] && seed) {
                    res[x] = true;
                    stack.push_back(x);
                }
            }
            while (!stack.empty()) {
                size_t x = stack.back();
                stack.pop_back();
                for (const auto &[y, tag] : forward ? out[x] : in[x])
                    if (!res[y]) {
                        res[y] = true;
                        stack.push_back(y);
                    }
            }
            return res;
        }

    private:
        void push(std::vector<size_t> &work, size_t x) const {
            if (alive[x] && kinds[x] == NodeKind::scc)
                work.push_back(x);
        }

        void remove(size_t x) {
            for (const auto &[y, tag] : out[x])
                in[y].erase({x, tag});
            for (const auto &[y, tag] : in[x])
                out[y].erase({x, tag});
            out[x].clear();
            in[x].clear();
            alive[x] = false;
        }
    };

//    Strongly connected components (Tarjan, iterative) of a graph on 0..n-1.
    std::vector<size_t> StronglyConnected(const std::vector<std::vector<size_t>> &succ, size_t &count) {
        size_t n = succ.size();
        const size_t UNSET = size_t(-1);
        std::vector<size_t> index(n, UNSET), low(n, 0), comp(n, UNSET);
        std::vector<bool> on_stack(n, false);
        std::vector<size_t> stack;
        std::vector<std::pair<size_t, size_t>> calls;  // (vertex, next successor position)
        size_t counter = 0;
        count = 0;
        for (size_t start = 0; start < n; start++) {
            if (index[start] != UNSET)
                continue;
            calls.emplace_back(start, 0);
            index[start] = low[start] = counter++;
            stack.push_back(start);
            on_stack[start] = true;
            while (!calls.empty()) {
                size_t v = calls.back().first;
                size_t &pos = calls.back().second;
                if (pos < succ[v].size()) {
                    size_t w = succ[v][pos++];
                    if (index[w] == UNSET) {
                        index[w] = low[w] = counter++;
                        stack.push_back(w);
                        on_stack[w] = true;
                        calls.emplace_back(w, 0);
                    } else if (on_stack[w]) {
                        low[v] = std::min(low[v], index[w]);
                    }
                    continue;
                }
                if (low[v] == index[v]) {
                    while (true) {
                        size_t w = stack.back();
                        stack.pop_back();
                        on_stack[w] = false;
                        comp[w] = count;
                        if (w == v)
                            break;
                    }
                    count++;
                }
                calls.pop_back();
                if (!calls.empty()) {
                    size_t parent = calls.back().first;
                    low[parent] = std::min(low[parent], low[v]);
                }
            }
        }
        return comp;
    }

    struct BreakCount {
        size_t dead_sides = 0;
        size_t vertices = 0;
        double length = 0;
    };

//    Obvious genome breaks: long vertices with about one copy of coverage (0.5-1.5) whose out-side reaches no other long
//    vertex through short vertices. Counted over both strands, so every dead out-side is one end of a genome copy.
    BreakCount CountObviousBreaks(const AssemblyGraph &graph, size_t min_long, double unit_coverage) {
//        Short vertices from which a long vertex can be reached: backwards from every long vertex.
        std::unordered_set<const Vertex *> reach_long;
        std::vector<const Vertex *> stack;
        for (const Vertex &w : graph.vertices())
            if (w.size() >= min_long)
                stack.push_back(&w);
        while (!stack.empty()) {
            const Vertex *x = stack.back();
            stack.pop_back();
            for (const Edge &e : x->incoming()) {
                const Vertex *u = &e.getStart();
                if (u->size() < min_long && reach_long.insert(u).second)
                    stack.push_back(u);
            }
        }
        BreakCount res;
        for (const Vertex &v : graph.vertices()) {
            if (v.size() < min_long || !v.hasCoverageInfo())
                continue;
            double copies = v.getSPGCoverage() / unit_coverage;
            if (copies < 0.5 || copies > 1.5)
                continue;
            res.vertices++;
            res.length += double(v.size());
            bool alive = false;
            for (const Edge &e : v) {
                const Vertex &x = e.getFinish();
                if (x.size() >= min_long || reach_long.count(&x))
                    alive = true;
            }
            if (!alive)
                res.dead_sides++;
        }
        return res;
    }

//    -log P(the genome ends right after a stretch of this length), for a break rate per bp; P is capped at 1/2.
    double EndCost(double break_rate, double length) {
        return -std::log(std::min(0.5, break_rate * length));
    }

    struct Context {
        const CoverageModel &model;
        const CoverageModel &error_model;
        const MultiplicityParams &params;
        double break_rate;
//        End and break costs after a long vertex of mean length: for the ambiguity cap.
        double end_cost;
        double break_cost;

//        Cost of one genome end (start) at the side of a vertex of this length: at a topological dead end or elsewhere.
        double endCost(double length) const {return EndCost(break_rate, length);}
        double breakCost(double length) const {return params.break_cost_factor * endCost(length);}

//        P(real) and the cost of one genome end at a dead end whose stretch has this length and coverage.
        std::pair<double, double> deadEndCost(double length, double coverage) const {
            double real = std::log(params.dead_end_prior) + model.logDensity(coverage, length, 1);
            double error = std::log(1 - params.dead_end_prior) + error_model.logDensity(coverage, length, 1);
            double log_p_real = real - LogSumExp({real, error});
            return {std::exp(log_p_real), endCost(length) - log_p_real};
        }

        std::vector<double> longCosts(const Vertex &v) const {
            std::vector<double> costs;
            if (!v.hasCoverageInfo()) {
                for (size_t m = 1; m <= params.max_multiplicity; m++)
                    costs.push_back(-model.logPrior(m));
                return costs;
            }
            double cov = v.getSPGCoverage();
            double expected = cov / model.unit_coverage;
            size_t size = std::min(params.max_multiplicity, size_t(std::ceil(1.5 * expected)) + 3);
            for (size_t m = 1; m <= size; m++)
                costs.push_back(model.cost(cov, double(v.size()), m));
            return costs;
        }
    };

//    The unbranching stretch of short vertices leading into the dead end v (forward) or out of it (backward), up to the
//    last branching: its spelled length and the length-weighted coverage of its vertices with coverage information.
    std::pair<size_t, double> DeadEndStretch(const Vertex &v, bool forward,
                                             const std::function<bool(const Vertex &)> &is_short) {
        size_t length = v.size();
        double coverage_sum = 0;
        double weight = 0;
        std::unordered_set<const Vertex *> visited;
        auto add = [&](const Vertex &x) {
            visited.insert(&x);
            if (x.hasCoverageInfo()) {
                coverage_sum += x.getSPGCoverage() * double(x.size());
                weight += double(x.size());
            }
        };
        add(v);
        const Vertex *cur = &v;
        while (forward ? cur->inDeg() == 1 : cur->outDeg() == 1) {
            const Edge &e = forward ? cur->incFront() : cur->front();
            const Vertex &next = forward ? e.getStart() : e.getFinish();
            bool unbranching = forward ? next.outDeg() == 1 : next.inDeg() == 1;
            if (!is_short(next) || !unbranching || visited.count(&next))
                break;
//            Either way the edge spells its full sequence, which overlaps cur's label by cur's size.
            length += e.fullSize() - cur->size();
            add(next);
            cur = &next;
        }
        return {length, weight > 0 ? coverage_sum / weight : 0};
    }

    RegionResult ProcessRegion(const Region &region, const Context &ctx) {
        const MultiplicityParams &params = ctx.params;
        RegionResult res;
        std::vector<size_t> lengths;
        for (const Vertex *v : region.vertices)
            if (v->size() >= params.min_long)
                lengths.push_back(v->size());
        VERIFY(!lengths.empty());
        std::sort(lengths.rbegin(), lengths.rend());
        size_t threshold = params.min_long;
        if (lengths.size() > params.max_long_per_problem)
            threshold = std::max(params.min_long, lengths[params.max_long_per_problem - 1]);
        bool last_level = threshold == params.min_long;
        res.threshold = threshold;

        std::unordered_map<const Vertex *, size_t> local;  // region vertex -> index in region.vertices
        for (size_t i = 0; i < region.vertices.size(); i++)
            local[region.vertices[i]] = i;
        auto isLong = [threshold](const Vertex *v) {return v->size() >= threshold;};

//        Short part: strongly connected components and weakly connected components (tags).
        std::vector<size_t> short_index(region.vertices.size(), size_t(-1));
        std::vector<const Vertex *> shorts;
        for (size_t i = 0; i < region.vertices.size(); i++)
            if (!isLong(region.vertices[i])) {
                short_index[i] = shorts.size();
                shorts.push_back(region.vertices[i]);
            }
        auto shortId = [&](const Vertex *v) -> size_t {
            auto it = local.find(v);
            return it == local.end() ? size_t(-1) : short_index[it->second];
        };
        std::vector<std::vector<size_t>> succ(shorts.size());
        std::vector<size_t> parent(shorts.size());
        std::iota(parent.begin(), parent.end(), 0);
        std::function<size_t(size_t)> find = [&](size_t x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };
        for (size_t i = 0; i < shorts.size(); i++)
            for (const Edge &e : *shorts[i]) {
                size_t j = shortId(&e.getFinish());
                if (j != size_t(-1)) {
                    succ[i].push_back(j);
                    parent[find(i)] = find(j);
                }
            }
        size_t scc_count = 0;
        std::vector<size_t> scc = StronglyConnected(succ, scc_count);
        std::vector<size_t> tag(shorts.size());
        for (size_t i = 0; i < shorts.size(); i++)
            tag[i] = find(i);

//        Gadget nodes
        Gadget gadget;
        std::unordered_map<const Vertex *, std::pair<size_t, size_t>> long_nodes;  // (in, out)
        std::vector<const Vertex *> longs;
        for (const Vertex *v : region.vertices)
            if (isLong(v)) {
                size_t in = gadget.addNode(NodeKind::long_in, NO_TAG);
                size_t out = gadget.addNode(NodeKind::long_out, NO_TAG);
                long_nodes[v] = {in, out};
                longs.push_back(v);
            }
        std::map<std::pair<const Vertex *, bool>, size_t> port_nodes;
        for (const Port &port : region.ports)
            port_nodes[{port.vertex, port.out_side}] =
                    gadget.addNode(port.out_side ? NodeKind::port_out : NodeKind::port_in, NO_TAG);
        std::vector<size_t> scc_nodes(scc_count, size_t(-1));
        for (size_t i = 0; i < shorts.size(); i++)
            if (scc_nodes[scc[i]] == size_t(-1))
                scc_nodes[scc[i]] = gadget.addNode(NodeKind::scc, tag[i]);
//        Node where flow leaves (out) / enters (in) a vertex; region vertices and ports
        auto outNode = [&](const Vertex *v) -> size_t {
            size_t s = shortId(v);
            if (s != size_t(-1))
                return scc_nodes[scc[s]];
            auto it = long_nodes.find(v);
            if (it != long_nodes.end())
                return it->second.second;
            auto port = port_nodes.find({v, true});
            VERIFY_MSG(port != port_nodes.end(), "Region edge leads to a vertex that is neither in the region nor a port");
            return port->second;
        };
        auto inNode = [&](const Vertex *v) -> size_t {
            size_t s = shortId(v);
            if (s != size_t(-1))
                return scc_nodes[scc[s]];
            auto it = long_nodes.find(v);
            if (it != long_nodes.end())
                return it->second.first;
            auto port = port_nodes.find({v, false});
            VERIFY_MSG(port != port_nodes.end(), "Region edge leads to a vertex that is neither in the region nor a port");
            return port->second;
        };
        auto edgeTag = [&](const Vertex *u, const Vertex *x) -> size_t {
            size_t su = shortId(u);
            if (su != size_t(-1))
                return tag[su];
            size_t sx = shortId(x);
            return sx != size_t(-1) ? tag[sx] : NO_TAG;
        };
        for (const Vertex *u : region.vertices) {
            for (const Edge &e : *u) {
                const Vertex *x = &e.getFinish();
                gadget.addArc(outNode(u), inNode(x), edgeTag(u, x));
            }
            for (const Edge &e : u->incoming()) {
                const Vertex *x = &e.getStart();
                if (local.find(x) == local.end())
                    gadget.addArc(outNode(x), inNode(u), edgeTag(x, u));
            }
        }
//        Short dead ends that may be real: the genome may end (start) there.
        struct DeadEnd {
            const Vertex *vertex;
            bool forward;
            size_t node;
            size_t length;
            double coverage;
        };
        std::vector<DeadEnd> dead_ends;
        auto isShortHere = [&](const Vertex &v) {return shortId(&v) != size_t(-1);};
        for (size_t i = 0; i < shorts.size(); i++) {
            const Vertex &s = *shorts[i];
            if (!s.hasCoverageInfo() || s.getSPGCoverage() < params.min_dead_end_coverage)
                continue;
            for (bool forward : {true, false}) {
                if (forward ? s.outDeg() != 0 : s.inDeg() != 0)
                    continue;
                size_t node = gadget.addNode(forward ? NodeKind::dead_end : NodeKind::dead_start, tag[i]);
                if (forward)
                    gadget.addArc(scc_nodes[scc[i]], node, tag[i]);
                else
                    gadget.addArc(node, scc_nodes[scc[i]], tag[i]);
                auto [length, coverage] = DeadEndStretch(s, forward, isShortHere);
                dead_ends.push_back({&s, forward, node, length, coverage});
            }
        }
//        A port is at a topological dead end if it connects to no other port (dead ends of the graph do not count).
        std::vector<bool> from_port_source = gadget.reach(true, true);
        std::vector<bool> to_port_sink = gadget.reach(false, true);
        gadget.trim();
        gadget.contract();

//        Flow network
        MultiplicityNetwork network;
        std::vector<size_t> net_node(gadget.kinds.size(), size_t(-1));
        for (size_t x = 0; x < gadget.kinds.size(); x++)
            if (gadget.alive[x])
                net_node[x] = network.addNode();
        res.nodes = network.nodeCount();
        std::vector<size_t> long_arcs;
        for (const Vertex *v : longs) {
            auto [in, out] = long_nodes.at(v);
            long_arcs.push_back(network.addLongArc(net_node[in], net_node[out], ctx.longCosts(*v)));
        }
        struct TaggedArc {size_t arc; size_t from; size_t to; size_t tag;};
        std::vector<TaggedArc> free_arcs;
        for (size_t x = 0; x < gadget.kinds.size(); x++)
            for (const auto &[y, t] : gadget.out[x])
                free_arcs.push_back({network.addFreeArc(net_node[x], net_node[y]), x, y, t});
//        The vertex of every long half and port: end costs depend on its length.
        std::unordered_map<size_t, const Vertex *> node_vertex;
        for (const auto &[v, halves] : long_nodes) {
            node_vertex[halves.first] = v;
            node_vertex[halves.second] = v;
        }
        for (const auto &[key, node] : port_nodes)
            node_vertex[node] = key.first;
        std::vector<std::pair<size_t, bool>> end_arcs;  // (arc, at a dead end)
        for (size_t x = 0; x < gadget.kinds.size(); x++) {
            if (!gadget.alive[x] || !(IsPortSink(gadget.kinds[x]) || IsPortSource(gadget.kinds[x])))
                continue;
            double length = double(node_vertex.at(x)->size());
            if (IsPortSink(gadget.kinds[x])) {
                bool dead = !from_port_source[x];
                double cost = dead ? ctx.endCost(length) : ctx.breakCost(length);
                end_arcs.emplace_back(network.addStart(net_node[x], cost), dead);
            } else {
                bool dead = !to_port_sink[x];
                double cost = dead ? ctx.endCost(length) : ctx.breakCost(length);
                end_arcs.emplace_back(network.addEnd(net_node[x], cost), dead);
            }
        }
        std::vector<size_t> dead_end_arcs(dead_ends.size(), size_t(-1));
        std::vector<std::pair<double, double>> dead_end_costs;  // (P(real), cost)
        for (size_t i = 0; i < dead_ends.size(); i++) {
            const DeadEnd &d = dead_ends[i];
            dead_end_costs.push_back(ctx.deadEndCost(double(d.length), d.coverage));
            if (!gadget.alive[d.node])
                continue;
            double cost = dead_end_costs.back().second;
            dead_end_arcs[i] = d.forward ? network.addEnd(net_node[d.node], cost)
                                         : network.addStart(net_node[d.node], cost);
        }
        for (const Port &port : region.ports) {
            size_t node = net_node[port_nodes.at({port.vertex, port.out_side})];
            if (port.out_side)
                network.addFixedStart(node, port.flow);
            else
                network.addFixedEnd(node, port.flow);
        }
        bool solved = network.solve();
        VERIFY_MSG(solved, "Multiplicity flow problem has no solution");
        for (const auto &[arc, dead] : end_arcs)
            (dead ? res.end_flow : res.break_flow) += network.flow(arc);
        for (size_t i = 0; i < dead_ends.size(); i++) {
            const DeadEnd &d = dead_ends[i];
            size_t flow = dead_end_arcs[i] == size_t(-1) ? 0 : network.flow(dead_end_arcs[i]);
            res.dead_end_flow += flow;
            res.dead_ends.push_back({{d.vertex, d.forward}, DeadEndEstimate{d.forward, d.length, d.coverage,
                    dead_end_costs[i].first, dead_end_costs[i].second, flow, region.level}});
        }

        double cap = params.ambiguity_cap_factor * ctx.break_cost;
        res.long_count = longs.size();
        for (size_t i = 0; i < longs.size(); i++) {
            const Vertex &v = *longs[i];
            MultiplicityEstimate est;
            est.multiplicity = network.multiplicity(long_arcs[i]);
            est.increase_cost = network.increaseCost(long_arcs[i], cap);
            est.decrease_cost = network.decreaseCost(long_arcs[i], cap);
            est.expected = v.hasCoverageInfo() ? v.getSPGCoverage() / ctx.model.unit_coverage : -1;
            est.level = region.level;
            est.threshold = threshold;
            res.estimates.emplace_back(&v, est);
        }
        if (last_level)
            return res;

//        Children: components of short vertices that still contain vertices to estimate, with the flow the solution
//        sends through their boundary.
        std::map<std::pair<size_t, size_t>, size_t> tag_out_flow;  // (gadget node, tag) -> flow leaving into the tag
        std::map<std::pair<size_t, size_t>, size_t> tag_in_flow;
        for (const TaggedArc &arc : free_arcs) {
            if (arc.tag == NO_TAG)
                continue;
            size_t f = network.flow(arc.arc);
            tag_out_flow[{arc.from, arc.tag}] += f;
            tag_in_flow[{arc.to, arc.tag}] += f;
        }
        std::unordered_map<size_t, std::vector<const Vertex *>> components;
        for (size_t i = 0; i < shorts.size(); i++)
            components[tag[i]].push_back(shorts[i]);
        for (auto &[t, vertices] : components) {
            bool needed = std::any_of(vertices.begin(), vertices.end(),
                                      [&](const Vertex *v) {return v->size() >= params.min_long;});
            if (!needed)
                continue;
            Region child;
            child.level = region.level + 1;
            child.mirrored = region.mirrored;
            std::set<std::pair<const Vertex *, bool>> seen;
            for (const Vertex *v : vertices) {
                for (const Edge &e : *v) {
                    const Vertex *x = &e.getFinish();
                    if (shortId(x) == size_t(-1) && seen.emplace(x, false).second) {
                        auto it = tag_in_flow.find({inNode(x), t});
                        child.ports.push_back({x, false, it == tag_in_flow.end() ? 0 : it->second});
                    }
                }
                for (const Edge &e : v->incoming()) {
                    const Vertex *x = &e.getStart();
                    if (shortId(x) == size_t(-1) && seen.emplace(x, true).second) {
                        auto it = tag_out_flow.find({outNode(x), t});
                        child.ports.push_back({x, true, it == tag_out_flow.end() ? 0 : it->second});
                    }
                }
            }
            child.vertices = std::move(vertices);
            res.children.push_back(std::move(child));
        }
        return res;
    }

    std::string CostStr(double cost) {
        if (std::isinf(cost))
            return "inf";
        std::ostringstream os;
        os << std::fixed << std::setprecision(1) << cost;
        return os.str();
    }
}

double CoverageModel::logDensity(double coverage, double length, size_t multiplicity) const {
    double r = 1 / (dispersion * dispersion);
    double n = coverage * length / sample_length;
    double mu = double(multiplicity) * unit_coverage * length / sample_length;
    return std::lgamma(n + r) - std::lgamma(r) - std::lgamma(n + 1) + r * std::log(r / (r + mu)) +
           n * std::log(mu / (r + mu)) + std::log(length / sample_length);
}

double CoverageModel::logPrior(size_t multiplicity) const {
    if (multiplicity == 0 || multiplicity > prior.size())
        return std::log(PRIOR_FLOOR);
    return std::log(std::max(prior[multiplicity - 1], PRIOR_FLOOR));
}

double CoverageModel::cost(double coverage, double length, size_t multiplicity) const {
    return -logDensity(coverage, length, multiplicity) - logPrior(multiplicity);
}

double CoverageModel::logLikelihood(double coverage, double length) const {
    std::vector<double> values;
    for (size_t m = 1; m <= prior.size(); m++)
        values.push_back(logPrior(m) + logDensity(coverage, length, m));
    return LogSumExp(values);
}

size_t CoverageModel::mostLikely(double coverage, double length) const {
    size_t best = 1;
    for (size_t m = 2; m <= prior.size(); m++)
        if (cost(coverage, length, m) < cost(coverage, length, best))
            best = m;
    return best;
}

std::string CoverageModel::str() const {
    std::ostringstream os;
    os << "unit coverage " << unit_coverage << ", sample length " << sample_length << ", dispersion " << dispersion
       << ", prior";
    for (size_t m = 1; m <= prior.size() && m <= 8; m++)
        os << " " << m << ":" << std::setprecision(3) << prior[m - 1];
    return os.str();
}

CoverageModel CoverageModel::Fit(logging::Logger &logger, const std::vector<std::pair<double, double>> &data,
                                 size_t max_multiplicity, double unit_hint) {
    VERIFY(!data.empty());
    double start = unit_hint > 0 ? unit_hint : CoverageMode(data);
//    Vertices far below one copy can not be explained by any multiplicity >= 1: they are errors and would only distort
//    the fit.
    double min_coverage = start / 4;
    std::vector<std::pair<double, double>> filtered;
    double excluded_length = 0;
    for (const auto &[cov, len] : data) {
        if (cov >= min_coverage)
            filtered.emplace_back(cov, len);
        else
            excluded_length += len;
    }
    logger.info() << (unit_hint > 0 ? "Given unit coverage " : "Read-mass weighted coverage mode of long vertices: ")
                  << start << ". Fitting " << filtered.size() << " vertices; " << data.size() - filtered.size()
                  << " vertices of total length " << size_t(excluded_length) << " with coverage below "
                  << min_coverage << " excluded" << std::endl;
    if (filtered.empty())
        filtered = data;
    auto [model, score] = FitEM(filtered, start, unit_hint > 0, max_multiplicity);
    logger.info() << "Coverage model: " << model.str() << ", mean log-likelihood " << score << std::endl;
    return model;
}

const MultiplicityEstimate *LongVertexMultiplicities::get(const Vertex &vertex) const {
    auto it = estimates.find(&vertex);
    return it == estimates.end() ? nullptr : &it->second;
}

LongVertexMultiplicities LongVertexMultiplicities::Estimate(logging::Logger &logger, size_t threads,
                                                            const AssemblyGraph &graph,
                                                            const MultiplicityParams &params) {
    LongVertexMultiplicities res;
    std::vector<std::pair<double, double>> data;
    size_t long_count = 0;
    for (const Vertex &v : graph.vertices()) {
        if (v.size() < params.min_long)
            continue;
        long_count++;
        if (v.isCanonical() && v.hasCoverageInfo() && v.getSPGCoverage() > 0)
            data.emplace_back(v.getSPGCoverage(), double(v.size()));
    }
    logger.info() << "Estimating multiplicities of " << long_count << " vertices of length at least "
                  << params.min_long << std::endl;
    if (data.empty()) {
        logger.info() << "No long vertices with coverage. Skipping multiplicity estimation." << std::endl;
        return res;
    }
    res.model = CoverageModel::Fit(logger, data, params.max_multiplicity, params.unit_coverage);
    BreakCount breaks = CountObviousBreaks(graph, params.min_long, res.model.unit_coverage);
//    One pseudo-break keeps the rate positive when no break is seen.
    res.break_rate = double(std::max<size_t>(breaks.dead_sides, 1)) / std::max(breaks.length, 1.0);
    double mean_length = breaks.vertices > 0 ? breaks.length / double(breaks.vertices) : double(params.min_long);
    res.end_cost = EndCost(res.break_rate, mean_length);
    res.break_cost = params.break_cost_factor * res.end_cost;
    logger.info() << "Obvious genome breaks: " << breaks.dead_sides << " dead out-sides of " << breaks.vertices
                  << " long vertices with about one copy (both strands, total length " << size_t(breaks.length)
                  << "): " << res.break_rate * 1e6 << " breaks per Mb. Genome end cost after a long vertex of mean "
                  << "length " << size_t(mean_length) << ": " << res.end_cost << " at dead ends, " << res.break_cost
                  << " elsewhere" << std::endl;
//    Error model: one "copy" at the length-weighted mean coverage of vertices below a quarter of a copy, with strong
//    overdispersion (dispersion 1: the coverage of an error depends on how many reads happen to share it).
    double error_sum = 0, error_length = 0;
    for (const Vertex &v : graph.vertices()) {
        if (!v.isCanonical() || !v.hasCoverageInfo())
            continue;
        double cov = v.getSPGCoverage();
        if (cov > 0 && cov < res.model.unit_coverage / 4) {
            error_sum += cov * double(v.size());
            error_length += double(v.size());
        }
    }
    res.error_model.unit_coverage = error_length > 0 ? error_sum / error_length : 1;
    res.error_model.sample_length = res.model.sample_length;
    res.error_model.dispersion = 1;
    res.error_model.prior = {1};
    logger.info() << "Error coverage model for dead ends: " << res.error_model.str() << std::endl;
    Context ctx{res.model, res.error_model, params, res.break_rate, res.end_cost, res.break_cost};

//    Top level: weakly connected components, one per rc pair unless the component is its own rc.
    std::vector<Region> current;
    std::unordered_set<const Vertex *> visited;
    for (const Vertex &start : graph.vertices()) {
        if (visited.count(&start))
            continue;
        Region region;
        std::vector<const Vertex *> stack = {&start};
        visited.insert(&start);
        while (!stack.empty()) {
            const Vertex *v = stack.back();
            stack.pop_back();
            region.vertices.push_back(v);
            for (const Edge &e : *v)
                if (visited.insert(&e.getFinish()).second)
                    stack.push_back(&e.getFinish());
            for (const Edge &e : v->incoming())
                if (visited.insert(&e.getStart()).second)
                    stack.push_back(&e.getStart());
        }
        std::unordered_set<const Vertex *> members(region.vertices.begin(), region.vertices.end());
        region.mirrored = members.count(&start.rc()) == 0;
        if (region.mirrored)
            for (const Vertex *v : region.vertices)
                visited.insert(&v->rc());
        bool needed = std::any_of(region.vertices.begin(), region.vertices.end(),
                                  [&](const Vertex *v) {return v->size() >= params.min_long;});
        if (needed)
            current.push_back(std::move(region));
    }

    while (!current.empty()) {
        std::vector<RegionResult> results(current.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (size_t i = 0; i < current.size(); i++)
            results[i] = ProcessRegion(current[i], ctx);
        std::vector<Region> next;
        size_t level = current.front().level;
        size_t estimated = 0, end_flow = 0, break_flow = 0, max_nodes = 0, dead_ends = 0, dead_end_flow = 0;
        size_t min_threshold = size_t(-1), max_threshold = 0;
        for (size_t i = 0; i < results.size(); i++) {
            RegionResult &r = results[i];
            for (const auto &[v, est] : r.estimates) {
                res.estimates[v] = est;
                if (current[i].mirrored)
                    res.estimates[&v->rc()] = est;
            }
//            A later level judges a dead end inside a smaller problem: it replaces the earlier judgement.
            for (const auto &[key, est] : r.dead_ends) {
                res.dead_ends[key] = est;
                if (current[i].mirrored) {
                    DeadEndEstimate mirrored = est;
                    mirrored.forward = !est.forward;
                    res.dead_ends[{&key.first->rc(), !key.second}] = mirrored;
                }
            }
            estimated += r.long_count;
            end_flow += r.end_flow;
            break_flow += r.break_flow;
            dead_ends += r.dead_ends.size();
            dead_end_flow += r.dead_end_flow;
            max_nodes = std::max(max_nodes, r.nodes);
            min_threshold = std::min(min_threshold, r.threshold);
            max_threshold = std::max(max_threshold, r.threshold);
            for (Region &child : r.children)
                next.push_back(std::move(child));
        }
        logger.info() << "Multiplicity level " << level << ": " << current.size() << " flow problems, "
                      << estimated << " long vertices, length thresholds " << min_threshold << "-" << max_threshold
                      << ", largest network " << max_nodes << " nodes, " << dead_ends
                      << " short dead ends kept, genome ends used: " << end_flow << " at long dead ends, "
                      << dead_end_flow << " at short dead ends, " << break_flow << " at breaks" << std::endl;
        current = std::move(next);
    }

    std::map<size_t, std::pair<size_t, size_t>> histogram;  // multiplicity -> (vertices, total length)
    size_t asymmetric = 0, ambiguous = 0;
    for (const auto &[v, est] : res.estimates) {
        if (!v->isCanonical())
            continue;
        histogram[est.multiplicity].first++;
        histogram[est.multiplicity].second += v->size();
        const MultiplicityEstimate *rc = res.get(v->rc());
        if (rc != nullptr && rc->multiplicity != est.multiplicity)
            asymmetric++;
        if (std::min(est.increase_cost, est.decrease_cost) < res.end_cost)
            ambiguous++;
    }
    for (const auto &[m, stat] : histogram)
        logger.info() << "Multiplicity " << m << ": " << stat.first << " vertices, total length " << stat.second
                      << std::endl;
    logger.info() << "Vertices with a change of multiplicity cheaper than one genome end: " << ambiguous
                  << "; rc pairs with different multiplicities: " << asymmetric << std::endl;
    return res;
}

void LongVertexMultiplicities::printReport(const std::experimental::filesystem::path &path,
                                           const ReferenceCounter &reference) const {
    std::vector<const Vertex *> vertices;
    for (const auto &[v, est] : estimates)
        vertices.push_back(v);
    std::sort(vertices.begin(), vertices.end(), [](const Vertex *a, const Vertex *b) {return *a < *b;});
    std::ofstream os(path);
    os << "vertex\tlength\tcoverage\texpected\tmultiplicity\tincrease_cost\tdecrease_cost\tlevel\tthreshold";
    if (reference)
        os << "\ttrue_multiplicity\ttraversals_at_reference_ends\tmin_covered_fraction";
    os << "\n";
    for (const Vertex *v : vertices) {
        const MultiplicityEstimate &est = estimates.at(v);
        os << v->getInnerId() << "\t" << v->size() << "\t"
           << (v->hasCoverageInfo() ? std::to_string(v->getSPGCoverage()) : "-") << "\t" << est.expected << "\t"
           << est.multiplicity << "\t" << CostStr(est.increase_cost) << "\t" << CostStr(est.decrease_cost) << "\t"
           << est.level << "\t" << est.threshold;
        if (reference) {
            ReferenceCount count = reference(*v);
            os << "\t" << count.traversals << "\t" << count.at_ends << "\t" << count.min_covered;
        }
        os << "\n";
    }
}

LongVertexMultiplicities::ReferenceCounter LongVertexMultiplicities::FromPaths(DLLAlignmentStorage &paths) {
    auto counts = std::make_shared<std::unordered_map<ConstVertexId, ReferenceCount>>();
    for (const auto &[vertex, traversals] : paths.vertexTraversals()) {
        ReferenceCount &count = (*counts)[vertex];
        for (const DLLAlignmentStorage::VertexTraversal &t : traversals) {
            (t.at_end ? count.at_ends : count.traversals)++;
            count.min_covered = std::min(count.min_covered, double(t.covered) / double(vertex->size()));
        }
    }
    return [counts](const Vertex &v) {
        auto it = counts->find(v.getId());
        return it == counts->end() ? ReferenceCount() : it->second;
    };
}

void LongVertexMultiplicities::compareWithReference(logging::Logger &logger,
                                                    const ReferenceCounter &reference) const {
    struct Stat {
        size_t vertices = 0;
        size_t length = 0;
        void add(const Vertex &v) {vertices++; length += v.size();}
    };
    std::map<std::pair<size_t, size_t>, Stat> table;  // (estimated, true) -> stat, vertices away from reference ends
    Stat correct, high, low, flagged, at_ends, absent;
    const std::vector<double> covered_bounds = {0.1, 0.5, 0.9, 1.0};
    std::vector<Stat> by_covered(covered_bounds.size() + 1);  // compared vertices by smallest covered fraction
    for (const auto &[v, est] : estimates) {
        if (!v->isCanonical())
            continue;
        ReferenceCount count = reference(*v);
        size_t full = count.traversals;
        if (count.at_ends > 0) {
            at_ends.add(*v);
            continue;
        }
        if (full == 0)
            absent.add(*v);
        else
            by_covered[std::upper_bound(covered_bounds.begin(), covered_bounds.end(), count.min_covered + 1e-9) -
                       covered_bounds.begin()].add(*v);
        table[{est.multiplicity, full}].add(*v);
        if (est.multiplicity == full) {
            correct.add(*v);
            continue;
        }
        (est.multiplicity > full ? high : low).add(*v);
        double towards_truth = est.multiplicity > full ? est.decrease_cost : est.increase_cost;
        if (towards_truth < end_cost)
            flagged.add(*v);
    }
    auto str = [](const Stat &s) {return std::to_string(s.vertices) + " vertices (" + std::to_string(s.length) + " bp)";};
    logger.info() << "Multiplicities against reference: correct " << str(correct) << ", too high " << str(high)
                  << ", too low " << str(low) << "; wrong ones flagged as ambiguous: " << str(flagged)
                  << "; absent from reference: " << str(absent) << "; touching a reference end (not compared): "
                  << str(at_ends) << std::endl;
    logger.info() << "Smallest fraction of a vertex spelled by one reference traversal: <0.1 " << str(by_covered[0])
                  << ", 0.1-0.5 " << str(by_covered[1]) << ", 0.5-0.9 " << str(by_covered[2]) << ", 0.9-1 "
                  << str(by_covered[3]) << ", 1 " << str(by_covered[4]) << std::endl;
    for (const auto &[key, stat] : table)
        logger.info() << "Estimated " << key.first << ", true " << key.second << ": " << str(stat) << std::endl;
}

void LongVertexMultiplicities::printDeadEnds(const std::experimental::filesystem::path &path,
                                             const ReferenceCounter &reference) const {
    std::vector<std::pair<const Vertex *, bool>> keys;
    for (const auto &[key, est] : dead_ends)
        keys.push_back(key);
    std::sort(keys.begin(), keys.end(), [](const auto &a, const auto &b) {
        return *a.first < *b.first || (*a.first == *b.first && a.second < b.second);
    });
    std::ofstream os(path);
    os << "vertex\ttype\tlength\tcoverage\texpected\tp_real\tcost\tflow\tlevel";
    if (reference)
        os << "\treference_traversals";
    os << "\n";
    for (const auto &key : keys) {
        const DeadEndEstimate &est = dead_ends.at(key);
        os << key.first->getInnerId() << "\t" << (est.forward ? "end" : "start") << "\t" << est.length << "\t"
           << est.coverage << "\t" << est.coverage / model.unit_coverage << "\t" << est.p_real << "\t"
           << CostStr(est.cost) << "\t" << est.flow << "\t" << est.level;
        if (reference) {
            ReferenceCount count = reference(*key.first);
            os << "\t" << count.traversals + count.at_ends;
        }
        os << "\n";
    }
}

VertexInfo LongVertexMultiplicities::vertexInfo() const {
    const std::unordered_map<const Vertex *, MultiplicityEstimate> *map = &estimates;
    const std::map<std::pair<const Vertex *, bool>, DeadEndEstimate> *ends = &dead_ends;
    return VertexInfo::Labeler([map, ends](const Vertex &v) -> std::string {
        std::string res;
        auto it = map->find(&v);
        if (it != map->end()) {
            const MultiplicityEstimate &est = it->second;
            res = "Mult:" + std::to_string(est.multiplicity) + " (+" + CostStr(est.increase_cost) + "/-" +
                  CostStr(est.decrease_cost) + ")";
        }
        for (bool forward : {false, true}) {
            auto end = ends->find({&v, forward});
            if (end == ends->end())
                continue;
            std::ostringstream os;
            os << (forward ? "End" : "Start") << " p=" << std::setprecision(2) << end->second.p_real << " used="
               << end->second.flow;
            res += (res.empty() ? "" : "\\n") + os.str();
        }
        return res;
    });
}
