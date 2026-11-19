#include "fold.hpp"

#include <cstdio>
#include <numeric>
#include <ostream>
#include <regex>
#include <sstream>
#include <unordered_set>

namespace dot_viewer {

    namespace {
        class Dsu {
        public:
            explicit Dsu(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }

            size_t find(size_t x) {
                while (parent[x] != x) {
                    parent[x] = parent[parent[x]];
                    x = parent[x];
                }
                return x;
            }

            void unite(size_t a, size_t b) { parent[find(a)] = find(b); }

        private:
            std::vector<size_t> parent;
        };

        bool extractNumber(const std::regex &re, const std::string &text, double &value) {
            std::smatch m;
            if (!std::regex_search(text, m, re) || m.size() < 2)
                return false;
            try {
                value = std::stod(m[1].str());
            } catch (const std::exception &) {
                return false;
            }
            return true;
        }

        std::string fmtLen(double len) {
            char buf[64];
            snprintf(buf, sizeof(buf), "%.0f", len);
            return buf;
        }

        std::string fmtCov(double cov) {
            char buf[64];
            snprintf(buf, sizeof(buf), "%.2f", cov);
            return buf;
        }

        /// Multi-line description of aggregated stats; `sep` is the dot line separator ("\\n") or "\n".
        std::string statsText(const FoldStats &s, const std::string &sep) {
            std::string res = std::to_string(s.vertices) + " vertices";
            if (s.long_vertices > 0)
                res += " (" + std::to_string(s.long_vertices) + " long)";
            if (s.has_spelled)
                res += sep + "len " + fmtLen(s.spelled);
            res += sep + "sum len " + fmtLen(s.sum_len) + ", max " + fmtLen(s.max_len);
            res += sep + "cov " + fmtCov(s.coverage());
            return res;
        }

        const std::vector<DotAttr> CONTEXT_STYLE = {
                makeAttr("style", "filled,dashed"),
                makeAttr("fillcolor", "#eeeeee"),
                makeAttr("fontcolor", "#555555"),
        };
    }

    void FoldStats::add(const FoldStats &other) {
        vertices += other.vertices;
        long_vertices += other.long_vertices;
        sum_len += other.sum_len;
        max_len = std::max(max_len, other.max_len);
        sum_len_cov += other.sum_len_cov;
        sum_cov += other.sum_cov;
    }

    double FoldStats::coverage() const {
        if (sum_len > 0)
            return sum_len_cov / sum_len;
        return vertices == 0 ? 0 : sum_cov / double(vertices);
    }

    FoldedGraph::FoldedGraph(const DotGraph &graph_, const FoldParams &params) : graph(graph_) {
        computeMetrics(params);
        buildLevel1();
        buildLevel2();
    }

    void FoldedGraph::computeMetrics(const FoldParams &params) {
        std::regex len_re(params.length_regex), cov_re(params.coverage_regex), trunc_re(params.trunc_regex);
        long_threshold = params.long_threshold;
        vertex_stats.resize(graph.nodeCount());
        for (size_t v = 0; v < graph.nodeCount(); ++v) {
            const DotAttr *label = findAttr(graph.node(v).attrs, "label");
            double len = 0, cov = 0;
            if (label) {
                extractNumber(len_re, label->value, len);
                extractNumber(cov_re, label->value, cov);
            }
            FoldStats &s = vertex_stats[v];
            s.vertices = 1;
            s.long_vertices = len >= params.long_threshold ? 1 : 0;
            s.sum_len = s.max_len = len;
            s.sum_len_cov = len * cov;
            s.sum_cov = cov;
            long_count += s.long_vertices;
            max_vertex_length = std::max(max_vertex_length, len);
        }
        has_trunc.assign(graph.edgeCount(), false);
        edge_trunc.assign(graph.edgeCount(), 0);
        for (size_t e = 0; e < graph.edgeCount(); ++e) {
            const DotAttr *label = findAttr(graph.edge(e).attrs, "label");
            double t = 0;
            if (label && extractNumber(trunc_re, label->value, t)) {
                has_trunc[e] = true;
                edge_trunc[e] = t;
            }
        }
    }

    void FoldedGraph::buildLevel1() {
        const size_t n = graph.nodeCount();
        auto isLong = [this](size_t v) { return vertex_stats[v].long_vertices > 0; };
        // Cut-graph node ids: v is the vertex itself (or the start half of a long vertex), n + v is the end half.
        auto fromSide = [&](size_t e) { size_t a = graph.edge(e).from; return isLong(a) ? n + a : a; };

        // Hanging pieces: connected sets of short vertices whose only long neighbour is a single long vertex,
        // touched at both its start and its end. Cutting would not separate the two sides of that long vertex,
        // so each piece is attached to one side only and its edges to the other side are detached.
        Dsu short_dsu(n);
        for (size_t e = 0; e < graph.edgeCount(); ++e)
            if (!isLong(graph.edge(e).from) && !isLong(graph.edge(e).to))
                short_dsu.unite(graph.edge(e).from, graph.edge(e).to);
        const size_t NONE = size_t(-1), MANY = size_t(-2);
        std::vector<size_t> long_neighbour(n, NONE);
        std::vector<bool> touches_start(n, false), touches_end(n, false);
        auto note = [&](size_t piece, size_t long_vertex, bool start) {
            size_t &nb = long_neighbour[piece];
            nb = (nb == NONE || nb == long_vertex) ? long_vertex : MANY;
            (start ? touches_start : touches_end)[piece] = true;
        };
        for (size_t e = 0; e < graph.edgeCount(); ++e) {
            size_t a = graph.edge(e).from, b = graph.edge(e).to;
            if (!isLong(a) && isLong(b))
                note(short_dsu.find(a), b, true);
            else if (isLong(a) && !isLong(b))
                note(short_dsu.find(b), a, false);
        }
        // The hanging piece of an edge between a piece and its long vertex, or NONE.
        auto hangingPiece = [&](size_t e) {
            size_t a = graph.edge(e).from, b = graph.edge(e).to;
            if (isLong(a) == isLong(b))
                return NONE;
            size_t piece = short_dsu.find(isLong(a) ? b : a);
            bool hanging = long_neighbour[piece] < MANY && touches_start[piece] && touches_end[piece];
            return hanging ? piece : NONE;
        };

        Dsu dsu(2 * n);
        for (size_t e = 0; e < graph.edgeCount(); ++e)
            if (hangingPiece(e) == NONE)
                dsu.unite(fromSide(e), graph.edge(e).to);
        // Attach every piece to the side of its long vertex whose component already has more short vertices
        // (rc-symmetric except for ties, which go to the end side).
        std::vector<size_t> short_count(2 * n, 0);
        for (size_t v = 0; v < n; ++v)
            if (!isLong(v))
                ++short_count[dsu.find(v)];
        std::vector<bool> piece_to_end(n, false);
        for (size_t v = 0; v < n; ++v) {
            size_t l = long_neighbour[v];
            if (!isLong(v) && short_dsu.find(v) == v && l < MANY)
                piece_to_end[v] = short_count[dsu.find(n + l)] >= short_count[dsu.find(l)];
        }
        std::vector<bool> detached(graph.edgeCount(), false);
        for (size_t e = 0; e < graph.edgeCount(); ++e) {
            size_t piece = hangingPiece(e);
            if (piece == NONE)
                continue;
            bool edge_at_end = isLong(graph.edge(e).from); // long -> short leaves from the long vertex's end
            if (edge_at_end == piece_to_end[piece])
                dsu.unite(fromSide(e), graph.edge(e).to);
            else
                detached[e] = true;
        }
        std::vector<long long> comp_of_root(2 * n, -1);
        auto compOf = [&](size_t x) -> Component & {
            size_t r = dsu.find(x);
            if (comp_of_root[r] < 0) {
                comp_of_root[r] = (long long) components.size();
                components.emplace_back();
            }
            return components[size_t(comp_of_root[r])];
        };
        for (size_t v = 0; v < n; ++v) {
            if (isLong(v)) {
                compOf(v).starts.push_back(v);
                compOf(n + v).ends.push_back(v);
            } else {
                Component &c = compOf(v);
                c.members.push_back(v);
                c.stats.add(vertex_stats[v]);
                // Part of a piece hanging on a single long vertex (at one or both of its ends).
                if (long_neighbour[short_dsu.find(v)] < MANY)
                    c.has_hanging = true;
            }
        }
        for (size_t e = 0; e < graph.edgeCount(); ++e) {
            if (detached[e]) {
                size_t piece_vertex = isLong(graph.edge(e).from) ? graph.edge(e).to : graph.edge(e).from;
                compOf(piece_vertex).edges.push_back(e);
            } else {
                compOf(graph.edge(e).to).edges.push_back(e);
            }
        }

        // Level-1 elements. Components that are not folded are dropped from `components` afterwards.
        std::vector<long long> elem_of_vertex(n, -1);
        for (size_t v = 0; v < n; ++v) {
            if (isLong(v)) {
                elem_of_vertex[v] = (long long) l1.size();
                l1.push_back({false, v});
            }
        }
        std::vector<Component> folded;
        std::vector<long long> folded_elem(components.size(), -1);
        for (size_t c = 0; c < components.size(); ++c) {
            // A hanging piece is always folded into its side: kept as a bare vertex it would add edges to the
            // long vertex and turn it into a fork.
            if (components[c].members.size() == 1 && !components[c].has_hanging) {
                size_t v = components[c].members[0];
                elem_of_vertex[v] = (long long) l1.size();
                l1.push_back({false, v});
            } else if (!components[c].members.empty()) {
                folded_elem[c] = (long long) l1.size();
                l1.push_back({true, folded.size()});
                folded.push_back(components[c]);
            }
        }
        for (size_t c = 0; c < components.size(); ++c) {
            const Component &comp = components[c];
            if (folded_elem[c] < 0) {
                for (size_t e: comp.edges)
                    l1_edges.push_back({size_t(elem_of_vertex[graph.edge(e).from]),
                                        size_t(elem_of_vertex[graph.edge(e).to]), (long long) e});
            } else {
                size_t f = size_t(folded_elem[c]);
                for (size_t v: comp.starts)
                    l1_edges.push_back({f, size_t(elem_of_vertex[v]), -1});
                for (size_t v: comp.ends)
                    l1_edges.push_back({size_t(elem_of_vertex[v]), f, -1});
            }
        }
        components = std::move(folded);
        l1_out.assign(l1.size(), {});
        l1_in.assign(l1.size(), {});
        for (size_t i = 0; i < l1_edges.size(); ++i) {
            l1_out[l1_edges[i].from].push_back(i);
            l1_in[l1_edges[i].to].push_back(i);
        }
    }

    void FoldedGraph::buildLevel2() {
        std::vector<bool> fork(l1.size());
        for (size_t x = 0; x < l1.size(); ++x)
            fork[x] = l1_out[x].size() > 1 || l1_in[x].size() > 1;
        Dsu dsu(l1.size());
        for (const L1Edge &e: l1_edges)
            if (!fork[e.from] && !fork[e.to])
                dsu.unite(e.from, e.to);
        std::vector<std::vector<size_t>> members(l1.size());
        for (size_t x = 0; x < l1.size(); ++x)
            if (!fork[x])
                members[dsu.find(x)].push_back(x);
        group_of_l1.assign(l1.size(), 0);
        for (size_t x = 0; x < l1.size(); ++x) {
            if (fork[x] || members[dsu.find(x)].size() == 1) {
                group_of_l1[x] = groups.size();
                groups.push_back({false, x});
                continue;
            }
            std::vector<size_t> &group = members[dsu.find(x)];
            if (group.front() != x)
                continue; // the path was created when its first member was visited
            std::unordered_set<size_t> in_group(group.begin(), group.end());
            Path path;
            size_t head = group.front();
            path.cycle = true;
            for (size_t y: group) {
                if (l1_in[y].empty() || !in_group.count(l1_edges[l1_in[y][0]].from)) {
                    head = y;
                    path.cycle = false;
                    break;
                }
            }
            std::unordered_set<size_t> visited;
            for (size_t y = head; visited.insert(y).second;) {
                path.elems.push_back(y);
                if (l1_out[y].empty() || !in_group.count(l1_edges[l1_out[y][0]].to))
                    break;
                y = l1_edges[l1_out[y][0]].to;
            }
            // Spelled length is only defined when every element is an original vertex joined by original edges.
            bool spellable = true;
            double spelled = path.cycle ? 0 : vertex_stats[l1[path.elems[0]].id].sum_len;
            for (size_t i = 0; i < path.elems.size(); ++i) {
                size_t y = path.elems[i];
                path.stats.add(l1Stats(y));
                if (l1[y].folded) {
                    spellable = false;
                    continue;
                }
                bool last = i + 1 == path.elems.size();
                if (last && !path.cycle)
                    continue;
                long long e = l1_edges[l1_out[y][0]].orig_edge;
                if (e < 0 || !has_trunc[size_t(e)])
                    spellable = false;
                else
                    spelled += edge_trunc[size_t(e)];
            }
            path.stats.has_spelled = spellable;
            path.stats.spelled = spelled;
            size_t p = paths.size();
            paths.push_back(std::move(path));
            size_t g = groups.size();
            groups.push_back({true, p});
            for (size_t y: group)
                group_of_l1[y] = g;
        }
    }

    const FoldStats &FoldedGraph::l1Stats(size_t x) const {
        return l1[x].folded ? components[l1[x].id].stats : vertex_stats[l1[x].id];
    }

    std::string FoldedGraph::l1Name(size_t x) const {
        if (l1[x].folded)
            return quoteDotId("__fold_C" + std::to_string(l1[x].id));
        return quoteDotId(graph.node(l1[x].id).name);
    }

    std::string FoldedGraph::groupName(size_t g) const {
        if (groups[g].is_path)
            return quoteDotId("__fold_P" + std::to_string(groups[g].id));
        return l1Name(groups[g].id);
    }

    std::string FoldedGraph::componentTitle(size_t c) const {
        return "Component C" + std::to_string(c) + ": " + statsText(components[c].stats, ", ");
    }

    std::string FoldedGraph::pathTitle(size_t p) const {
        const Path &path = paths[p];
        return std::string(path.cycle ? "Cycle P" : "Path P") + std::to_string(p) + ": " +
               std::to_string(path.elems.size()) + " elements, " + statsText(path.stats, ", ");
    }

    void FoldedGraph::writeL1Node(std::ostream &os, size_t x, bool context) const {
        if (!l1[x].folded) {
            graph.writeNode(os, l1[x].id, context ? CONTEXT_STYLE : std::vector<DotAttr>{});
            return;
        }
        size_t c = l1[x].id;
        std::string text = "C" + std::to_string(c) + "\\n" + statsText(components[c].stats, "\\n");
        std::vector<DotAttr> attrs = {
                {"label", "\"" + text + "\"", text},
                makeAttr("shape", "doubleoctagon"),
                makeAttr("style", context ? "filled,dashed" : "filled"),
                makeAttr("fillcolor", context ? "#f6e0d2" : "#f4b183"),
                makeAttr("tooltip", "Level 1: folded component C" + std::to_string(c) +
                                    " (click to show it in the original graph)"),
                makeAttr("id", "c" + std::to_string(c)),
        };
        os << l1Name(x) << " " << attrListString(attrs) << "\n";
    }

    void FoldedGraph::writePathNode(std::ostream &os, size_t p) const {
        const Path &path = paths[p];
        size_t folded_count = 0;
        for (size_t y: path.elems)
            folded_count += l1[y].folded ? 1 : 0;
        std::string text = std::string(path.cycle ? "cycle P" : "P") + std::to_string(p) + "\\n" +
                           std::to_string(path.elems.size()) + " elements";
        if (folded_count > 0)
            text += " (" + std::to_string(folded_count) + " folded)";
        text += "\\n" + statsText(path.stats, "\\n");
        std::vector<DotAttr> attrs = {
                {"label", "\"" + text + "\"", text},
                makeAttr("shape", "box"),
                makeAttr("style", "filled,rounded,bold"),
                makeAttr("peripheries", "2"),
                makeAttr("fillcolor", "#9dc3e6"),
                makeAttr("tooltip", "Level 2: folded unbranching path P" + std::to_string(p) + " (click to open)"),
                makeAttr("id", "p" + std::to_string(p)),
        };
        os << quoteDotId("__fold_P" + std::to_string(p)) << " " << attrListString(attrs) << "\n";
    }

    void FoldedGraph::writeL1Edge(std::ostream &os, const L1Edge &edge, const std::string &from,
                                  const std::string &to) const {
        if (edge.orig_edge >= 0) {
            graph.writeEdge(os, size_t(edge.orig_edge), from, to, {});
            return;
        }
        os << from << (graph.isDirected() ? " -> " : " -- ") << to << " [\"color\"=\"#7f6000\" \"style\"=bold]\n";
    }

    std::string FoldedGraph::topDot(bool hide_isolated, size_t &hidden) const {
        std::vector<bool> has_edge(groups.size(), false);
        for (const L1Edge &e: l1_edges) {
            size_t a = group_of_l1[e.from], b = group_of_l1[e.to];
            if (a != b || !groups[a].is_path)
                has_edge[a] = has_edge[b] = true;
        }
        hidden = 0;
        std::ostringstream os;
        os << graph.dotHeader();
        for (size_t g = 0; g < groups.size(); ++g) {
            bool original_vertex = !groups[g].is_path && !l1[groups[g].id].folded;
            if (hide_isolated && !has_edge[g] && original_vertex) {
                ++hidden;
                continue;
            }
            if (groups[g].is_path)
                writePathNode(os, groups[g].id);
            else
                writeL1Node(os, groups[g].id, false);
        }
        for (const L1Edge &e: l1_edges) {
            size_t a = group_of_l1[e.from], b = group_of_l1[e.to];
            if (a == b && groups[a].is_path)
                continue; // inside a path
            writeL1Edge(os, e, groupName(a), groupName(b));
        }
        os << "}\n";
        return os.str();
    }

    std::vector<size_t> FoldedGraph::componentVertices(size_t c) const {
        const Component &comp = components[c];
        std::vector<size_t> res = comp.members;
        std::unordered_set<size_t> attached;
        for (const std::vector<size_t> *side: {&comp.starts, &comp.ends})
            for (size_t v: *side)
                if (attached.insert(v).second)
                    res.push_back(v);
        return res;
    }

    std::string FoldedGraph::pathDot(size_t p) const {
        const Path &path = paths[p];
        std::unordered_set<size_t> members(path.elems.begin(), path.elems.end());
        std::unordered_set<size_t> context;
        std::vector<size_t> edges_to_draw;
        for (size_t y: path.elems) {
            for (const std::vector<size_t> *adj: {&l1_in[y], &l1_out[y]}) {
                for (size_t i: *adj) {
                    size_t other = l1_edges[i].from == y ? l1_edges[i].to : l1_edges[i].from;
                    if (!members.count(other))
                        context.insert(other);
                    // An edge between two members is reached from both ends; draw it from its source only.
                    if (!members.count(other) || l1_edges[i].from == y)
                        edges_to_draw.push_back(i);
                }
            }
        }
        std::ostringstream os;
        os << graph.dotHeader();
        for (size_t y: path.elems)
            writeL1Node(os, y, false);
        for (size_t y: context)
            writeL1Node(os, y, true);
        for (size_t i: edges_to_draw)
            writeL1Edge(os, l1_edges[i], l1Name(l1_edges[i].from), l1Name(l1_edges[i].to));
        os << "}\n";
        return os.str();
    }
}
