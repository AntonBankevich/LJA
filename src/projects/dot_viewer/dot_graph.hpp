#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dot_viewer {

    /// One attribute of a node/edge/graph exactly as written in the source dot file.
    /// `raw` is the value token verbatim (including quotes if it was quoted), so it can be
    /// re-emitted without any loss; `value` is the unquoted/unescaped value used for searching.
    struct DotAttr {
        std::string key;
        std::string raw;
        std::string value;
    };

    struct DotNode {
        std::string name;
        std::vector<DotAttr> attrs;
        std::vector<size_t> out_edges;
        std::vector<size_t> in_edges;
    };

    struct DotEdge {
        size_t from;
        size_t to;
        std::vector<DotAttr> attrs;
    };

    /// Element of the graph that a query token can refer to.
    struct DotElement {
        bool is_edge;
        size_t idx;
        bool operator==(const DotElement &other) const { return is_edge == other.is_edge && idx == other.idx; }
    };

    /// Result of a query: a set of vertices to draw plus edges to emphasize.
    struct Selection {
        std::vector<size_t> nodes;
        std::vector<size_t> highlighted_edges;
        bool truncated = false;
        std::string message;
    };

    /// In-memory copy of a (possibly huge) dot file that supports drawing induced subgraphs.
    ///
    /// Supported input: the flat "strict? (di)graph name? { ... }" layout produced by ag::Printer and
    /// most other tools: node statements, edge chains (a -> b -> c), graph attribute assignments,
    /// graph/node/edge default attribute statements and subgraphs (subgraph boundaries are flattened,
    /// defaults inside subgraphs are ignored). Node ports ("a:p") are accepted and dropped.
    ///
    /// The object is immutable after construction and therefore safe to query from many threads.
    class DotGraph {
    public:
        explicit DotGraph(const std::string &path);

        size_t nodeCount() const { return nodes.size(); }
        size_t edgeCount() const { return edges.size(); }
        const DotNode &node(size_t idx) const { return nodes[idx]; }
        const DotEdge &edge(size_t idx) const { return edges[idx]; }

        /// Elements whose name, tooltip, or first label word/line equals the token.
        const std::vector<DotElement> &lookup(const std::string &token) const;

        /// All vertices within `radius` undirected steps of the vertices matched by the tokens in `query`.
        /// Matched edges contribute both endpoints. Stops growing once `limit` vertices are collected.
        Selection selectAround(const std::string &query, size_t radius, size_t limit) const;

        /// Resolves `query` (a sequence of vertex/edge tokens) to a directed walk, using adjacency to
        /// disambiguate tokens that match several elements, then adds a neighbourhood of `radius`.
        /// Consecutive vertex tokens must be joined by an edge; every edge between them is highlighted.
        Selection selectPath(const std::string &query, size_t radius, size_t limit) const;

        /// `current` plus all neighbours (in and out) of `vertex`, capped at `limit` new vertices.
        Selection expand(const std::vector<size_t> &current, size_t vertex, size_t limit) const;

        /// Vertices of `nodes_set` that have at least one neighbour outside of it.
        std::vector<size_t> frontier(const std::vector<size_t> &nodes_set) const;

        /// Dot text of the subgraph induced by `nodes_set` with all original attributes. Frontier vertices
        /// are filled yellow, highlighted edges get a thicker pen. Every node receives id="v<idx>" so that
        /// the browser can map SVG elements back to vertices.
        std::string inducedDot(const std::vector<size_t> &nodes_set, const std::vector<size_t> &highlighted) const;

    private:
        bool directed = true;
        std::string graph_name;
        std::vector<std::string> graph_statements; // top-level "k = v", "graph [...]", "node [...]", "edge [...]"
        std::vector<DotNode> nodes;
        std::vector<DotEdge> edges;
        std::unordered_map<std::string, size_t> name_to_node;
        std::unordered_map<std::string, std::vector<DotElement>> index;

        friend class DotParser;

        size_t getOrCreateNode(const std::string &name);
        void buildIndex();
        void addIndexKey(const std::string &key, DotElement element);
        std::vector<std::vector<DotElement>> resolveTokens(const std::string &query, std::string &error) const;
        void growBall(std::vector<size_t> &result, std::unordered_set<size_t> &taken, size_t radius, size_t limit,
                      bool &truncated) const;
    };

    std::string quoteDotId(const std::string &id);
}
