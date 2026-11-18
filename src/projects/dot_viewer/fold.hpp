#pragma once

#include "dot_graph.hpp"

#include <string>
#include <vector>

namespace dot_viewer {

    struct FoldParams {
        double long_threshold = 50000;
        /// Applied to the decoded "label" attribute of vertices; capture group 1 is the number.
        std::string length_regex = R"(\\n\s*(\d+))";
        std::string coverage_regex = R"(Cov:\s*([0-9.eE+-]+))";
        /// Applied to edge labels: sequence an edge adds beyond its start vertex (truncSize in ag::Printer dumps).
        std::string trunc_regex = R"(^\S+\s+(\d+))";
    };

    /// Aggregated length/coverage of a set of original vertices.
    struct FoldStats {
        size_t vertices = 0;
        size_t long_vertices = 0;
        double sum_len = 0;
        double max_len = 0;
        double sum_len_cov = 0;
        double sum_cov = 0;
        /// Length of the sequence spelled by a path of original vertices: |first| + sum of edge truncs.
        bool has_spelled = false;
        double spelled = 0;

        void add(const FoldStats &other);
        /// Length-weighted mean coverage (plain mean if no lengths are known).
        double coverage() const;
    };

    /// Two-level folding of a DotGraph.
    ///
    /// Level 1: every long vertex (length >= threshold) is cut into a start half (keeps in-edges) and an end half
    /// (keeps out-edges). Each connected component of the cut graph with >= 2 short vertices becomes a component
    /// vertex C; components with <= 1 short vertex stay unfolded (their original vertices/edges are kept).
    /// Long vertices are reattached: C -> long if long's start half is in C, long -> C if its end half is.
    /// Hanging pieces are connected sets of short vertices whose only long neighbour is a single long vertex.
    /// A component containing one is always folded, even with a single short vertex: kept bare, the piece would
    /// add an edge to the long vertex and turn it into a fork. A piece touching both ends of its long vertex
    /// would glue the two halves back together, so it is attached to one side only (the one whose component has
    /// more short vertices; end on ties) and its edges to the other side are kept as internal component edges.
    ///
    /// Level 2: in the level-1 graph, vertices with indegree or outdegree > 1 are forks. Connected components
    /// of the non-fork vertices (unbranching paths or cycles) with >= 2 elements become path vertices P.
    class FoldedGraph {
    public:
        FoldedGraph(const DotGraph &graph, const FoldParams &params);

        size_t componentCount() const { return components.size(); }
        size_t pathCount() const { return paths.size(); }
        size_t longCount() const { return long_count; }
        size_t topSize() const { return groups.size(); }

        /// The fully folded graph. With `hide_isolated`, vertices without edges are left out and counted in `hidden`.
        std::string topDot(bool hide_isolated, size_t &hidden) const;
        /// Original vertices of component c: its short vertices followed by the attached long vertices.
        std::vector<size_t> componentVertices(size_t c) const;
        /// Level-1 elements of path p in order with their outside neighbours as dashed context.
        std::string pathDot(size_t p) const;

        std::string componentTitle(size_t c) const;
        std::string pathTitle(size_t p) const;

    private:
        struct Component {
            std::vector<size_t> members; // short original vertices
            std::vector<size_t> starts;  // long vertices whose start half is in the component
            std::vector<size_t> ends;    // long vertices whose end half is in the component
            std::vector<size_t> edges;   // original edges inside the cut graph component (incl. detached ones)
            bool has_hanging = false;    // contains a piece of short vertices hanging on a single long vertex
            FoldStats stats;
        };

        /// Element of the level-1 graph: an original vertex or a folded component.
        struct L1Elem {
            bool folded;
            size_t id;
        };

        struct L1Edge {
            size_t from;
            size_t to;
            long long orig_edge; // -1 for synthetic component <-> long vertex edges
        };

        struct Path {
            std::vector<size_t> elems; // level-1 elements in walk order
            bool cycle = false;
            FoldStats stats;
        };

        /// Vertex of the top graph: a single level-1 element or a path.
        struct Group {
            bool is_path;
            size_t id;
        };

        const DotGraph &graph;
        std::vector<FoldStats> vertex_stats;
        std::vector<bool> has_trunc;
        std::vector<double> edge_trunc;
        size_t long_count = 0;

        std::vector<Component> components;
        std::vector<L1Elem> l1;
        std::vector<L1Edge> l1_edges;
        std::vector<std::vector<size_t>> l1_out, l1_in;
        std::vector<Path> paths;
        std::vector<Group> groups;
        std::vector<size_t> group_of_l1;

        void computeMetrics(const FoldParams &params);
        void buildLevel1();
        void buildLevel2();

        const FoldStats &l1Stats(size_t x) const;
        std::string l1Name(size_t x) const;
        std::string groupName(size_t g) const;
        void writeL1Node(std::ostream &os, size_t x, bool context) const;
        void writePathNode(std::ostream &os, size_t p) const;
        void writeL1Edge(std::ostream &os, const L1Edge &edge, const std::string &from, const std::string &to) const;
    };
}
