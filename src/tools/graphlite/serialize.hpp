// Copyright 2021 Guohao Dou
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef GRAPHLITE_SERIALIZE_H
#define GRAPHLITE_SERIALIZE_H

#include "graphlite.hpp"
#include <deque>
#include "optional"

namespace graph_lite {
    namespace detail {
        // serialize properties if they are map/unordered_map
        template<typename M>
        void serialize_map_like(std::ostream& os, const M& m) {
            static_assert(is_either_map_v<M>);
            static_assert(is_streamable_v<typename M::key_type>
                          and is_streamable_v<typename M::mapped_type>,
                          "both key and value in property pair key=value should be serializable");
            typename M::size_type size = 0;
            for (const auto& p: m) {
                os << p.first << '=' << '"' << p.second << '"';
                ++size;
                if (size!=m.size()) {
                    os << ", ";
                }
            }
        }

        template <typename PT>
        struct formatter {
            std::optional<std::function<std::string(const PT&)>> fmt;
        };
        template<>
        struct formatter<void> {};
    }

    template<typename NodeType, typename NodePropType, typename EdgePropType,
            EdgeDirection direction, MultiEdge multi_edge, SelfLoop self_loop,
            Map adj_list_spec, Container neighbors_container_spec>
    class Serializer {
    public:
        using GType = Graph<NodeType, NodePropType, EdgePropType, direction, multi_edge, self_loop, adj_list_spec, neighbors_container_spec>;
    private:  // formatting settings
        // no line wrap by default
        size_t max_num_nodes_per_line = (size_t)-1;
        size_t max_num_edges_per_line = (size_t)-1;
    private:
        void add_indent(std::ostream& os, size_t indent) const {
            for (size_t i = 0; i < indent; ++i) {
                os << '\t';
            }
        }

        auto get_neighbors(const NodeType& node) const {
            if constexpr(direction==EdgeDirection::UNDIRECTED) {
                return graph.neighbors(node);
            } else {
                return graph.out_neighbors(node);
            }
        }

        template<bool is_node>
        void message() const {
            using PT = std::conditional_t<is_node, NodePropType, EdgePropType>;
            std::string node_or_edge = is_node ? "node" : "edge";
            const auto& fmt = [this]() -> auto& {
                if constexpr(is_node) {
                    return node_fmt;
                } else {
                    return edge_fmt;
                }
            }();
            if constexpr(std::is_void_v<PT>) {
                // std::cerr << "Serializer: " << "no " << node_or_edge << " property needed at all\n";
            } else if (fmt.fmt.has_value()) {
                // std::cerr << "Serializer: " << "using " << node_or_edge << " formatter provided\n";
            } else if constexpr(detail::is_either_map_v<PT>) {
                // std::cerr << "Serializer: " << node_or_edge << " prop is a map/unordered_map; using map serializer\n";
            } else if constexpr(detail::is_streamable_v<PT>) {
                // std::cerr << "Serializer: " << node_or_edge + " is by itself serializable; populating the field \"label\"\n";
            } else {
                throw std::runtime_error("failed to serialize " + node_or_edge + " properties");
            }
        }

        void serialize_nodes(std::ostream& os) const {
            if (graph.size()) {
                add_indent(os, 1);
            } else { return; }
            size_t node_count = 0;
            typename GType::ConstIterator begin = graph.begin();
            typename GType::ConstIterator end = graph.end();
            for (auto it=begin; it!=end; ++it) {
                const auto& node = *it;
                // case-by-case on NodePropType
                if constexpr(std::is_void_v<NodePropType>) {
                    os << node << "; ";
                } else if (node_fmt.fmt.has_value()) {
                    os << node << '[' << node_fmt.fmt.value()(graph.node_prop(it));
                    os << "]; ";
                } else if constexpr(detail::is_either_map_v<NodePropType>) {
                    os << node << '[';
                    detail::serialize_map_like(os, graph.node_prop(it));
                    os << "]; ";
                } else if constexpr(detail::is_streamable_v<NodePropType>) {
                    os << node << "[label=\"" << graph.node_prop(it) << "\"]; ";
                }
                ++node_count;
                // start a new line with indent if not the last node in the graph
                if (node_count && !(node_count % max_num_nodes_per_line) && node_count!=graph.size()) {
                    os << '\n';
                    add_indent(os, 1);
                }
            }
            os << '\n';
        }

        const NodeType& get_neighbor_node(typename GType::NeighborsConstIterator n_it) const {
            if constexpr(std::is_void_v<EdgePropType>) {
                return *n_it;
            } else {
                return n_it->first;
            }
        }
        void print_edge(std::ostream& os,
                        const NodeType& curr,
                        const std::string& edge_op,
                        typename GType::NeighborsConstIterator n_it,
                        size_t& edge_count) const {
            os << curr << edge_op << get_neighbor_node(n_it);
            if constexpr(std::is_void_v<EdgePropType>) {
                os << "; ";
            } else if (edge_fmt.fmt.has_value()) {
                os << '[' << edge_fmt.fmt.value()(n_it->second.prop());
                os << "]; ";
            } else if constexpr(detail::is_either_map_v<EdgePropType>){
                os << '[';
                detail::serialize_map_like(os, n_it->second.prop());
                os << "]; ";
            } else if constexpr(detail::is_streamable_v<EdgePropType>) {
                os << "[label=\"" << n_it->second.prop() << "\"]; ";
            }
            ++edge_count;
            if (edge_count && !(edge_count % max_num_edges_per_line) && edge_count!=graph.num_edges()) {
                os << '\n';
                add_indent(os, 1);
            }
        }
        struct hash_node_pair {
            std::size_t operator() (const std::pair<NodeType, NodeType> &p) const {
                std::hash<NodeType> hasher;
                std::size_t seed = hasher(p.first);
                seed ^= hasher(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
                return seed;
            }
        };

        // basically, do a BFS...
        void serialize_edges(std::ostream& os) const {
            if (!graph.size()) {
                return;
            }
            using NSet = std::conditional_t<adj_list_spec==Map::MAP, std::set<NodeType>, std::unordered_set<NodeType>>;
            NSet visited_nodes;
            using ESet = std::conditional_t<adj_list_spec==Map::MAP,
                    std::set<std::pair<NodeType, NodeType>>,
                    std::unordered_set<std::pair<NodeType, NodeType>, hash_node_pair>>;
            ESet visited_edges;
            std::deque<NodeType> queue;
            std::string edge_op = direction==EdgeDirection::UNDIRECTED ? "--": "->";
            typename GType::ConstIterator begin = graph.begin();
            typename GType::ConstIterator end = graph.end();
            add_indent(os, 1);
            size_t edge_count = 0;
            for (auto it=begin; it!=end; ++it) {
                const auto& root = *it;
                if (visited_nodes.count(root)) { continue; }  // skip nodes from visited connected components
                queue.emplace_back(root);
                visited_nodes.insert(root);
                while (!queue.empty()) {
                    NodeType curr = queue.front();
                    queue.pop_front();
                    const auto [n_begin, n_end] = get_neighbors(curr);
                    for (auto n_it = n_begin; n_it != n_end ; ++n_it) {
                        const NodeType& neighbor = get_neighbor_node(n_it);
                        // print edge between curr and neighbor
                        if (!visited_nodes.count(neighbor)) {
                            // unseen neighbor; edge must be new
                            print_edge(os, curr, edge_op, n_it, edge_count);
                            if constexpr(direction==EdgeDirection::UNDIRECTED) {
                                visited_edges.insert({curr, neighbor});
                            }
                            visited_nodes.insert(neighbor);
                            queue.emplace_back(neighbor);
                        } else {
                            // neighbor has already been seen; many possibilities
                            if constexpr(direction==EdgeDirection::DIRECTED) {
                                // always new for a directed graph
                                print_edge(os, curr, edge_op, n_it, edge_count);
                            } else {
                                // for an undirected graph, make sure that this edge (u, v) has not been added by (v, u)
                                // also check for self-loop so that it works for multi-self-loop
                                if (curr==neighbor or !visited_edges.count({neighbor, curr})) {
                                    print_edge(os, curr, edge_op, n_it, edge_count);
                                    visited_edges.insert({curr, neighbor});
                                }
                            }
                        }
                    }
                }
            }
            os << '\n';
        }

    private:  // user-provided prop formatters
        detail::formatter<NodePropType> node_fmt;
        detail::formatter<EdgePropType> edge_fmt;

        const GType& graph;
    public:
        Serializer() = delete;

        template<typename G>
        explicit Serializer(const G& graph): graph{graph} {
            assert(max_num_nodes_per_line > 0 && max_num_edges_per_line > 0);
        }

        void set_max_num_nodes_per_line(size_t num) { max_num_nodes_per_line = num; }
        void unset_max_num_nodes_per_line() { max_num_nodes_per_line = (size_t)-1; }
        void set_max_num_edges_per_line(size_t num) { max_num_edges_per_line = num; }
        void unset_max_num_edges_per_line() { max_num_edges_per_line = (size_t)-1; }

        // F has type NodePropType -> string
        template<typename F>
        void register_node_formatter(F func) {
            static_assert(not std::is_void_v<NodePropType>, "no node prop needed at all");
            static_assert(std::is_assignable_v<decltype(node_fmt.fmt), F>);
            node_fmt.fmt = func;
        }
        void delete_node_formatter() {
            static_assert(not std::is_void_v<NodePropType>, "no node prop needed at all");
            node_fmt.fmt = std::nullopt;
        }

        template<typename F>
        void register_edge_formatter(F func) {
            static_assert(not std::is_void_v<EdgePropType>, "no edge prop needed at all");
            static_assert(std::is_assignable_v<decltype(edge_fmt.fmt), F>);
            edge_fmt.fmt = func;
        }
        void delete_edge_formatter() {
            static_assert(not std::is_void_v<EdgePropType>, "no edge prop needed at all");
            edge_fmt.fmt = std::nullopt;
        }

        // serialize to the dot language
        void serialize_to_dot(std::ostream& os) const {
            message<true>();
            message<false>();
            if constexpr(multi_edge==MultiEdge::DISALLOWED) {
                os << "strict ";
            }
            if constexpr(direction==EdgeDirection::DIRECTED) {
                os << "digraph {\n";
            } else {
                static_assert(direction==EdgeDirection::UNDIRECTED);
                os << "graph {\n";
            }
            serialize_nodes(os);
            serialize_edges(os);
            os << '}' << std::endl;
        }
    };
    // deduction guide for serializer
    template<typename G>
    explicit Serializer(const G& graph) -> Serializer<
            typename G::node_type, typename G::node_prop_type, typename G::edge_prop_type,
            G::DIRECTION, G::MULTI_EDGE, G::SELF_LOOP,
            G::ADJ_LIST_SPEC, G::NEIGHBORS_CONTAINER_SPEC
    >;
}

#endif //GRAPHLITE_SERIALIZE_H
