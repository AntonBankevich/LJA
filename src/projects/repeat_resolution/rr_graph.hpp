//
// Created by Andrey Bzikadze on 10/14/21.
//

#pragma once

#include <functional>
#include <graphlite/graphlite.hpp>
#include <graphlite/serialize.hpp>

namespace repeat_resolution {
struct RRVertexType {
  uint64_t len {0};
  bool frozen {false};
};
//
//    inline bool operator==(const RRVertexType & lhs, const RRVertexType & rhs)
//    { return lhs.id == rhs.id; } inline bool operator!=(const RRVertexType &
//    lhs, const RRVertexType & rhs) { return !operator==(lhs, rhs); } inline
//    bool operator< (const RRVertexType & lhs, const RRVertexType & rhs) {
//    return lhs.id < rhs.id; } inline bool operator> (const RRVertexType & lhs,
//    const RRVertexType & rhs) { return  operator< (rhs, lhs); } inline bool
//    operator<=(const RRVertexType & lhs, const RRVertexType & rhs) { return
//    !operator> (lhs, rhs); } inline bool operator>=(const RRVertexType & lhs,
//    const RRVertexType & rhs) { return !operator< (lhs, rhs); }
//
    std::ostream& operator<<(std::ostream & os, const RRVertexType &
    vertex_property) {
        size_t new_len = 7;
        if (vertex_property.id.front() == '-') {
            ++new_len;
        }
        os << vertex_property.id.substr(0, new_len);
        return os;
    }
} // namespace repeat_resolution
//
// namespace std {
//    template<> struct hash<repeat_resolution::RRVertexType> {
//        size_t operator()(const repeat_resolution::RRVertexType & s) const
//        noexcept {
//            return hash<string>()(s.id);
//        }
//    };
//} // End namespace std

namespace repeat_resolution {

using EdgeIndexType = uint64_t;

//    struct RREdgeProperty {
//        std::vector<uint64_t> dbg_edge_ids;
//        RREdgeProperty& operator+=(const RREdgeProperty& rhs) {
//            dbg_edge_ids.insert(dbg_edge_ids.end(), rhs.dbg_edge_ids.cbegin(),
//            rhs.dbg_edge_ids.cend()); return *this;
//        }
//    };
//    inline RREdgeProperty operator+(RREdgeProperty lhs, const RREdgeProperty
//    &rhs) {
//        lhs += rhs;
//        return lhs;
//    }
//
//    std::ostream& operator<<(std::ostream & os, const RREdgeProperty &
//    edge_property) {
//        for (const uint64_t id : edge_property.dbg_edge_ids) {
//            os << id << " ";
//        }
//        os << "\n";
//        return os;
//    }
//
//    class RRGraph :
//        public graph_lite::Graph<
//            /*typename NodeType=*/RRVertexType,
//            /*typename NodePropType=*/void,
//            /*typename EdgePropType=*/RREdgeProperty,
//            /*EdgeDirection direction=*/graph_lite::EdgeDirection::DIRECTED,
//            /*MultiEdge multi_edge=*/graph_lite::MultiEdge::ALLOWED,
//            /*SelfLoop self_loop=*/graph_lite::SelfLoop::ALLOWED,
//            /*Map adj_list_spec=*/graph_lite::Map::UNORDERED_MAP,
//            /*Container
//            neighbors_container_spec=*/graph_lite::Container::MULTISET> {
//        bool node_has_loop(const node_type & node) const {
//            auto [begin, end] = in_neighbors(node);
//            for (auto it = begin; it != end; ++it) {
//                const node_type & neighbor = it->first;
//                if (neighbor == node) {
//                    return true;
//                }
//            }
//            return false;
//        }
//    public:
//        explicit RRGraph(dbg::SparseDBG & dbg) {
//            uint64_t cnt = 0;
//            for (auto it = dbg.edges().begin(); it != dbg.edges().end(); ++it)
//            {
//                const Edge & edge = *it;
//                const std::string start_id = edge.start()->getId();
//                const std::string end_id = edge.end()->getId();
//                add_nodes(RRVertexType{start_id});
//                add_nodes(RRVertexType{end_id});
//                RREdgeProperty edge_property { {cnt} };
//                add_edge_with_prop(
//                        RRVertexType{start_id},
//                        RRVertexType{end_id},
//                        edge_property);
//                ++cnt;
//            }
//        }
//
//        void serialize_to_dot(const std::experimental::filesystem::path &
//        path) const {
//            graph_lite::Serializer serializer(*this);
//            std::ofstream dot_os(path);
//            serializer.serialize_to_dot(dot_os);
//        }
//
//        void resolve_graph() {
//            std::vector<node_type> nodes_to_remove;
//            for (auto it = begin(); it != end(); ++it) {
//                node_type node { *it };
//                int indegree = count_in_neighbors(node);
//                int outdegree = count_out_neighbors(node);
//                VERIFY(indegree != 1 or outdegree != 1); // a node cannot be
//                on a non-branching path if (node_has_loop(node)) {
//                    // TODO. Loops need special care. For now, skip
//                }
//                if (indegree == 0 or outdegree == 0) {
//                    /* isolated vertex or a (multi-)tip should be skipped
//                     * multi-tip (indegree = 0, outdegree > 1 or vice versa)
//                     technically can be resolved
//                     * but that will not increase the length of contigs and
//                     thus is meaningless.
//                     * It will also invalidate iterator unless we switch to
//                     ordered containers
//                     */
//                } else if (indegree == 1) {
//                    auto [ibegin, iend] = in_neighbors(it);
//                    auto & [inode, iedge] = *ibegin;
//                    auto [obegin, oend] = out_neighbors(it);
//                    for (auto oit = obegin; oit != oend; ++oit) {
//                        // inode --[iedge]> node --[oedge]> onode
//                        auto &[onode, oedge] = *oit;
//                        RREdgeProperty new_prop = iedge.prop() + oedge.prop();
//                        add_edge_with_prop(inode, onode, new_prop);
//                    }
//                    nodes_to_remove.emplace_back(std::move(node));
//                } else if (outdegree == 1) {
//                    auto [ibegin, iend] = in_neighbors(it);
//                    auto [obegin, oend] = out_neighbors(it);
//                    auto & [onode, oedge] = *obegin;
//                    for (auto iit = ibegin; iit != iend; ++iit) {
//                        // inode --[iedge]> node --[oedge]> onode
//                        auto &[inode, iedge] = *iit;
//                        RREdgeProperty new_prop = iedge.prop() + oedge.prop();
//                        add_edge_with_prop(inode, onode, new_prop);
//                    }
//                    nodes_to_remove.emplace_back(std::move(node));
//                }
//            }
//            for (const node_type & node : nodes_to_remove) {
//                remove_nodes(node);
//            }
//        }
//    };

} // End namespace repeat_resolution