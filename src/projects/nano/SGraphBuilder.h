//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#ifndef DR_SGRAPHBUILDER_H
#define DR_SGRAPHBUILDER_H

#include <lja/multi_graph.hpp>
#include <nano/GraphContig.hpp>

namespace nano {

    class SGraphBuilder {
    public:
        SGraphBuilder(const multigraph::MultiGraph &mg,
                      const std::experimental::filesystem::path &unique_edges,
                      bool by_vertex = true)
                      :mg_(mg), by_vertex_(by_vertex) {
            LoadUEdges(unique_edges);
        }

        void LoadAlignments(const std::unordered_map<std::string, nano::GraphContig> &alignments);
        void PrintSgraph();

    private:
        void LoadUEdges(const std::experimental::filesystem::path &unique_edges);
        void AddEdge(const std::string &prev_edge_id, const std::string &edge_id);
        std::pair<std::string, std::string> ResolveByVertex(const std::vector<std::string> &subpath, const Contig &read_str);
        const multigraph::Edge* GetEdgeByNuc(const char nuc, const multigraph::Edge *edge, bool is_inEdge = true);

        const multigraph::MultiGraph &mg_;
        std::unordered_set<std::string> uedges_;
        std::unordered_map<std::string, std::unordered_map<std::string, int>> sgraph_;
        bool by_vertex_;
    };
}

#endif //DR_SGRAPHBUILDER_H
