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
                      const std::unordered_map<int, int> &new_edges_map,
                      bool by_vertex = true)
                      :mg_(mg), by_vertex_(by_vertex) {
            LoadUEdges(unique_edges, new_edges_map);
        }

        void LoadAlignments(const std::unordered_map<std::string, nano::GraphContig> &alignments,
                            const size_t threads);

        void LoadSGraphEdges(const std::experimental::filesystem::path &sgraph_filename);

        std::unordered_map<int, std::unordered_map<int, std::pair<int, std::vector<int>>>>& GetSGraph() {
            return sgraph_;
        };

        void PrintSGraph();

        void SaveSGraph(const std::experimental::filesystem::path &sgraph_filename);

    private:
        void LoadUEdges(const std::experimental::filesystem::path &unique_edges,
                        const std::unordered_map<int, int> &new_edges_map);
        void AddEdge(const std::string &prev_edge_id, const std::string &edge_id, const std::vector<std::string> &subpath);
        std::pair<std::string, std::string> ResolveByVertex(const std::vector<std::string> &subpath, const Contig &read_str);
        const multigraph::Edge* GetEdgeByNuc(const char nuc, const multigraph::Edge *edge, bool is_inEdge = true);

        const multigraph::MultiGraph &mg_;
        std::unordered_set<std::string> uedges_;
        std::unordered_map<int, std::unordered_map<int, std::pair<int, std::vector<int>>>> sgraph_;
        bool by_vertex_;
    };
}

#endif //DR_SGRAPHBUILDER_H
