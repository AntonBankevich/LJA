//
// Created by Tatiana Dvorkina on 04.05.2022.
//

#ifndef DR_GRAPHSIMPLIFICATOR_H
#define DR_GRAPHSIMPLIFICATOR_H

#include <lja/multi_graph.hpp>

namespace nano {
    class GraphSimplificator {
    public:
        GraphSimplificator(const std::unordered_map<int, std::unordered_map<int, std::pair<int, std::vector<int>>>> &sgraph):
        sgraph_(sgraph) {}

        void Simplify(multigraph::MultiGraph &mg, const std::experimental::filesystem::path &dir);

        void ResolveWithMajor(multigraph::MultiGraph &mg);

    private:
        std::unordered_map<int, int> chains_;
        std::unordered_set<int> inside_edges_;
        const std::unordered_map<int, std::unordered_map<int, std::pair<int, std::vector<int>>>> sgraph_;
    };
}

#endif //DR_GRAPHSIMPLIFICATOR_H
