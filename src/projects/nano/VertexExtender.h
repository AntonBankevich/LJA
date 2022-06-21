//
// Created by Tatiana Dvorkina on 21.06.2022.
//

#ifndef DR_VERTEXEXTENDER_H
#define DR_VERTEXEXTENDER_H

#include <lja/multi_graph.hpp>

namespace nano {
    class VertexExtender {
    public:
        std::unordered_map<int, int> ExtendVertices(multigraph::MultiGraph &mg);

    private:
        int IsDimerStructure(const int v_id, multigraph::MultiGraph &mg);
        bool IsDimerString(const std::string &str, int sz);
        bool IsNearlyDimerString(const std::string &str);
    };
}


#endif //DR_VERTEXEXTENDER_H
