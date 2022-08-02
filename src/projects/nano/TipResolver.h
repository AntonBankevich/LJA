//
// Created by Tatiana Dvorkina on 05.07.2022.
//

#ifndef DR_TIPRESOLVER_H
#define DR_TIPRESOLVER_H

#include <lja/multi_graph.hpp>
#include <nano/GraphContig.hpp>

namespace nano {
class TipResolver {
public:
    explicit TipResolver(multigraph::MultiGraph &mg) : mg_(mg) {
        ExtractTipVertices();
    };

    void LoadPaths(std::unordered_map<std::string, std::vector<GraphContig>> &alignments);

    void ResolveWithPaths();

private:

    void ExtractTipVertices();

    bool AlignedOnTipVertex(const std::vector<GraphContig> &r_alignments);

    std::unordered_set<std::string> IdentifyConnectedVertices(const std::vector<GraphContig> &r_alignments);

    std::unordered_map<std::string, std::vector<GraphContig>> tipReads_;
    std::set<int> tipvertices_;
    multigraph::MultiGraph &mg_;
};

}


#endif //DR_TIPRESOLVER_H
