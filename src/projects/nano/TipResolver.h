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

    std::unordered_map<int, int> ResolveWithPaths();

private:

    void ExtractTipVertices();

    bool AlignedOnTipVertex(const std::vector<GraphContig> &r_alignments);

    std::pair<int, int> GetBestEdgePair(const std::string &key, int v1, int v2, bool rc);

    std::unordered_map<int, std::unordered_map<int,std::vector<std::string> > > GetTransitions
                                        (const std::vector<std::string> &count, int v1, int v2);

    std::string GetConsensus(std::vector<std::string> &str_lst
                            , multigraph::Edge *e_in
                            , multigraph::Edge *e_out);

    int AlignRead2Tip(std::string& read_str, std::string& edge_str);

    std::unordered_set<std::string> IdentifyConnectedVertices(const std::vector<GraphContig> &r_alignments);

    std::unordered_map<std::string, std::vector<GraphContig>> tipReads_;
    std::set<int> tipvertices_;
    multigraph::MultiGraph &mg_;
};

}


#endif //DR_TIPRESOLVER_H
