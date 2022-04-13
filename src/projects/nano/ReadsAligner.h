//
// Created by Tatiana Dvorkina on 12.04.2022.
//

#ifndef DR_READSALIGNER_H
#define DR_READSALIGNER_H

#include <nano/GraphContig.hpp>
#include <lja/multi_graph.hpp>

class ReadsAlignerGA {
public:
    ReadsAlignerGA(const multigraph::MultiGraph &mg): mg_(mg) {}

    std::vector<GraphContig> Align(const std::vector<StringContig> &sequences,
                                   );

private:
    const multigraph::MultiGraph &mg_;

};


#endif //DR_READSALIGNER_H
