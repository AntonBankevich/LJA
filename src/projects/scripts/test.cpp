#include <array>
#include <vector>
#include <unordered_map>
#include <experimental/filesystem>
#include "dbg/multi_graph.hpp"
#include "../lja/subdataset_processing.hpp"
#include "sequences/seqio.hpp"

int main(int argc, char **argv) {
    multigraph::MultiGraph mg;
    //multigraph::MultiGraph mg = mmg.DBG();

    multigraph::Vertex &v = mg.addVertex(Sequence("AT"), 23, "vertex");
    std::cerr <<v.getId() << " "<< mg.getVertexById(23) <<std::endl;
    auto &e1 = mg.addEdge(v, v, Sequence("ATATAT"), 239, "e1");
    auto &e2 = mg.addEdge(v, v, Sequence("ATCAT"), 240, "e2");
    std::cerr << mg.getVertexById(23) <<std::endl;
    for (multigraph::Edge &edge : *mg.getVertexById(23)){
        std::cerr<<"EEEdge " <<edge.getId() <<std::endl;
    }
    multigraph::MultiGraphHelper::printEdgeGFA(mg, "bd.gfa");
    multigraph::MultiGraphHelper::printVertexGFA(mg, "vertex.gfa");
    std::cout<<mg.size() << " v/e " << mg.edgeNumber() <<std::endl;
    std::cout.flush();
    std::cout <<"in/out degs" << mg.getVertexById(23)->inDeg() << " " <<mg.getVertexById(23)->outDeg() <<std::endl;
    std::cout << mg.getEdgeById(239)->start().getId() << " " << mg.getEdgeById(239)->end().getId() <<std::endl;
    multigraph::MultiGraphHelper::printExtractedContigs(mg, "edges.fasta", false);
    mg.internalRemoveEdge(*mg.getEdgeById(239));
    multigraph::MultiGraphHelper::printEdgeGFA(mg, "ad.gfa");
}
