#include "dbg/multi_graph.hpp"
#include "dbg/id_index.hpp"
#include "sequences/seqio.hpp"
#include <experimental/filesystem>
#include <array>
#include <vector>
#include <unordered_map>

int main(int argc, char **argv) {
    multigraph::MultiGraph mg;
    //multigraph::MultiGraph mg = mmg.DBG();

    multigraph::MGVertex &v = mg.addVertex(Sequence("AT"), multigraph::MGVertexData("getVertex"), 23);
    auto &e1 = v.addEdge(v, Sequence("ATATAT"), multigraph::MGEdgeData("e1"));
    auto &e2 = v.addEdge(v, Sequence("ATCAT"), multigraph::MGEdgeData("e2"));
    IdIndex<multigraph::Vertex> index(mg.vertices().begin(), mg.vertices().end());
    std::cerr <<v.getId() << " "<< index.getById(23).getId() <<std::endl;
    std::cerr << index.getById(23).getId() <<std::endl;
    for (multigraph::MGEdge &edge : index.getById(23)){
        std::cerr<<"EEEdge " <<edge.getId() <<std::endl;
    }
    multigraph::MultiGraphHelper::printEdgeGFA(mg, "bd.gfa");
    multigraph::MultiGraphHelper::printVertexGFA(mg, "getVertex.gfa");
    IdIndex<multigraph::Edge> eindex(mg.edges().begin(), mg.edges().end());
    std::cout << mg.size() << " v/e " << mg.edgeCount() << std::endl;
    std::cout.flush();
    std::cout <<"in/out degs" << index.getById(23).inDeg() << " " << index.getById(23).outDeg() <<std::endl;
    std::cout << eindex.getById(e1.getInnerId()).getStart().getId() << " " << eindex.getById(e1.getInnerId()).getFinish().getId() << std::endl;
    multigraph::MultiGraphHelper::printExtractedContigs(mg, "edges.fasta", false);
    eindex.getById(e1.getInnerId()).getStart().removeEdge(eindex.getById(e1.getInnerId()));
    multigraph::MultiGraphHelper::printEdgeGFA(mg, "ad.gfa");
}
