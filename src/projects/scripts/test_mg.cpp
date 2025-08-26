#include "dbg/multi_graph.hpp"
#include "common/id_index.hpp"
#include "sequences/seqio.hpp"
#include <experimental/filesystem>
#include <array>
#include <vector>
#include <unordered_map>
#include <assembly_graph/visualization.hpp>

int main(int argc, char **argv) {
    multigraph::MultiGraph mg;
    //multigraph::MultiGraph mg = mmg.DBG();

    multigraph::Vertex &v = mg.addVertex(Sequence("AT"), 23);
    auto &e1 = mg.addEdge(v, v, Sequence("ATATAT"));
    auto &e2 = mg.addEdge(v, v, Sequence("ATCAT"));
    IdIndex<multigraph::Vertex> index(mg.vertices().begin(), mg.vertices().end());
    std::cerr <<v.getId() << " "<< index.getById(23).getId() <<std::endl;
    std::cerr << index.getById(23).getId() <<std::endl;
    for (multigraph::Edge &edge : index.getById(23)){
        std::cerr<<"EEEdge " <<edge.getId() <<std::endl;
    }
    ag::Printer().printGFA("bd.gfa", mg);
    IdIndex<multigraph::Edge> eindex(mg.edges().begin(), mg.edges().end());
    std::cout << mg.size() << " v/e " << mg.edgeCount() << std::endl;
    std::cout.flush();
    std::cout <<"in/out degs" << index.getById(23).inDeg() << " " << index.getById(23).outDeg() <<std::endl;
    std::cout << eindex.getById(e1.getInnerId()).getStart().getId() << " " << eindex.getById(e1.getInnerId()).getFinish().getId() << std::endl;
    multigraph::MultiGraphHelper::printExtractedContigs(mg, "edges.fasta", false);
    mg.removeEdge(eindex.getById(e1.getInnerId()));
    ag::Printer().printGFA("ad.gfa", mg);
}
