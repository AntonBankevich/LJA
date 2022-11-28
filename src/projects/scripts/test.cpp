#include <array>
#include <vector>
#include <unordered_map>
#include "../lja/multi_graph.hpp"
#include "../lja/subdataset_processing.hpp"
#include "sequences/seqio.hpp"

int main(int argc, char **argv) {
    multigraph::MultiGraph mg;
    //multigraph::MultiGraph mg = mmg.DBG();

    multigraph::Vertex &v = mg.addVertex(Sequence("AT"), 23, "vertex");
    std::cerr <<&v << " "<< &mg.vertices[23] <<std::endl;
    auto &e1 = mg.addEdge(v, v, Sequence("ATATAT"), 239, "e1");
    auto &e2 = mg.addEdge(v, v, Sequence("ATCAT"), 240, "e2");
    std::cerr << &mg.vertices[23] <<std::endl;
    for (size_t i = 0; i < mg.vertices[23].outgoing.size() ; i++ ){
        std::cerr<<"EEEdge " <<mg.vertices[23].outgoing[i]->getId() <<std::endl;
    }
    mg.printEdgeGFA("bd.gfa");
    mg.printVertexGFA("vertex.gfa");
    std::cout<<mg.vertices.size() << " v/e " << mg.edges.size() <<std::endl;
    std::cout.flush();
    std::cout <<"in/out degs" << mg.vertices[23].inDeg() << " " <<mg.vertices[23].outDeg() <<std::endl;
    std::cout << mg.edges[239].start->id << " " << mg.edges[239].end->id <<std::endl;
    mg.printEdges("edges.fasta", false);
    mg.deleteEdgeById(239);
    mg.printEdgeGFA("ad.gfa");

}
