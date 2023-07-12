#include <common/cl_parser.hpp>
#include <dbg/multi_graph.hpp>
#include "graph_polishing.hpp"
using namespace multigraph;

int main(int argc, char **argv) {
    GraphPolishingPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram polishing("GraphPolishing", std::move(phase), std::move(parser), "Starting graph polishing procedure", "Finished graph polishing procedure");
    polishing.run(oneline::initialize<std::string, char*>(argv, argv + argc));
}

std::unordered_map<std::string, std::experimental::filesystem::path>
RunGraphPolishing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                  const io::Library &graph_gfa, const io::Library &corrected_edges, bool debug) {
    logger.info() << "Loading graph" << std::endl;
    MultiGraph mg = MultiGraphHelper::LoadGFA(graph_gfa.front(), true);
    mg.checkConsistency();
    logger.info() << "Preparing graph" << std::endl;
    mg = MultiGraphHelper::TransformToEdgeGraph(mg, 5001);
    mg.checkConsistency();
    logger.info() << "Linking positions" << std::endl;
    DisjointSet<EdgePosition> linked_positions;
    std::vector<EdgePosition> queue;
    std::unordered_set<EdgePosition> visited;
    for(Edge &edge : mg.edges()) {
        queue.emplace_back(edge, 0);
        queue.emplace_back(edge, edge.size());
    }
    while(!queue.empty()) {
        EdgePosition ep = queue.back();
        queue.pop_back();
        if(visited.find(ep) != visited.end())
            continue;
        visited.insert(ep);
        if(ep.getPos() <= ep.contig().start().size()) {
            Vertex &start = ep.contig().start();
            for(Edge &rcedge2: start.rc()) {
                Edge &edge2 = rcedge2.rc();
                EdgePosition other(edge2, edge2.size() - start.size() + ep.getPos());
                if(visited.find(other) == visited.end()) {
                    if(ep.getPos() < start.size())
                        linked_positions.link(ep, other);
                    queue.emplace_back(other);
                }
            }
            for(Edge &edge2: start) {
                EdgePosition other(edge2, ep.getPos());
                if(visited.find(other) == visited.end()) {
                    if(ep.getPos() < start.size())
                        linked_positions.link(ep, other);
                    queue.emplace_back(other);
                }
            }
        }
        if(ep.getPos() >= ep.contig().size() - ep.contig().end().size()) {
            Vertex &end = ep.contig().end();
            for(Edge &edge2 : end) {
                EdgePosition other(edge2, ep.getPos() - ep.contig().size() + end.size());
                if(visited.find(other) == visited.end()) {
                    if(ep.getPos() < ep.contig().size())
                        linked_positions.link(ep, other);
                    queue.emplace_back(other);
                }
            }
            for(Edge &rcedge2 : end.rc()) {
                Edge &edge2 = rcedge2.rc();
                EdgePosition other(edge2, edge2.size() - ep.contig().size() + ep.getPos());
                if(visited.find(other) == visited.end()) {
                    if(ep.getPos() < ep.contig().size())
                        linked_positions.link(ep, other);
                    queue.emplace_back(other);
                }
            }
        }
    }

    logger.info() << "Preparing output" << std::endl;
    std::vector<EdgePosition> positions(visited.begin(), visited.end());
    std::unordered_map<EdgePosition, size_t> len;
    std::sort(positions.begin(), positions.end());
    std::vector<size_t> res;
    EdgePosition prev;
    for(EdgePosition pos: positions) {
        if(pos != positions.front() && pos.contig() == prev.contig()) {
            size_t tmp = pos.getPos() - prev.getPos();
            VERIFY(tmp > 0);
            size_t &tmp1 = len[prev];
            VERIFY(tmp1 == 0|| tmp1 == tmp);
            tmp1 = tmp;
            res.back()++;
        } else {
            res.emplace_back(1);
        }
        prev = pos;
    }
    std::ofstream os1;
    os1.open("edges.txt");
    std::sort(res.begin(), res.end());
    size_t cnt = 1;
    for(size_t i =1; i < res.size(); i++) {
        if(res[i] == res[i-1]) {
            cnt++;
        } else {
            os1 << res[i-1] << " " << cnt  << std::endl;
            cnt = 1;
        }
    }
    os1 << res.back() << " " << cnt << std::endl;
    os1.close();
    std::ofstream os2;
    os2.open("segments.txt");
    for(auto &it : linked_positions.nontrivialSubsets()) {
        os2 << it.second.size() << " " << len[it.second.front()] << std::endl;
        if(len[it.second.front()] > 100000) {
            std::cout << it.second.size() << " " << len[it.second.front()] << std::endl;
            std::cout << it.second << std::endl;
        }
    }
    os2.close();
    return {};
}
