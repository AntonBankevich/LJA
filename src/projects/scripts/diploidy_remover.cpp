#include "dbg/graph_algorithms.hpp"
#include "dbg/multi_graph.hpp"
#include <common/dir_utils.hpp>
#include "common/disjoint_sets.hpp"
#include <experimental/filesystem>
#include <unordered_map>
#include <common/pipeline_tools.hpp>

using namespace multigraph;



class ConsensusPhase : public Stage {
public:
    ConsensusPhase() : Stage(AlgorithmParameters({}, {}, ""), {"graph"}, {"graph", "assembly"}) {
    }

private:
    void addGoodVertex(MGVertex &v, std::unordered_set<VertexId> &good, std::vector<EdgeId> &candidates) {
        if(good.find(v.getId()) != good.end())
            return;
        std::cout << "Good vertex " << v.getId() << std::endl;
        good.insert(v.getId());
        good.insert(v.rc().getId());
        for(MGEdge &e : v) {
            if(good.find(e.getFinish().getId()) == good.end())
                candidates.emplace_back(e.getId());
        }
        if(v != v.rc())
            for(MGEdge &e : v.rc()) {
                if(good.find(e.getFinish().getId()) == good.end())
                    candidates.emplace_back(e.getId());
            }
    }

    std::vector<EdgeId> component(MGEdge &initial) {
        std::unordered_set<VertexId> visited;
        std::unordered_set<EdgeId> res;
        std::vector<VertexId> queue;
        queue.emplace_back(initial.getStart().getId());
        queue.emplace_back(initial.getFinish().rc().getId());
        visited.insert(initial.getStart().rc().getId());
        visited.insert(initial.getFinish().getId());
        while(!queue.empty()) {
            MGVertex &v = *queue.back();
            queue.pop_back();
            if(visited.find(v.getId()) != visited.end())
                continue;
            visited.insert(v.getId());
            for(MGEdge &e: v) {
                if(visited.find(e.getFinish().getId()) == visited.end())
                    queue.emplace_back(e.getFinish().getId());
                if(visited.find(e.getFinish().rc().getId()) == visited.end())
                    queue.emplace_back(e.getFinish().rc().getId());
            }
        }
        for(VertexId v: visited) {
            if(v != initial.getStart().rc().getId() && v != initial.getFinish().getId())
                for(MGEdge &e: *v) {
                    if(e != initial && e != initial.rc()) {
                        if (e.getSeq() <= e.rc().getSeq())
                            res.insert(e.getId());
                        else
                            res.insert(e.rc().getId());
                    }
                }
        }
        return {res.begin(), res.end()};
    }

    void CollapseSimpleBulges(logging::Logger &logger, size_t threads, MultiGraph &mg) {
        std::unordered_set<VertexId> good;
        std::vector<EdgeId> candidates;
        std::unordered_set<ConstEdgeId> to_remove;
        DisjointSet<VertexId> paths;
        for(MGVertex &v: mg.vertices()) {
            paths.link(v.getId(), v.rc().getId());
            if(v.outDeg() > 0 && v.outDeg() <= 2) {
                bool ok = true;
                for(MGEdge &e : v) {
                    if(e.getFinish() != v.front().getFinish()) {
                        ok = false;
                        break;
                    }
                }
                if(ok && v.outDeg() == v.front().getFinish().inDeg())
                    paths.link(v.getId(), v.front().getFinish().getId());
            }
        }
        std::unordered_map<VertexId, std::vector<VertexId>> split = paths.nontrivialSubsets();
        for(auto &it : split) {
            size_t len = 0;
            for(VertexId v : it.second) {
                if(v->outDeg() > 0 && paths.get(v->begin()->getFinish().getId()) == it.first) {
                    len += 2 * std::min(v->front().fullSize(), v->back().fullSize()) - v->front().getFinish().size() - v->size();
                }
            }
            len /= 2;
            if(len > 200000) {
                for(VertexId v : it.second) {
                    addGoodVertex(*v, good, candidates);
                }
            }
        }
        while(!candidates.empty()) {
            EdgeId e = candidates.back();
            candidates.pop_back();
            if(to_remove.find(e) !=to_remove.end())
                continue;
            if(e->getStart().outDeg() == 1 && e->getFinish().inDeg() == 1) {
                addGoodVertex(e->getFinish(), good, candidates);
                continue;
            }
            if(e->getStart().outDeg() > 2 || e->getFinish().inDeg() > 2) {
                continue;
            }
            std::cout << "Checking edge " << e->getId() << std::endl;
            std::vector<EdgeId> comp = component(*e);
            size_t size = 0;
            for(EdgeId e1: comp) {
                size+= 2 * e1->fullSize() - e1->getStart().size() - e1->getFinish().size();
            }
            size_t elen = 2* e->fullSize() - e->getStart().size() - e->getFinish().size();
            std::cout << elen << " " << comp.size() << " " << size << std::endl;
            if((comp.size() == 1 && e->getStart().outDeg() == 2 && e->getFinish().inDeg() == 2 && size <= elen) || size < elen * 1.5 + 100000) {
                for(EdgeId e1 : comp) {
                    std::cout << "Deleting edges " << e1->getId() << " and " << e1->rc().getId() << " based on analysis of edge " << e->getId() << std::endl;
                    to_remove.emplace(e1->getId());
                    to_remove.emplace(e1->rc().getId());
                }
                addGoodVertex(e->getFinish(), good, candidates);
            }
        }
        mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()});
        std::cout << "cleaned " << mg.size() << " " << mg.edgeCount() << std::endl;
        ag::MergeAll(logger, threads, mg);
        std::cout << "merged " << mg.size() << " " << mg.edgeCount() << std::endl;
    }

    void ChooseShortcuts(logging::Logger &logger, size_t threads, MultiGraph &mg) {
        DisjointSet<VertexId> small_components;
        for(MGEdge &edge : mg.edges()) {
            if(edge.fullSize() < 1000000)
                small_components.link(edge.getStart().getId(), edge.getFinish().getId());
        }
        std::unordered_set<VertexId> all_vert;
        for(MGVertex &vertex : mg.vertices()) {
            all_vert.insert(vertex.getId());
        }
        auto cmap = small_components.subsets(all_vert);
        std::unordered_set<EdgeId> to_remove;
        for(auto & it:cmap) {
            std::unordered_set<VertexId> cvert(it.second.begin(), it.second.end());
            if(cvert.size() > 6)
                continue;
            bool good = true;
            VertexId start;
            VertexId end;
            std::unordered_set<EdgeId> edges;
            for(VertexId v : cvert) {
                for(MGEdge &e: *v) {
                    if(cvert.find(e.getFinish().getId()) == cvert.end()) {
                        if(end.valid()) {
                            good = false;
                        }
                        end = v;
                    } else {
                        edges.insert(e.getId());
                    }
                }
                for(MGEdge &e: v->rc()) {
                    if(cvert.find(e.getFinish().rc().getId()) == cvert.end()) {
                        if(start.valid()) {
                            good = false;
                        }
                        start = v;
                    }
                }
            }
            std::cout << "Component " << cvert.size() << " " << edges.size() << std::endl;
            if(edges.empty() || !start.valid() || !end.valid() || !good) {
                continue;
            }
            std::cout << "Processing component" << std::endl;
            for(EdgeId e : edges) {
                std::cout << e->getId() << " ";
            }
            std::cout << std::endl;
            EdgeId shortcut;
            for(EdgeId e : edges) {
                if(e->getStart().getId() == start && e->getFinish().getId() == end && (!shortcut.valid() ||
                                                                                       shortcut->fullSize() <
                                                                                       e->fullSize())) {
                    shortcut = e;
                }
            }
            size_t size = 0;
            for(EdgeId e : edges) {
                size += e->fullSize() * 2 - e->getStart().size() - e->getFinish().size();
            }
            if(size <= 1000000 && (to_remove.find(shortcut) == to_remove.end()) && (shortcut.valid() || start == end)) {
                if(shortcut.valid())
                    std::cout << "Found shortcut " << shortcut->getId() << std::endl;
                else
                    std::cout << "Found tip " << std::endl;
                for(EdgeId e : edges) {
                    if(e != shortcut && e->rc().getId() != shortcut) {
                        to_remove.insert(e);
                        to_remove.insert(e->rc().getId());
                    }
                }
            }
        }
        mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()});
        std::cout << "cleaned " << mg.size() << " " << mg.edgeCount() << std::endl;
        ag::MergeAll<MGTraits>(logger, threads, mg);
//    MultiGraphHelper::MergeAllPaths(mg, true);
        std::cout << "merged " << mg.size() << " " << mg.edgeCount() << std::endl;
    }

    void PassAcyclic(logging::Logger &logger, size_t threads, MultiGraph &mg) {
        DisjointSet<VertexId> small_components;
        for(MGEdge &edge : mg.edges()) {
            if(edge.fullSize() < 1000000)
                small_components.link(edge.getStart().getId(), edge.getFinish().getId());
        }
        std::unordered_set<VertexId> all_vert;
        for(MGVertex &vertex : mg.vertices()) {
            all_vert.insert(vertex.getId());
        }
        auto cmap = small_components.subsets(all_vert);
        std::unordered_set<ConstEdgeId> to_remove;
        for(auto & it:cmap) {
            std::unordered_set<VertexId> cvert(it.second.begin(), it.second.end());
            if(cvert.size() > 6)
                continue;
            bool good = true;
            VertexId start;
            VertexId end;
            std::unordered_set<EdgeId> edges;
            for(VertexId v : cvert) {
                for(MGEdge &e: *v) {
                    if(cvert.find(e.getFinish().getId()) == cvert.end()) {
                        if(end.valid()) {
                            good = false;
                        }
                        end = v;
                    } else {
                        edges.emplace(e.getId());
                    }
                }
                for(MGEdge &e: v->rc()) {
                    if(cvert.find(e.getFinish().rc().getId()) == cvert.end()) {
                        if(start.valid()) {
                            good = false;
                        }
                        start = v;
                    }
                }
            }
            std::cout << "Component " << cvert.size() << " " << edges.size() << std::endl;
            if(edges.empty() || !start.valid() || !end.valid() || start == end || !good) {
                continue;
            }
            std::cout << "Component " << cvert.size() << " " << edges.size() << std::endl;
            for(EdgeId e : edges) {
                std::cout << e->getId() << " ";
            }
            std::cout << std::endl;
            size_t size = 0;
            for(EdgeId e : edges) {
                size += e->fullSize() * 2 - e->getStart().size() - e->getFinish().size();
            }
            if(size > 1000000)
                continue;
            for(EdgeId e: edges) {
                if(to_remove.find(ConstEdgeId(e)) != to_remove.end())
                    continue;
                std::vector<EdgeId> path = {e};
                while(path.size() < 10 && path.back()->getFinish().outDeg() == 1 && path.back()->getFinish().getId() != end) {
                    path.push_back(path.back()->getFinish().front().getId());
                }
                while(path.size() < 10 && path.front()->getStart().inDeg() == 1 && path.front()->getStart().getId() != start) {
                    path.insert(path.begin(), path.front()->getStart().rc().front().rc().getId());
                }
                if(path.front()->getStart().getId() == start && path.back()->getFinish().getId() == end) {
                    std::cout << "Found path that must must be correct";
                    std::unordered_set<EdgeId> tmp(path.begin(), path.end());
                    for(EdgeId e1 : edges) {
                        if(tmp.find(e1) == tmp.end()) {
                            to_remove.insert(e1);
                            to_remove.insert(e1->rc().getId());
                        }
                    }
                    break;
                }
            }
        }
        mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()});
        std::cout << "cleaned " << mg.size() << " " << mg.edgeCount() << std::endl;
        ag::MergeAll(logger, threads, mg);
//        mg = MultiGraphHelper::MergeAllPaths(mg, true);
        std::cout << "merged " << mg.size() << " " << mg.edgeCount() << std::endl;
    }

    void RemoveSmall(logging::Logger &logger, size_t threads, MultiGraph &mg) {
        DisjointSet<EdgeId> components;
        for(MGVertex &vertex: mg.vertices()) {
            for(MGEdge &e1 : vertex)
                for(MGEdge &e2: vertex.rc()) {
                    components.link(e1.getId(), e2.rc().getId());
                }
        }
        std::unordered_set<EdgeId>all;
        for(MGEdge &edge : mg.edges()) {
            all.insert(edge.getId());
        }
        std::unordered_set<EdgeId> to_remove;
        auto comps  = components.subsets(all);
        for(auto &it : comps) {
            size_t size = 0;
            for(EdgeId e : it.second) {
                size += e->fullSize();
                if(e->getStart().inDeg() != 0)
                    size -= e->getStart().size();
            }
            if(size < 1000000) {
                for(EdgeId e : it.second) {
                    to_remove.insert(e);
                    to_remove.insert(e->rc().getId());
                }
            }
        }
        std::unordered_set<VertexId> vert_to_remove;
        for(MGVertex &vertex : mg.vertices()) {
            if(vertex.inDeg() == 0 && vertex.outDeg() == 0 && vertex.size() < 1000000) {
                vert_to_remove.insert(vertex.getId());
                vert_to_remove.insert(vertex.rc().getId());
            }
        }
        mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()}, {vert_to_remove.begin(), vert_to_remove.end()});
        std::cout << "cleaned " << mg.size() << " " << mg.edgeCount() << std::endl;
        ag::MergeAll(logger, threads, mg);
//        mg = MultiGraphHelper::MergeAllPaths(mg, true);
        std::cout << "merged " << mg.size() << " " << mg.edgeCount() << std::endl;
    }

    void RemoveInversions(logging::Logger &logger, size_t threads, MultiGraph &mg) {
        std::unordered_set<EdgeId> to_remove;
        for(MGVertex &v : mg.vertices()) {
            if(!v.isCanonical())
                continue;
            if(v.inDeg()== 2 && v.outDeg() == 2 && v.front() == v.back().rc() && v.front().fullSize() < 300000) {
                MGEdge &in1 = v.rc().front().rc();
                MGEdge &in2 = v.rc().front().rc();
                MGEdge &b1 = v.front();
                MGEdge &b2 = v.back();
                if(to_remove.find(b1.getId()) != to_remove.end()|| to_remove.find(b2.getId()) != to_remove.end())
                    continue;
                Sequence seq = in1.getSeq() + b1.getSeq().Subseq(v.size()) + in2.rc().getSeq().Subseq(v.size());
                mg.addEdge(in1.getStart(), in2.getFinish().rc(), seq, MGEdgeData(""));
                to_remove.insert(in1.getId());
                to_remove.insert(in2.getId());
                to_remove.insert(b1.getId());
                to_remove.insert(b2.getId());
            }
        }
        mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()});
        std::cout << "cleaned " << mg.size() << " " << mg.edgeCount() << std::endl;
        ag::MergeAll(logger, threads, mg);
//        mg = MultiGraphHelper::MergeAllPaths(mg, true);
        std::cout << "merged " << mg.size() << " " << mg.edgeCount() << std::endl;
    }

protected:
    std::unordered_map<std::string, std::experimental::filesystem::path>
    innerRun(logging::Logger &logger, size_t threads,
             const std::experimental::filesystem::path &dir, bool debug,
             const AlgorithmParameterValues &parameterValues,
             const std::unordered_map<std::string, io::Library> &input) override {
        const std::experimental::filesystem::path &graph_path = input.find("graph")->second.front();
        MultiGraph mg = MultiGraphHelper::LoadGFA(graph_path, true);
        mg = MultiGraphHelper::TransformToEdgeGraph(mg, 501);
        logger << "DBG parameters: " << mg.size() << " vertices, " << mg.edgeCount() << " edges" << std::endl;
        CollapseSimpleBulges(logger, threads, mg);
        ChooseShortcuts(logger, threads, mg);
        RemoveSmall(logger, threads, mg);
        PassAcyclic(logger, threads, mg);
        RemoveInversions(logger, threads, mg);
        MultiGraphHelper::printEdgeGFA(mg, dir / "consensus.gfa");
        MultiGraphHelper::printExtractedContigs(mg, dir / "consensus.fasta", true);
        return {{"graph",    dir / "consensus.gfa"},
                {"assembly", dir / "consensus.fasta"}};
    }
};

int main(int argc, char **argv) {
    ConsensusPhase phase;
    AlgorithmParameters params = phase.getStandaloneParameters();
    CLParser parser(params, {"o=output-dir", "t=threads"});
    LoggedProgram consensus("Consensus", std::move(phase), std::move(parser), "Starting consensus assembly construction", "Finished consensus construction");
    consensus.run(oneline::initialize<std::string, char*>(argv, argv + argc));
    return 0;
}
