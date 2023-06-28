#include "dbg/multi_graph.hpp"
#include <common/dir_utils.hpp>
#include "common/disjoint_sets.hpp"
#include <experimental/filesystem>
#include <unordered_map>
using namespace multigraph;

void addGoodVertex(Vertex &v, std::unordered_set<VertexId> &good, std::vector<EdgeId> &candidates) {
    if(good.find(v.getId()) != good.end())
        return;
    std::cout << "Good vertex " << v.getId() << std::endl;
    good.insert(v.getId());
    good.insert(v.rc().getId());
    for(Edge &e : v) {
        if(good.find(e.end().getId()) == good.end())
            candidates.emplace_back(e.getId());
    }
    if(v != v.rc())
        for(Edge &e : v.rc()) {
            if(good.find(e.end().getId()) == good.end())
                candidates.emplace_back(e.getId());
        }
}

std::vector<EdgeId> component(Edge &initial) {
    std::unordered_set<VertexId> visited;
    std::unordered_set<EdgeId> res;
    std::vector<VertexId> queue;
    queue.emplace_back(initial.start().getId());
    queue.emplace_back(initial.end().rc().getId());
    visited.insert(initial.start().rc().getId());
    visited.insert(initial.end().getId());
    while(!queue.empty()) {
        Vertex &v = *queue.back();
        queue.pop_back();
        if(visited.find(v.getId()) != visited.end())
            continue;
        visited.insert(v.getId());
        for(Edge &e: v) {
            if(visited.find(e.end().getId()) == visited.end())
                queue.emplace_back(e.end().getId());
            if(visited.find(e.end().rc().getId()) == visited.end())
                queue.emplace_back(e.end().rc().getId());
        }
    }
    for(VertexId v: visited) {
        if(v != initial.start().rc().getId() && v != initial.end().getId())
            for(Edge &e: *v) {
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

void CollapseSimpleBulges(MultiGraph &mg) {
    std::unordered_set<VertexId> good;
    std::vector<EdgeId> candidates;
    std::unordered_set<ConstEdgeId> to_remove;
    DisjointSet<VertexId> paths;
    for(Vertex &v: mg.vertices()) {
        paths.link(v.getId(), v.rc().getId());
        if(v.outDeg() > 0 && v.outDeg() <= 2) {
            bool ok = true;
            for(Edge &e : v) {
                if(e.end() != v[0].end()) {
                    ok = false;
                    break;
                }
            }
            if(ok && v.outDeg() == v[0].end().inDeg())
                paths.link(v.getId(), v[0].end().getId());
        }
    }
    std::unordered_map<VertexId, std::vector<VertexId>> split = paths.nontrivialSubsets();
    for(auto &it : split) {
        size_t len = 0;
        for(VertexId v : it.second) {
            if(v->outDeg() > 0 && paths.get(v->begin()->end().getId()) == it.first) {
                len += 2 * std::min(v->front().size(), v->back().size()) - v->front().end().size() - v->size();
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
        if(e->start().outDeg() == 1 && e->end().inDeg() == 1) {
            addGoodVertex(e->end(), good, candidates);
            continue;
        }
        if(e->start().outDeg() > 2 || e->end().inDeg() > 2) {
            continue;
        }
        std::cout << "Checking edge " << e->getId() << std::endl;
        std::vector<EdgeId> comp = component(*e);
        size_t size = 0;
        for(EdgeId e1: comp) {
            size+= 2 * e1->size() - e1->start().size() - e1->end().size();
        }
        size_t elen = 2*e->size() - e->start().size() - e->end().size();
        std::cout << elen << " " << comp.size() << " " << size << std::endl;
        if((comp.size() == 1 && e->start().outDeg()==2 && e->end().inDeg()==2 && size <= elen) || size< elen * 1.5+100000) {
            for(EdgeId e1 : comp) {
                std::cout << "Deleting edges " << e1->getId() << " and " << e1->rc().getId() << " based on analysis of edge " << e->getId() << std::endl;
                to_remove.emplace(e1->getId());
                to_remove.emplace(e1->rc().getId());
            }
            addGoodVertex(e->end(), good, candidates);
        }
    }
    mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()});
    std::cout << "cleaned " << mg.size() << " " << mg.edgeNumber() << std::endl;
    mg = MultiGraphHelper::MergeAllPaths(mg, true);
    std::cout << "merged " << mg.size() << " " << mg.edgeNumber() << std::endl;
}

void ChooseShortcuts(MultiGraph &mg) {
    DisjointSet<VertexId> small_components;
    for(Edge &edge : mg.edges()) {
        if(edge.size() < 1000000)
            small_components.link(edge.start().getId(), edge.end().getId());
    }
    std::unordered_set<VertexId> all_vert;
    for(Vertex &vertex : mg.vertices()) {
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
            for(Edge &e: *v) {
                if(cvert.find(e.end().getId()) == cvert.end()) {
                    if(end.valid()) {
                        good = false;
                    }
                    end = v;
                } else {
                    edges.insert(e.getId());
                }
            }
            for(Edge &e: v->rc()) {
                if(cvert.find(e.end().rc().getId()) == cvert.end()) {
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
            if(e->start().getId() == start && e->end().getId() == end && (!shortcut.valid() || shortcut->size() < e->size())) {
                shortcut = e;
            }
        }
        size_t size = 0;
        for(EdgeId e : edges) {
            size += e->size() * 2 - e->start().size() - e->end().size();
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
    std::cout << "cleaned " << mg.size() << " " << mg.edgeNumber() << std::endl;
    mg = MultiGraphHelper::MergeAllPaths(mg, true);
    std::cout << "merged " << mg.size() << " " << mg.edgeNumber() << std::endl;
}

void PassAcyclic(MultiGraph &mg) {
    DisjointSet<VertexId> small_components;
    for(Edge &edge : mg.edges()) {
        if(edge.size() < 1000000)
            small_components.link(edge.start().getId(), edge.end().getId());
    }
    std::unordered_set<VertexId> all_vert;
    for(Vertex &vertex : mg.vertices()) {
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
            for(Edge &e: *v) {
                if(cvert.find(e.end().getId()) == cvert.end()) {
                    if(end.valid()) {
                        good = false;
                    }
                    end = v;
                } else {
                    edges.emplace(e.getId());
                }
            }
            for(Edge &e: v->rc()) {
                if(cvert.find(e.end().rc().getId()) == cvert.end()) {
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
            size += e->size() * 2 - e->start().size() - e->end().size();
        }
        if(size > 1000000)
            continue;
        for(EdgeId e: edges) {
            if(to_remove.find(ConstEdgeId(e)) != to_remove.end())
                continue;
            std::vector<EdgeId> path = {e};
            while(path.size() < 10 && path.back()->end().outDeg() == 1 && path.back()->end().getId() != end) {
                path.push_back(path.back()->end().front().getId());
            }
            while(path.size() < 10 && path.front()->start().inDeg() == 1 && path.front()->start().getId() != start) {
                path.insert(path.begin(), path.front()->start().rc().front().rc().getId());
            }
            if(path.front()->start().getId() == start && path.back()->end().getId() == end) {
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
    std::cout << "cleaned " << mg.size() << " " << mg.edgeNumber() << std::endl;
    mg = MultiGraphHelper::MergeAllPaths(mg, true);
    std::cout << "merged " << mg.size() << " " << mg.edgeNumber() << std::endl;
}

void RemoveSmall(MultiGraph &mg) {
    DisjointSet<EdgeId> components;
    for(Vertex &vertex: mg.vertices()) {
        for(Edge &e1 : vertex)
            for(Edge &e2: vertex.rc()) {
                components.link(e1.getId(), e2.rc().getId());
            }
    }
    std::unordered_set<EdgeId>all;
    for(Edge &edge : mg.edges()) {
        all.insert(edge.getId());
    }
    std::unordered_set<EdgeId> to_remove;
    auto comps  = components.subsets(all);
    for(auto &it : comps) {
        size_t size = 0;
        for(EdgeId e : it.second) {
            size += e->size();
            if(e->start().inDeg() != 0)
                size -= e->start().size();
        }
        if(size < 1000000) {
            for(EdgeId e : it.second) {
                to_remove.insert(e);
                to_remove.insert(e->rc().getId());
            }
        }
    }
    std::unordered_set<VertexId> vert_to_remove;
    for(Vertex &vertex : mg.vertices()) {
        if(vertex.inDeg() == 0 && vertex.outDeg() == 0 && vertex.size() < 1000000) {
            vert_to_remove.insert(vertex.getId());
            vert_to_remove.insert(vertex.rc().getId());
        }
    }
    mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()}, {vert_to_remove.begin(), vert_to_remove.end()});
    std::cout << "cleaned " << mg.size() << " " << mg.edgeNumber() << std::endl;
    mg = MultiGraphHelper::MergeAllPaths(mg, true);
    std::cout << "merged " << mg.size() << " " << mg.edgeNumber() << std::endl;
}

void RemoveInversions(MultiGraph &mg) {
    std::unordered_set<EdgeId> to_remove;
    for(Vertex &v : mg.vertices()) {
        if(!v.isCanonical())
            continue;
        if(v.inDeg()== 2 && v.outDeg() == 2 && v[0] == v[1].rc() && v[0].size() < 300000) {
            Edge &in1 = v.rc()[0].rc();
            Edge &in2 = v.rc()[1].rc();
            Edge &b1 = v[0];
            Edge &b2 = v[1];
            if(to_remove.find(b1.getId()) != to_remove.end()|| to_remove.find(b2.getId()) != to_remove.end())
                continue;
            Sequence seq = in1.getSeq() + b1.getSeq().Subseq(v.size()) + in2.rc().getSeq().Subseq(v.size());
            mg.addEdge(in1.start(), in2.end().rc(), seq);
            to_remove.insert(in1.getId());
            to_remove.insert(in2.getId());
            to_remove.insert(b1.getId());
            to_remove.insert(b2.getId());
        }
    }
    mg = MultiGraphHelper::Delete(mg, {to_remove.begin(), to_remove.end()});
    std::cout << "cleaned " << mg.size() << " " << mg.edgeNumber() << std::endl;
    mg = MultiGraphHelper::MergeAllPaths(mg, true);
    std::cout << "merged " << mg.size() << " " << mg.edgeNumber() << std::endl;
}

int main(int argc, char **argv) {
    MultiGraph mg = MultiGraphHelper::LoadGFA(argv[1], true);
    std::experimental::filesystem::path dir = argv[2];
    ensure_dir_existance(dir);
    std::cout << "dbg " << mg.size() << " " << mg.edgeNumber() << std::endl;
    mg = MultiGraphHelper::TransformToVertexGraph(mg, 501);
    std::cout << "initial " << mg.size() << " " << mg.edgeNumber() << std::endl;
    CollapseSimpleBulges(mg);
    ChooseShortcuts(mg);
    RemoveSmall(mg);
    PassAcyclic(mg);
    RemoveInversions(mg);


    MultiGraphHelper::printEdgeGFA(mg, dir / "consensus.gfa");
    MultiGraphHelper::printExtractedContigs(mg, dir / "consensus.fasta", true);

    return 0;
}
