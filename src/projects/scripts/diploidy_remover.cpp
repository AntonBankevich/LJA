#include "common/disjoint_sets.hpp"
#include <common/cl_parser.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <common/dir_utils.hpp>
#include "dbg/multi_graph.hpp"
using namespace multigraph;

void
addGoodVertex(Vertex &v, std::unordered_set<multigraph::Vertex *> &good, std::vector<multigraph::Edge *> &candidates) {
    if(good.find(&v) != good.end())
        return;
    std::cout << "Good vertex " << v.id << std::endl;
    good.insert(&v);
    good.insert(v.rc);
    for(Edge *e : v.outgoing) {
        if(good.find(e->end) == good.end())
            candidates.emplace_back(e);
    }
    if(&v != v.rc)
        for(Edge *e : v.rc->outgoing) {
            if(good.find(e->end) == good.end())
                candidates.emplace_back(e);
        }
}

std::vector<multigraph::Edge *> component(Edge *initial) {
    std::unordered_set<Vertex *> visited;
    std::unordered_set<Edge *> res;
    std::vector<Vertex *> queue;
    queue.emplace_back(initial->start);
    queue.emplace_back(initial->end->rc);
    visited.insert(initial->start->rc);
    visited.insert(initial->end);
    while(!queue.empty()) {
        Vertex *v = queue.back();
        queue.pop_back();
        if(visited.find(v) != visited.end())
            continue;
        visited.insert(v);
        for(Edge *e: v->outgoing) {
            if(visited.find(e->end) == visited.end())
                queue.emplace_back(e->end);
            if(visited.find(e->end->rc) == visited.end())
                queue.emplace_back(e->end->rc);
        }
    }
    for(Vertex *v: visited) {
        if(v != initial->start->rc && v != initial->end)
            for(Edge *e: v->outgoing) {
                if(e != initial && e != initial->rc) {
                    if (e->getSeq() <= e->rc->getSeq())
                        res.insert(e);
                    else
                        res.insert(e->rc);
                }
            }
    }
    return {res.begin(), res.end()};
}

void CollapseSimpleBulges(MultiGraph &mg) {
    std::unordered_set<Vertex *> good;
    std::vector<Edge *> candidates;
    std::unordered_set<Edge *> to_remove;
    DisjointSet<Vertex*> paths;
    for(auto &it: mg.vertices) {
        Vertex &v = it.second;
        paths.link(&v, v.rc);
        if(v.outDeg() > 0 && v.outDeg() <= 2) {
            bool ok = true;
            for(Edge *e : v.outgoing) {
                if(e->end != v.outgoing[0]->end) {
                    ok = false;
                    break;
                }
            }
            if(ok && v.outDeg() == v.outgoing[0]->end->inDeg())
                paths.link(&v, v.outgoing[0]->end);
        }
    }
    std::unordered_map<Vertex *, std::vector<Vertex*>> split = paths.nontrivialSubsets();
    for(auto it : split) {
        size_t len = 0;
        for(Vertex *v : it.second) {
            if(v->outDeg() > 0 && paths.get(v->outgoing[0]->end) == it.first) {
                len += 2 * std::min(v->outgoing.front()->size(), v->outgoing.back()->size()) - v->outgoing.front()->end->size() - v->size();
            }
        }
        len /= 2;
        if(len > 200000) {
            for(Vertex *v : it.second) {
                addGoodVertex(*v, good, candidates);
            }
        }
    }
    while(!candidates.empty()) {
        Edge *e = candidates.back();
        candidates.pop_back();
        if(to_remove.find(e) !=to_remove.end())
            continue;
        if(e->start->outDeg() == 1 && e->end->inDeg() == 1) {
            addGoodVertex(*e->end, good, candidates);
            continue;
        }
        if(e->start->outDeg() > 2 || e->end->inDeg() > 2) {
            continue;
        }
        std::cout << "Checking edge " << e->getId() << std::endl;
        std::vector<Edge*> comp = component(e);
        size_t size = 0;
        for(Edge *e1: comp) {
            size+= 2 * e1->size() - e1->start->size() - e1->end->size();
        }
        size_t elen = 2*e->size() - e->start->size() - e->end->size();
        std::cout << elen << " " << comp.size() << " " << size << std::endl;
        if((comp.size() == 1 && e->start->outDeg()==2 && e->end->inDeg()==2 && size <= elen) || size< elen * 1.5+100000) {
            for(Edge *e1 : comp) {
                std::cout << "Deleting edges " << e1->getId() << " and " << e1->rc->getId() << " based on analysis of edge " << e->getId() << std::endl;
                to_remove.insert(e1);
                to_remove.insert(e1->rc);
            }
            addGoodVertex(*e->end, good, candidates);
        }
    }
    mg = mg.DeleteEdges({to_remove.begin(), to_remove.end()});
    std::cout << "cleaned " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.Merge(true);
    std::cout << "merged " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
}

void ChooseShortcuts(MultiGraph &mg) {
    DisjointSet<Vertex *> small_components;
    for(auto &it : mg.edges) {
        if(it.second.size() < 1000000)
            small_components.link(it.second.start, it.second.end);
    }
    std::unordered_set<Vertex *> all_vert;
    for(auto &it : mg.vertices) {
        all_vert.insert(&it.second);
    }
    auto cmap = small_components.subsets(all_vert);
    std::unordered_set<Edge *> to_remove;
    for(auto & it:cmap) {
        std::unordered_set<Vertex *> cvert(it.second.begin(), it.second.end());
        if(cvert.size() > 6)
            continue;
        bool good = true;
        Vertex * start = nullptr;
        Vertex * end = nullptr;
        std::unordered_set<Edge *> edges;
        for(Vertex * v : cvert) {
            for(Edge *e: v->outgoing) {
                if(cvert.find(e->end) == cvert.end()) {
                    if(end != nullptr) {
                        good = false;
                    }
                    end = v;
                } else {
                    edges.insert(e);
                }
            }
            for(Edge *e: v->rc->outgoing) {
                if(cvert.find(e->end->rc) == cvert.end()) {
                    if(start != nullptr) {
                        good = false;
                    }
                    start = v;
                }
            }
        }
        std::cout << "Component " << cvert.size() << " " << edges.size() << std::endl;
        if(edges.size() == 0 || start == nullptr || end == nullptr || !good) {
            continue;
        }
        std::cout << "Processing component" << std::endl;
        for(Edge *e : edges) {
            std::cout << e->getId() << " ";
        }
        std::cout << std::endl;
        Edge * shortcut = nullptr;
        for(Edge *e : edges) {
            if(e->start == start && e->end == end && (shortcut == nullptr || shortcut->size() < e->size())) {
                shortcut = e;
            }
        }
        size_t size = 0;
        for(Edge *e : edges) {
            size += e->size() * 2 - e->start->size() - e->end->size();
        }
        if(size <= 1000000 && (to_remove.find(shortcut) == to_remove.end()) && (shortcut != nullptr || start == end)) {
            if(shortcut != nullptr)
                std::cout << "Found shortcut " << shortcut->getId() << std::endl;
            else
                std::cout << "Found tip " << std::endl;
            for(Edge *e : edges) {
                if(e != shortcut && e->rc != shortcut) {
                    to_remove.insert(e);
                    to_remove.insert(e->rc);
                }
            }
        }
    }
    mg = mg.DeleteEdges({to_remove.begin(), to_remove.end()});
    std::cout << "cleaned " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.Merge(true);
    std::cout << "merged " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
}

void PassAcyclic(MultiGraph &mg) {
    DisjointSet<Vertex *> small_components;
    for(auto &it : mg.edges) {
        if(it.second.size() < 1000000)
            small_components.link(it.second.start, it.second.end);
    }
    std::unordered_set<Vertex *> all_vert;
    for(auto &it : mg.vertices) {
        all_vert.insert(&it.second);
    }
    auto cmap = small_components.subsets(all_vert);
    std::unordered_set<Edge *> to_remove;
    for(auto & it:cmap) {
        std::unordered_set<Vertex *> cvert(it.second.begin(), it.second.end());
        if(cvert.size() > 6)
            continue;
        bool good = true;
        Vertex * start = nullptr;
        Vertex * end = nullptr;
        std::unordered_set<Edge *> edges;
        for(Vertex * v : cvert) {
            for(Edge *e: v->outgoing) {
                if(cvert.find(e->end) == cvert.end()) {
                    if(end != nullptr) {
                        good = false;
                    }
                    end = v;
                } else {
                    edges.insert(e);
                }
            }
            for(Edge *e: v->rc->outgoing) {
                if(cvert.find(e->end->rc) == cvert.end()) {
                    if(start != nullptr) {
                        good = false;
                    }
                    start = v;
                }
            }
        }
        std::cout << "Component " << cvert.size() << " " << edges.size() << std::endl;
        if(edges.size() == 0 || start == nullptr || end == nullptr || start == end || !good) {
            continue;
        }
        std::cout << "Component " << cvert.size() << " " << edges.size() << std::endl;
        for(Edge *e : edges) {
            std::cout << e->getId() << " ";
        }
        std::cout << std::endl;
        size_t size = 0;
        for(Edge *e : edges) {
            size += e->size() * 2 - e->start->size() - e->end->size();
        }
        if(size > 1000000)
            continue;
        for(Edge * e: edges) {
            if(to_remove.find(e) != to_remove.end())
                continue;
            std::vector<Edge *> path = {e};
            while(path.size() < 10 && path.back()->end->outDeg() == 1 && path.back()->end!= end) {
                path.push_back(path.back()->end->outgoing[0]);
            }
            while(path.size() < 10 && path.front()->start->inDeg() == 1 && path.front()->start != start) {
                path.insert(path.begin(), path.front()->start->rc->outgoing[0]->rc);
            }
            if(path.front()->start == start && path.back()->end==end) {
                std::cout << "Found path that must must be correct";
                std::unordered_set<Edge *> tmp(path.begin(), path.end());
                for(Edge * e1 : edges) {
                    if(tmp.find(e1) == tmp.end()) {
                        to_remove.insert(e1);
                        to_remove.insert(e1->rc);
                    }
                }
                break;
            }
        }
    }
    mg = mg.DeleteEdges({to_remove.begin(), to_remove.end()});
    std::cout << "cleaned " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.Merge(true);
    std::cout << "merged " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
}

void RemoveSmall(MultiGraph &mg) {
    DisjointSet<Edge *> components;
    for(auto &it: mg.vertices) {
        for(Edge *e1 : it.second.outgoing)
            for(Edge *e2: it.second.rc->outgoing) {
                components.link(e1, e2->rc);
            }
    }
    std::unordered_set<Edge *>all;
    for(auto & it : mg.edges) {
        all.insert(&it.second);
    }
    std::unordered_set<Edge *> to_remove;
    auto comps  = components.subsets(all);
    for(auto it : comps) {
        size_t size = 0;
        for(Edge * e : it.second) {
            size += e->size();
            if(e->start->inDeg() != 0)
                size -= e->start->size();
        }
        if(size < 1000000) {
            for(Edge * e : it.second) {
                to_remove.insert(e);
                to_remove.insert(e->rc);
            }
        }
    }
    std::unordered_set<Vertex *> vert_to_remove;
    for(auto &it : mg.vertices) {
        if(it.second.inDeg() == 0 && it.second.outDeg() == 0 && it.second.size() < 1000000) {
            vert_to_remove.insert(&it.second);
            vert_to_remove.insert(it.second.rc);
        }
    }
    mg = mg.Delete({vert_to_remove.begin(), vert_to_remove.end()}, {to_remove.begin(), to_remove.end()});
    std::cout << "cleaned " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.Merge(true);
    std::cout << "merged " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
}

void RemoveInversions(MultiGraph &mg) {
    std::unordered_set<Edge *> to_remove;
    for(auto &it : mg.vertices) {
        Vertex *v = &it.second;
        if(!v->isCanonical())
            continue;
        if(v->inDeg()== 2 && v->outDeg() == 2 && v->outgoing[0] == v->outgoing[1]->rc && v->outgoing[0]->size() < 300000) {
            Edge *in1 = v->rc->outgoing[0]->rc;
            Edge *in2 = v->rc->outgoing[1]->rc;
            Edge *b1 = v->outgoing[0];
            Edge *b2 = v->outgoing[1];
            if(to_remove.find(b1) != to_remove.end()|| to_remove.find(b2) != to_remove.end())
                continue;
            Sequence seq = in1->getSeq() + b1->getSeq().Subseq(v->size()) + in2->rc->getSeq().Subseq(v->size());
            mg.addEdge(*in1->start, *in2->end->rc, seq);
            to_remove.insert(in1);
            to_remove.insert(in2);
            to_remove.insert(b1);
            to_remove.insert(b2);
        }
    }
    mg = mg.DeleteEdges({to_remove.begin(), to_remove.end()});
    std::cout << "cleaned " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.Merge(true);
    std::cout << "merged " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
}

int main(int argc, char **argv) {
    MultiGraph mg;
    mg.LoadGFA(argv[1], true);
    std::experimental::filesystem::path dir = argv[2];
    ensure_dir_existance(dir);
    std::cout << "dbg " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    mg = mg.DBG(501);
    std::cout << "initial " << mg.vertices.size() << " " << mg.edges.size() << std::endl;
    CollapseSimpleBulges(mg);
    ChooseShortcuts(mg);
    RemoveSmall(mg);
    PassAcyclic(mg);
    RemoveInversions(mg);


    mg.printEdgeGFA(dir / "consensus.gfa");
    mg.printEdges(dir / "consensus.fasta", true);

    return 0;
}
