#pragma once

#include <sequences/sequence.hpp>
#include <unordered_set>
#include <unordered_map>
#include <experimental/filesystem>
#include <fstream>
#include <common/string_utils.hpp>
#include <sequences/contigs.hpp>

namespace multigraph {
    class Edge;
    struct Vertex {
        Sequence seq;
        int id;
        std::vector<Edge *> outgoing;
        Vertex *rc = nullptr;
        explicit Vertex(const Sequence &seq, int id) : seq(seq), id(id) {
            outgoing.reserve(4);
        }

        size_t inDeg() const {
            return rc->outgoing.size();
        }

        size_t outDeg() const {
            return outgoing.size();
        }
    };

    struct Edge {
    private:
        Sequence seq;
        int id;
        size_t sz;
        bool canonical;
    public:
        Vertex *start = nullptr;
        Vertex *end = nullptr;
        Edge *rc = nullptr;
        explicit Edge(const Sequence &seq, int id = 0) : seq(seq), id(id), sz(seq.size()), canonical(seq <= !seq) {
        }
        Sequence getSeq() const {
            if(seq.empty())
                return start->seq + end->seq.Subseq(sz);
            return seq;
        }

        int getId() const {
            return id;
        }

        size_t size() const {
            return sz;
        }

        size_t overlap() const {
            VERIFY(start->seq.size() + end->seq.size() > sz);
            return start->seq.size() + end->seq.size() - sz;
        }

        bool isCanonical() const {
            return canonical;
        }
    };

    struct MultiGraph {
        int maxVId = 0;
        int maxEId = 0;
        std::vector<Vertex *> vertices;
        std::vector<Edge *> edges;
        MultiGraph() = default;
        MultiGraph(MultiGraph &&other) = default;
        MultiGraph &operator=(MultiGraph &&other) = default;
        MultiGraph(const MultiGraph &) = delete;

        void LoadGFA(const std::experimental::filesystem::path &gfa_file) {
            std::ifstream is;
            is.open(gfa_file);
            std::unordered_map<std::string, Vertex*> vmap;
            for( std::string line; getline(is, line); ) {
                std::vector<std::string> tokens = ::split(line);
                if(tokens[0] == "S") {
                    Vertex &newV = addVertex(Sequence(tokens[2]));
                    std::string name = tokens[1];
                    VERIFY(vmap.find(name) == vmap.end());
                    vmap[name] = &newV;
                } else if(tokens[0] == "L") {
                    Vertex *v1 = vmap[tokens[1]];
                    Vertex *v2 = vmap[tokens[3]];
                    if(tokens[2] == "-")
                        v1 = v1->rc;
                    if(tokens[4] == "-")
                        v2 = v2->rc;
                    size_t overlap = std::stoull(tokens[5].substr(0, tokens[5].size() - 1));
                    if(v1->seq.Subseq(v1->seq.size() - overlap) != v2->seq.Subseq(0, overlap)) {
                        v1 = v1->rc;
                    }
                    VERIFY(v1->seq.Subseq(v1->seq.size() - overlap) == v2->seq.Subseq(0, overlap));
                    addEdge(*v1, *v2, v1->seq + v2->seq.Subseq(overlap));
                }
            }
            is.close();
        }

        MultiGraph DBG() const {
            MultiGraph dbg;
            std::unordered_map<Edge *, Vertex *> emap;
            for(Vertex * v : vertices) {
                if(v->outDeg() == 0 || emap.find(v->outgoing[0]) != emap.end()) {
                    continue;
                }
                Vertex *newv = &dbg.addVertex(v->seq.Subseq(v->seq.size() - v->outgoing[0]->overlap()));
                for(Edge * edge : v->outgoing) {
                    Vertex * right = edge->end;
                    for(Edge * edge1 : right->rc->outgoing) {
                        emap[edge1] = newv->rc;
                        emap[edge1->rc] = newv;
                    }
                }
            }
            for(Vertex * v : vertices) {
                if(!(v->seq <= !v->seq))
                    continue;
                Vertex * start = nullptr;
                Vertex * end = nullptr;
                if(v->inDeg() == 0) {
                    start = &dbg.addVertex(v->seq.Subseq(0, 4001));
                } else {
                    start = emap[v->rc->outgoing[0]]->rc;
                }
                if(v->outDeg() == 0) {
                    end = &dbg.addVertex(v->seq.Subseq(v->seq.size() - 4001));
                } else {
                    end = emap[v->outgoing[0]];
                }
                dbg.addEdge(*start, *end, v->seq, v->id);
            }
            return std::move(dbg);
        }

        ~MultiGraph() {
            for(Vertex * &v : vertices) {
                delete v;
                v = nullptr;
            }
            for(Edge * &e : edges) {
                delete e;
                e = nullptr;
            }
        }

        MultiGraph DeleteEdges(const std::unordered_set<Edge *> &to_delete) const {
            MultiGraph mg;
            std::unordered_map<Vertex *, Vertex *> vmap;
            std::unordered_set<Edge *> visited;
            for(Vertex * v : vertices) {
                if(vmap.find(v) != vmap.end())
                    continue;
                vmap[v] = &mg.addVertex(v->seq, v->id);
                vmap[v->rc] = vmap[v]->rc;
            }
            for(Edge * edge : edges) {
                if(!edge->isCanonical() || to_delete.find(edge) != to_delete.end())
                    continue;
                mg.addEdge(*vmap[edge->start], *vmap[edge->end], edge->getSeq(), edge->getId());
            }
            return std::move(mg);
        }

        MultiGraph BulgeSubgraph() const {
            std::unordered_set<Vertex *> good;
            std::unordered_set<Edge *> to_delete;
            size_t sz = 0;
            for(Vertex *v : vertices) {
                if(v->outDeg() != 2) {
                    continue;
                }
                if(v->outgoing[0]->end == v->outgoing[1]->end) {
                    good.emplace(v);
                    good.emplace(v->rc);
                }
            }
            size_t bulges = 0;
            for(Vertex *v : vertices) {
                if(v->outDeg() != 2) {
                    continue;
                }
                if(to_delete.find(v->outgoing[0]) != to_delete.end() ||
                        to_delete.find(v->outgoing[1]) != to_delete.end()) {
                    continue;
                }
                Edge * todel = nullptr;
                if(v->outgoing[0]->end == v->outgoing[1]->end) {
                    todel = v->outgoing[0];
                    bulges++;
                    sz += v->outgoing[0]->size();
                } else {
                    if((good.find(v) == good.end() || v->outgoing[0]->size() < 1000000) && v->outgoing[0]->end->outDeg() == 0) {
                        todel = v->outgoing[0];
                    } else if((good.find(v) == good.end() || v->outgoing[0]->size() < 1000000) && v->outgoing[1]->end->outDeg() == 0) {
                        todel = v->outgoing[1];
                    } else if(good.find(v->outgoing[0]->end)== good.end()) {
                        todel = v->outgoing[0];
                    } else if(good.find(v->outgoing[1]->end)== good.end()) {
                        todel = v->outgoing[1];
                    }
                }
                if(todel != nullptr) {
                    to_delete.emplace(todel);
                    to_delete.emplace(todel->rc);
                    sz += todel->size();
                }
            }
            std::cout << "Deleting " << sz << " " << to_delete.size() / 2 << std::endl;
            std::cout << "Bulges " << bulges << std::endl;
            return DeleteEdges(to_delete);
        }



        void checkConsistency() const {
            std::unordered_set<Edge const *> eset;
            std::unordered_set<Vertex const *> vset;
            for(Edge const * edge: edges) {
                eset.emplace(edge);
                VERIFY(edge->rc->start == edge->end->rc);
                VERIFY(edge->rc->rc == edge);
            }
            for(Edge const *edge : edges) {
                VERIFY(eset.find(edge->rc) != eset.end());
            }
            for(Vertex const * v: vertices) {
                vset.emplace(v);
                VERIFY(v->rc->rc == v);
                for(Edge *edge : v->outgoing) {
                    VERIFY(eset.find(edge) != eset.end());
                    VERIFY(edge->start == v);
                }
            }
            for(Vertex const *v : vertices) {
                VERIFY(vset.find(v->rc) != vset.end());
            }
        }

        Vertex &addVertex(const Sequence &seq, int id = 0) {
            if(id = 0) {
                id = maxVId + 1;
            }
            maxVId = std::max(std::abs(id), maxVId);
            vertices.emplace_back(new Vertex(seq, id));
            Vertex *res = vertices.back();
            Vertex *rc = res;
            if(seq != !seq) {
                vertices.emplace_back(new Vertex(!seq, -id));
                rc = vertices.back();
            }
            res->rc = rc;
            rc->rc = res;
            return *res;
        }

        Edge &addEdge(Vertex &from, Vertex &to, const Sequence &seq, int id = 0) {
            if(id == 0) {
                id = maxEId + 1;
            }
            maxEId = std::max(std::abs(id), maxEId);
            edges.emplace_back(new Edge(seq, id));
            Edge *res = edges.back();
            res->start = &from;
            res->end = &to;
            res->start->outgoing.emplace_back(res);
            if(seq != !seq) {
                edges.emplace_back(new Edge(!seq, -id));
                Edge *rc = edges.back();
                rc->start = to.rc;
                rc->end = from.rc;
                res->rc = rc;
                rc->rc = res;
                res->rc->start->outgoing.emplace_back(res->rc);
            } else {
                res->rc = res;
            }
            return *res;
        }

        std::vector<Edge *> uniquePathForward(Edge &edge) {
            std::vector<Edge *> res = {&edge};
            Vertex *cur = edge.end;
            while(cur != edge.start && cur->inDeg() == 1 && cur->outDeg() == 1) {
                res.emplace_back(cur->outgoing[0]);
                cur = res.back()->end;
            }
            return std::move(res);
        }

        std::vector<Edge *> uniquePath(Edge &edge) {
            std::vector<Edge *> path = uniquePathForward(*edge.rc);
            return uniquePathForward(*path.back()->rc);
        }

        MultiGraph Merge() {
            MultiGraph res;
            std::unordered_set<Edge *> used;
            std::unordered_map<Vertex *, Vertex *> old_to_new;
            for(Edge *edge : edges) {
                if(used.find(edge) != used.end())
                    continue;
                std::vector<Edge *> unique = uniquePath(*edge);
                for(Edge *e : unique) {
                    used.emplace(e);
                    used.emplace(e->rc);
                }
                Vertex *old_start = unique.front()->start;
                Vertex *old_end = unique.back()->end;
                Vertex *new_start = nullptr;
                Vertex *new_end = nullptr;
                if(old_to_new.find(old_start) == old_to_new.end()) {
                    old_to_new[old_start] = &res.addVertex(old_start->seq);
                    old_to_new[old_start->rc] = old_to_new[old_start]->rc;
                }
                new_start = old_to_new[old_start];
                if(old_to_new.find(old_end) == old_to_new.end()) {
                    old_to_new[old_end] = &res.addVertex(old_end->seq);
                    old_to_new[old_end->rc] = old_to_new[old_end]->rc;
                }
                new_end = old_to_new[old_end];
                Sequence new_seq;
                if(unique.size() == 1) {
                    new_seq = unique[0]->getSeq();
                } else {
                    SequenceBuilder sb;
                    sb.append(old_start->seq);
                    for(Edge *e : unique) {
                        sb.append(e->getSeq().Subseq(e->start->seq.size()));
                    }
                    new_seq = sb.BuildSequence();
                }
                res.addEdge(*new_start, *new_end, new_seq);
            }
            for(Vertex *vertex : vertices) {
                if(vertex->inDeg() == 0 && vertex->outDeg() == 0 && vertex->seq <= !vertex->seq) {
                    res.addVertex(vertex->seq);
                }
            }
            res.checkConsistency();
            return std::move(res);
        }

        std::vector<Contig> getCutEdges() {
            std::unordered_map<Vertex *, size_t> cut;
            for(Vertex *v : vertices) {
                if(v->seq <= !v->seq) {
                    if(v->outDeg() == 1) {
                        cut[v] = 0;
                    } else {
                        cut[v] = 1;
                    }
                    cut[v->rc] = 1 - cut[v];
                }
            }
            std::vector<Contig> res;
            size_t cnt = 1;
            for(Edge *edge : edges) {
                if(edge->isCanonical()) {
                    size_t cut_left = edge->start->seq.size() * cut[edge->start];
                    size_t cut_right = edge->end->seq.size() * (1 - cut[edge->end]);
                    if(cut_left + cut_right + 1000 >= edge->size()) {
                        continue;
                    }
                    res.emplace_back(edge->getSeq().Subseq(cut_left, edge->size() - cut_right), itos(cnt));
                    cnt++;
                }
            }
            return std::move(res);
        }

        void printCutEdges(const std::experimental::filesystem::path &f) {
            std::ofstream os;
            os.open(f);
            for(const Contig &contig : getCutEdges()) {
                os << ">" << contig.id << "\n" << contig.seq << "\n";
            }
            os.close();
        }

        void printDot(const std::experimental::filesystem::path &f) {
            std::ofstream os;
            os.open(f);
            os << "digraph {\nnodesep = 0.5;\n";
            std::unordered_map<Vertex *, int> vmap;
            for(Vertex *vertex : vertices) {
                if(vertex->seq <= !vertex->seq) {
                    vmap[vertex] = vmap.size() + 1;
                    if(vertex->rc != vertex)
                        vmap[vertex->rc] = -vmap[vertex];
                }
            }
            for(Vertex *vertex : vertices) {
                os << vmap[vertex] << " [label=\"" << vertex->seq.size() << "\" style=filled fillcolor=\"white\"]\n";
            }
            std::unordered_map<Edge *, std::string> eids;
            for (Edge *edge : edges) {
                os << "\"" << vmap[edge->start] << "\" -> \"" << vmap[edge->end] <<
                   "\" [label=\"" << edge->size() << "\" color = \"black\"]\n" ;
            }
            os << "}\n";
            os.close();
        }

        void printEdgeGFA(const std::experimental::filesystem::path &f, const std::vector<Vertex *> &component) const {
            std::ofstream os;
            os.open(f);
            os << "H\tVN:Z:1.0" << std::endl;
            std::unordered_map<Edge *, std::string> eids;
            for(Vertex *v : component)
                for (Edge *edge : v->outgoing) {
                    if (edge->isCanonical()) {
                        eids[edge] = itos(edge->getId());
                        eids[edge->rc] = itos(edge->getId());
                        os << "S\t" << eids[edge] << "\t" << edge->getSeq() << "\n";
                    }
                }
            for (Vertex *vertex : component) {
                if(!(vertex->seq <= !vertex->seq))
                    continue;
                for (Edge *out_edge : vertex->outgoing) {
                    std::string outid = eids[out_edge];
                    bool outsign = out_edge->isCanonical();
                    for (Edge *inc_edge : vertex->rc->outgoing) {
                        std::string incid = eids[inc_edge];
                        bool incsign = inc_edge->rc->isCanonical();
                        os << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                            << (outsign ? "+" : "-") << "\t" << vertex->seq.size() << "M" << "\n";
                    }
                }
            }
            os.close();
        }

        void printEdgeGFA(const std::experimental::filesystem::path &f) const {
            printEdgeGFA(f, vertices);
        }

        void printVertexGFA(const std::experimental::filesystem::path &f, const std::vector<Vertex *> &component) const {
            std::ofstream os;
            os.open(f);
            os << "H\tVN:Z:1.0" << std::endl;
            size_t cnt = 1;
            std::unordered_map<Vertex *, std::string> vids;
            for(Vertex *v : component)
                if(v->seq <= !v->seq) {
                    vids[v] = itos(v->id);
                    vids[v->rc] = itos(v->id);
                    os << "S\t" << vids[v] << "\t" << v->seq << "\n";
                    cnt++;
                }
            for(Vertex *v : component)
                for (Edge *edge : v->outgoing) {
                    if (edge->isCanonical()) {
                        bool incsign = v->seq <= !v->seq;
                        bool outsign = edge->end->seq <= !edge->end->seq;
                        os  << "L\t" << vids[v] << "\t" << (incsign ? "+" : "-") << "\t"
                            << vids[edge->end] << "\t" << (outsign ? "+" : "-") << "\t"
                            << (v->seq.size() + edge->end->seq.size() - edge->size()) << "M" << "\n";
                    }
                }
            os.close();
        }

        void printVertexGFA(const std::experimental::filesystem::path &f) const {
            printVertexGFA(f, vertices);
        }

        std::vector<std::vector<Vertex *>> split() const {
            std::vector<std::vector<Vertex *>> res;
            std::unordered_set<Vertex *> visited;
            for(Vertex *v : vertices) {
                if(visited.find(v) != visited.end())
                    continue;
                std::vector<Vertex *> stack = {v};
                res.emplace_back(std::vector<Vertex *>());
                while(!stack.empty()) {
                    Vertex *p = stack.back();
                    stack.pop_back();
                    if(visited.find(p) != visited.end())
                        continue;
                    visited.emplace(p);
                    visited.emplace(p->rc);
                    res.back().emplace_back(p);
                    res.back().emplace_back(p->rc);
                    for(Edge *e : p->outgoing)
                        stack.emplace_back(e->end);
                    for(Edge *e : p->rc->outgoing)
                        stack.emplace_back(e->end);
                }
            }
            return std::move(res);
        }
    };
}