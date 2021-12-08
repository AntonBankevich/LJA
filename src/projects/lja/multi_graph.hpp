#pragma once

#include <sequences/sequence.hpp>
#include <unordered_set>
#include <unordered_map>
#include <set>
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
        std::string label;
        Vertex *rc = nullptr;
        explicit Vertex(const Sequence &seq, int id, std::string label = "") : seq(seq), id(id), label(label) {
            outgoing.reserve(4);
        }

        bool isCanonical() const {
            return seq <= !seq;
        }

        size_t inDeg() const {
            return rc->outgoing.size();
        }

        size_t outDeg() const {
            return outgoing.size();
        }

        size_t k() const {
            return seq.size();
        };
    };

    struct Edge {
    private:
        Sequence seq;
        int id;
        size_t sz;
        bool canonical;
        std::string label;
    public:
        Vertex *start = nullptr;
        Vertex *end = nullptr;
        Edge *rc = nullptr;
        explicit Edge(const Sequence &seq, int id = 0, std::string label = "") : seq(seq), id(id), sz(seq.size()), canonical(seq <= !seq), label(label) {
        }
        Sequence getSeq() const {
            if(seq.empty())
                return start->seq + end->seq.Subseq(sz);
            return seq;
        }

        int getId() const {
            return id;
        }

        std::string getLabel() const {
            return label;
        }

        size_t size() const {
            return sz;
        }

        size_t overlap() const {
            VERIFY(start->seq.size() + end->seq.size() > sz);
            return start->seq.size() + end->seq.size() - sz;
        }

        int internalLength() const {
            return start->seq.size() + end->seq.size() - sz;
        }

        bool isCanonical() const {
            VERIFY(canonical == (id > 0));
            return canonical;
        }

        bool isTip() const {
            if (start->inDeg() == 0 && end->inDeg() > 1 && end->outDeg() > 0)
                return true;
            if (end->outDeg() == 0 && start->outDeg() > 1 && start->inDeg() > 0)
                return true;
            return false;
        }

//simplified check, works only for trivial cases
        bool isBridge() const {
            if (isTip())
                return false;
            for (auto alt_e: start->outgoing) {
                if (alt_e->getId() != getId() and !alt_e->isTip())
                    return false;
            }
            for (auto alt_e: rc->start->outgoing) {
                if (alt_e->getId() != rc->getId() and !alt_e->isTip())
                    return false;
            }
            return true;
        }

    };

    struct MultiGraph {
        int maxVId = 0;
        int maxEId = 0;
        std::unordered_map<int, Vertex *> vertices;
        std::unordered_map<int, Edge *> edges;
//To canonic ID
        std::unordered_map<std::string, int> label_to_id;
//        std::unordered_map<int, Edge* > id_to_edge;


        MultiGraph() = default;
        MultiGraph(MultiGraph &&other) = default;
        MultiGraph &operator=(MultiGraph &&other) = default;
        MultiGraph(const MultiGraph &) = delete;

        void FillMaps() {
            for (auto p: edges) {
                Edge * e = p.second;
                if (e->isCanonical()) {
                    label_to_id[e->getLabel()] = e->getId();
                }
//                id_to_edge[e->getId()] = e;
            }
        }

        MultiGraph &LoadGFA(const std::experimental::filesystem::path &gfa_file, bool int_ids) {
            std::ifstream is;
            is.open(gfa_file);
            std::unordered_map<std::string, Vertex*> vmap;
            for( std::string line; getline(is, line); ) {
                std::vector<std::string> tokens = ::split(line);
                if(tokens[0] == "S") {
                    std::string name = tokens[1];
                    Vertex &newV = int_ids ? addVertex(Sequence(tokens[2]), std::stoi(name), name): addVertex(Sequence(tokens[2]))

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
            return *this;
        }

//Is it obsolate?
        MultiGraph DBG() const {
            MultiGraph dbg;
            std::unordered_map<Edge *, Vertex *> emap;
            for(auto p : vertices) {
                Vertex* v = p.second;
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
            for(auto p : vertices) {
                Vertex* v = p.second;
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
                dbg.addEdge(*start, *end, v->seq, v->id, v->label);
            }
            std::cout <<"loaded V/E " << vertices.size() << " " << edges.size() << endl;
            std::cout <<"transformed V/E " << dbg.vertices.size() << " " << dbg.edges.size() << endl;
            return std::move(dbg);
        }

        ~MultiGraph() {
            for(auto p : vertices) {
                delete p.second;
                p.second = nullptr;
            }
            for(auto p : edges) {
                delete p.second;
                p.second = nullptr;
            }
        }

        MultiGraph DeleteEdges(const std::unordered_set<Edge *> &to_delete) const {
            MultiGraph mg;
            std::unordered_map<Vertex *, Vertex *> vmap;
            std::unordered_set<Edge *> visited;
            for(auto p : vertices) {
                Vertex* v = p.second;

                if(vmap.find(v) != vmap.end())
                    continue;
                vmap[v] = &mg.addVertex(v->seq, v->id);
                vmap[v->rc] = vmap[v]->rc;
            }
            for(auto p : edges) {
                Edge * edge = p.second;
                if(!edge->isCanonical() || to_delete.find(edge) != to_delete.end())
                    continue;
                mg.addEdge(*vmap[edge->start], *vmap[edge->end], edge->getSeq(), edge->getId(), edge->getLabel());
            }
            return std::move(mg);
        }

        MultiGraph BulgeSubgraph() const {
            std::unordered_set<Vertex *> good;
            std::unordered_set<Edge *> to_delete;
            size_t sz = 0;
            for(auto p : vertices) {
                Vertex* v = p.second;
                if(v->outDeg() != 2) {
                    continue;
                }
                if(v->outgoing[0]->end == v->outgoing[1]->end) {
                    good.emplace(v);
                    good.emplace(v->rc);
                }
            }
            size_t bulges = 0;
            for(auto p : vertices) {
                Vertex* v = p.second;

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
            for(const auto p: edges) {
                Edge * edge = p.second;
                eset.emplace(edge);
                VERIFY(edge->rc->start == edge->end->rc);
                VERIFY(edge->rc->rc == edge);
            }
            for(const auto p: edges) {
                Edge * edge = p.second;
                VERIFY(eset.find(edge->rc) != eset.end());
            }
            for(const auto p : vertices) {
                Vertex * v = p.second;
                vset.emplace(v);
                VERIFY(v->rc->rc == v);
                for(Edge *edge : v->outgoing) {
                    VERIFY(eset.find(edge) != eset.end());
                    VERIFY(edge->start == v);
                }
            }
            for(const auto p : vertices) {
                Vertex * v = p.second;
                VERIFY(vset.find(v->rc) != vset.end());
            }
        }


        Vertex &addVertex(const Sequence &seq, int id = 0, std::string label = "") {
            if (id == 0) {
                id = maxVId + 1;
            }
            maxVId = std::max(std::abs(id), maxVId);
            vertices[id] = new Vertex(seq, id, label);
            Vertex *res = vertices[id];
            Vertex *rc = res;
            if(seq != !seq) {
                vertices[-id] = new Vertex(!seq, -id, label);
                rc = vertices[-id];
            }
            res->rc = rc;
            rc->rc = res;
            return *res;
        }

        Edge &addEdge(Vertex &from, Vertex &to, const Sequence &seq, int id = 0, std::string label = "") {
            if(id == 0) {
                id = maxEId + 1;
            }
            maxEId = std::max(std::abs(id), maxEId);
            if(!(seq <= !seq)) {
                id = -id;
            }

            edges[id] = new Edge(seq, id, label);
            Edge *res = edges[id];
            res->start = &from;
            res->end = &to;
            res->start->outgoing.emplace_back(res);
            if(seq != !seq) {
                edges[-id] = new Edge(!seq, -id, label);
                Edge *rc = edges[-id];
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

        void internalRemoveEdge (int eid) {
            Edge * edge = edges[eid];
            Edge * rc_edge = edge->rc;
            int rcid = rc_edge->getId();
            std::unordered_map<Vertex* , Edge*> to_delete;
            to_delete[edge->start] = edge;
            to_delete[rc_edge->start] = rc_edge;
            for (auto p: to_delete){
                Vertex * d_vertex = p.first;
                Edge* d_edge = p.second;
                d_vertex->outgoing.erase(std::remove(d_vertex->outgoing.begin(), d_vertex->outgoing.end(), d_edge),
                        d_vertex->outgoing.end());
            }
            delete edge;
            delete rc_edge;
            edges.erase(eid);
            edges.erase(rcid);
        }

        void compressVertex(int vid) {
            if (vertices.find(vid) != vertices.end() && vertices[vid]->outDeg() == 1 && vertices[vid]->inDeg() == 1) {
//Do not compress 1-1 loops and RC loops
                if (vertices[vid]->outgoing[0] == vertices[vid]->rc->outgoing[0]->rc)
                    return;
                if (vertices[vid]->outgoing[0] == vertices[vid]->rc->outgoing[0])
                    return;
                int rcid = vertices[vid]->rc->id;
                std::set<int> edgeids_to_remove;
                Edge* e_out =  vertices[vid]->outgoing[0];
                Vertex* end_v = e_out->end;
                edgeids_to_remove.insert(e_out->getId());
//                edgeids_to_remove.insert(e_out->rc->id)
                Edge* e_in = vertices[vid]->rc->outgoing[0]->rc;
                Vertex* start_v = e_in->start;
                edgeids_to_remove.insert(e_in->getId());
//                edgeids_to_remove.insert(e_out->rc->id);
                size_t overlap = vertices[vid]->k();
                size_t pref = e_in->getSeq().size();
                if (pref >= overlap )
                        pref -= overlap;
                else {
                    pref = 0;
                    std::cerr << " wrong overlap " << overlap << " " << e_in->getSeq().size() << endl;
                }
                Sequence new_seq = e_in->getSeq().Prefix(pref) + e_out->getSeq();
                string new_label = e_in->getLabel()+ "_"+ e_out->getLabel();
                addEdge(*start_v, *end_v, new_seq, 0, new_label);
                for (auto eid: edgeids_to_remove){
                    internalRemoveEdge(eid);
                }
                delete vertices[vid];
                delete vertices[rcid];
                vertices.erase(vid);
                vertices.erase(rcid);
            }
        }

        void deleteEdgeById(int eid){
            Vertex* start_v = edges[eid]->start;
            Vertex* end_v = edges[eid]->end;
            internalRemoveEdge(eid);
            compressVertex(start_v->id);
            compressVertex(end_v->id);
        }


        void deleteEdgeByLabel(std::string label){
            if (label_to_id.find(label) != label_to_id.end()) {
                int eid = label_to_id[label];
                deleteEdgeById(eid);
            }
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
            for(auto p: edges) {
                Edge *edge = p.second;
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
            for(auto p : vertices) {
                Vertex* vertex = p.second;

                if(vertex->inDeg() == 0 && vertex->outDeg() == 0 && vertex->seq <= !vertex->seq) {
                    res.addVertex(vertex->seq);
                }
            }
            res.checkConsistency();
            return std::move(res);
        }

        std::vector<Contig> getEdges(bool cut_overlaps) {
            std::unordered_map<Vertex *, size_t> cut;
            for(auto p : vertices) {
                Vertex* v = p.second;

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
            for(auto p: edges) {
                Edge *edge = p.second;
                if(edge->isCanonical()) {
                    size_t cut_left = edge->start->seq.size() * cut[edge->start];
                    size_t cut_right = edge->end->seq.size() * (1 - cut[edge->end]);
                    if(!cut_overlaps) {
                        cut_left = 0;
                        cut_right = 0;
                    }
                    if(cut_left + cut_right >= edge->size()) {
                        continue;
                    }
                    res.emplace_back(edge->getSeq().Subseq(cut_left, edge->size() - cut_right), itos(edge->getId()));
                    cnt++;
                }
            }
            return std::move(res);
        }

        void printEdges(const std::experimental::filesystem::path &f, bool cut_overlaps) {
            std::ofstream os;
            os.open(f);
            for(const Contig &contig : getEdges(cut_overlaps)) {
                os << ">" << contig.id << "\n" << contig.seq << "\n";
            }
            os.close();
        }

        void printDot(const std::experimental::filesystem::path &f) {
            std::ofstream os;
            os.open(f);
            os << "digraph {\nnodesep = 0.5;\n";
            std::unordered_map<Vertex *, int> vmap;
            for(auto p : vertices) {
                Vertex* vertex = p.second;
                os << vertex->id << " [label=\"" << vertex->seq.size() << "\" style=filled fillcolor=\"white\"]\n";
            }
            std::unordered_map<Edge *, std::string> eids;
            for (auto p : edges) {
                Edge *edge = p.second;
                os << "\"" << edge->start->id << "\" -> \"" << edge->end->id <<
                   "\" [label=\"" << edge->getId() << "(" << edge->size() << ")\" color = \"black\"]\n" ;

            }
            os << "}\n";
            os.close();
        }

        void printEdgeGFA(const std::experimental::filesystem::path &f, const std::vector<Vertex *> &component, bool labels = false) const {
            std::ofstream os;
            os.open(f);
            os << "H\tVN:Z:1.0" << std::endl;
            std::unordered_map<Edge *, std::string> eids;
            for(Vertex *v : component)
                for (Edge *edge : v->outgoing) {
                    if (edge->isCanonical()) {
                        if (labels) {
                            eids[edge] = edge->getLabel();
                            eids[edge->rc] = edge->getLabel();
                        } else {
                            eids[edge] = itos(edge->getId());
                            eids[edge->rc] = itos(edge->getId());

                        }

                        os << "S\t" << eids[edge] << "\t" << edge->getSeq() << "\n";
                    }
                }
            for (Vertex *vertex : component) {
                if(!vertex->isCanonical())
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

        void printEdgeGFA(const std::experimental::filesystem::path &f, bool labels = false) const {
            std::vector<Vertex* > component;
            for (auto p: vertices)
                component.push_back(p.second);
            printEdgeGFA(f, component, labels);
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
            std::vector<Vertex* > component;
            for (auto p: vertices)
                component.push_back(p.second);
            printVertexGFA(f, component);
        }

        std::vector<std::vector<Vertex *>> split() const {
            std::vector<std::vector<Vertex *>> res;
            std::unordered_set<Vertex *> visited;
            for(auto p : vertices) {
                Vertex * v = p.second;
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
