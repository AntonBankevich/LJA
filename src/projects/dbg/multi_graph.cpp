#include "multi_graph.hpp"

namespace multigraph {
    TransformingIterator<std::vector<Edge *>::iterator, Edge> Vertex::begin() {
        return TransformingIterator<std::vector<Edge *>::iterator, Edge>::DereferencingIterator(outgoing.begin(),
                                                                                                outgoing.end());
    }

    TransformingIterator<std::vector<Edge *>::iterator, Edge> Vertex::end() {
        return TransformingIterator<std::vector<Edge *>::iterator, Edge>::DereferencingIterator(outgoing.end(), outgoing.end());
    }

    TransformingIterator<std::vector<Edge *>::const_iterator, const Edge> Vertex::begin() const {
        return TransformingIterator<std::vector<Edge *>::const_iterator, const Edge>::DereferencingConstIterator(outgoing.begin(), outgoing.end());
    }

    TransformingIterator<std::vector<Edge *>::const_iterator, const Edge>
    Vertex::end() const {
        return TransformingIterator<std::vector<Edge *>::const_iterator, const Edge>::DereferencingConstIterator(outgoing.end(), outgoing.end());
    }

    VertexId MultiGraph::getVertexById(int id) {
        auto it = vertices_map.find(id);
        if(it == vertices_map.end())
            return {};
        else
            return it->second.getId();
    }

    ConstVertexId MultiGraph::getVertexById(int id) const {
        auto it = vertices_map.find(id);
        if(it == vertices_map.end())
            return {};
        else
            return it->second.getId();
    }

    EdgeId MultiGraph::getEdgeById(int id) {
        auto it = edges_map.find(id);
        if(it == edges_map.end())
            return {};
        else
            return it->second.getId();
    }

    ConstEdgeId MultiGraph::getEdgeById(int id) const {
        auto it = edges_map.find(id);
        if(it == edges_map.end())
            return {};
        else
            return it->second.getId();
    }

    IterableStorage<TransformingIterator<std::unordered_map<int, Vertex>::iterator, Vertex>> MultiGraph::vertices() {
        typedef std::unordered_map<int, Vertex>::iterator iterator;
        typedef Vertex value_type;
        std::function<value_type &(iterator::reference &)> transform =
                [](iterator::reference value) -> value_type & { return value.second;};
        TransformingIterator<iterator, value_type> begin(vertices_map.begin(), vertices_map.end(), transform);
        TransformingIterator<iterator, value_type> end(vertices_map.end(), vertices_map.end(), transform);
        return {begin, end};
    }

    IterableStorage<TransformingIterator<std::unordered_map<int, Edge>::iterator, Edge>> MultiGraph::edges() {
        typedef std::unordered_map<int, Edge>::iterator iterator;
        typedef Edge value_type;
        std::function<value_type &(iterator::reference &)> transform =
                [](iterator::reference value) -> value_type & { return value.second;};
        TransformingIterator<iterator, value_type> begin(edges_map.begin(), edges_map.end(), transform);
        TransformingIterator<iterator, value_type> end(edges_map.end(), edges_map.end(), transform);
        return {begin, end};
    }

    IterableStorage<TransformingIterator<std::unordered_map<int, Vertex>::const_iterator, const Vertex>>
    MultiGraph::vertices() const {
        typedef std::unordered_map<int, Vertex>::const_iterator iterator;
        typedef const Vertex value_type;
        std::function<value_type &(iterator::reference &)> transform =
                [](iterator::reference value) -> value_type & { return value.second;};
        TransformingIterator<iterator, value_type> begin(vertices_map.begin(), vertices_map.end(), transform);
        TransformingIterator<iterator, value_type> end(vertices_map.end(), vertices_map.end(), transform);
        return {begin, end};
    }

    IterableStorage<TransformingIterator<std::unordered_map<int, Edge>::const_iterator, const Edge>>
    MultiGraph::edges() const {
        typedef std::unordered_map<int, Edge>::const_iterator iterator;
        typedef const Edge value_type;
        std::function<value_type &(iterator::reference &)> transform =
                [](iterator::reference value) -> value_type & { return value.second;};
        TransformingIterator<iterator, value_type> begin(edges_map.begin(), edges_map.end(), transform);
        TransformingIterator<iterator, value_type> end(edges_map.end(), edges_map.end(), transform);
        return {begin, end};
    }


    void MultiGraph::printStats(std::ostream &os) {
        size_t tips = 0;
        size_t total_length = 0;
        size_t total_clean_length = 0;
        for(Vertex &v : vertices()) {
            for(Edge &edge : v) {
                total_length += edge.size();
                total_clean_length += edge.size() - v.seq.size();
            }
            if(v.outDeg() == 1 && v.inDeg() == 0) {
                tips += 2;
                total_clean_length += v.seq.size();
            }
        }
        os << "Vertices: " << vertices_map.size() << "\nEdges: " << edges_map.size() << "\n";
        os << "Total getEdge length: " << total_length << "\nTotal clean getEdge length: " << total_clean_length << "\n";
        os << "Number of tips: " << tips << "\n";
    }

    MultiGraph MultiGraphHelper::TransformToEdgeGraph(const MultiGraph &mg, size_t tip_size) {
        MultiGraph dbg;
        std::unordered_map<ConstEdgeId, VertexId> emap;
        for(const Vertex &v : mg.vertices()) {
            if(v.outDeg() == 0 || emap.find(v.begin()->getId()) != emap.end()) {
                continue;
            }
            Vertex &newv = dbg.addVertex(v.getSeq().Subseq(v.size() - v[0].overlap()));
            for(const Edge &edge : v) {
                const Vertex &right = edge.getFinish();
                for(const Edge &edge1 : right.rc()) {
                    emap[edge1.getId()] = newv.rc().getId();
                    emap[edge1.rc().getId()] = newv.getId();
                }
            }
        }
        for(const Vertex &v : mg.vertices()) {
            if(!v.isCanonical())
                continue;
            VertexId start;
            VertexId end;
            if(v.inDeg() == 0) {
                start = dbg.addVertex(v.getSeq().Subseq(0, std::min(tip_size, v.size() - 1))).getId();
            } else {
                start = emap[v.rc().begin()->getId()]->rc().getId();
            }
            if(v.outDeg() == 0) {
                end = dbg.addVertex(v.getSeq().Subseq(v.size() - std::min(tip_size, v.size() - 1))).getId();
            } else {
                end = emap[v.begin()->getId()];
            }
            dbg.addEdge(*start, *end, v.getSeq(), v.getId().innerId(), v.getLabel());
        }
        return std::move(dbg);
    }

    MultiGraph MultiGraphHelper::Delete(const MultiGraph &initial, const std::unordered_set<ConstEdgeId> &to_delete,
                                  const std::unordered_set<ConstVertexId> &to_delete_vertices) {
        MultiGraph res;
        std::unordered_map<ConstVertexId , VertexId> vmap;
        std::unordered_set<ConstEdgeId> visited;
        for(const Vertex &v : initial.vertices()) {
            if(to_delete_vertices.find(v.getId()) != to_delete_vertices.end())
                continue;
            if(vmap.find(v.getId()) != vmap.end())
                continue;
            vmap[v.getId()] = res.addVertex(v.getSeq(), v.getId().innerId()).getId();
            vmap[v.rc().getId()] = vmap[v.getId()]->rc().getId();
        }
        for(const Edge &edge : initial.edges()) {
            if(!edge.isCanonical() || to_delete.find(edge.getId()) != to_delete.end())
                continue;
            if(to_delete_vertices.find(edge.getStart().getId()) == to_delete_vertices.end() ||
               to_delete_vertices.find(edge.getFinish().getId()) == to_delete_vertices.end())
                res.addEdge(*vmap[edge.getStart().getId()], *vmap[edge.getFinish().getId()], edge.getSeq(), edge.getId().innerId(), edge.getLabel());
        }
        return std::move(res);
    }

    void MultiGraph::checkConsistency() const {
        std::unordered_set<ConstEdgeId> eset;
        std::unordered_set<ConstVertexId> vset;
        for(const Edge &edge: edges()) {
            eset.emplace(edge.getId());
            VERIFY(edge.rc().getStart() == edge.getFinish().rc());
            VERIFY(edge.rc().rc() == edge);
        }
        for(const Edge &edge: edges()) {
            VERIFY(eset.find(edge.rc().getId()) != eset.end());
        }
        for(const Vertex &v : vertices()) {
            vset.emplace(v.getId());
            VERIFY(v.rc().rc() == v);
            for(const Edge &edge : v) {
                VERIFY(eset.find(edge.getId()) != eset.end());
                VERIFY(edge.getStart() == v);
            }
        }
        for(const Vertex &v : vertices()) {
            VERIFY(vset.find(v.getId()) != vset.end());
        }
    }

    Vertex &MultiGraph::addVertex(const Sequence &seq, int id, std::string label) {
        if (id == 0) {
            if(seq <= !seq)
                id = maxVId + 1;
            else
                id = -maxVId - 1;
        }
        maxVId = std::max(std::abs(id), maxVId);
        VertexId res = vertices_map.emplace(std::piecewise_construct, std::forward_as_tuple(id),
                                            std::forward_as_tuple(seq, id, label)).first->second.getId();
        VertexId rc = res;
        if(seq != !seq) {
            rc = vertices_map.emplace(std::piecewise_construct, std::forward_as_tuple(-id),
                                      std::forward_as_tuple(!seq, -id, label)).first->second.getId();
        }
        res->setRC(*rc);
        rc->setRC(*res);
        return *res;
    }

    Edge &MultiGraph::addEdge(Vertex &from, Vertex &to, Sequence seq, int id, std::string label) {
        if(id == 0) {
            if(seq <= !seq)
                id = maxEId + 1;
            else
                id = -maxEId - 1;
        }
        maxEId = std::max(std::abs(id), maxEId);
        if(!(seq <= !seq)) {
            id = -(std::abs(id));
        } else {
            id = std::abs(id);
        }

        EdgeId res = edges_map.emplace(std::piecewise_construct, std::forward_as_tuple(id),
                                       std::forward_as_tuple(from, to, std::move(seq), id, label)).first->second.getId();
        res->getStart().addOutgoing(*res);
        EdgeId rc = res;
        if(res->getSeq() != !res->getSeq()) {
            rc = edges_map.emplace(std::piecewise_construct, std::forward_as_tuple(-id),
                                   std::forward_as_tuple(to.rc(), from.rc(), !res->getSeq(), -id, res->getReverseLabel())).first->second.getId();
            rc->getStart().addOutgoing(*rc);
        }
        res->setRC(*rc);
        rc->setRC(*res);
        return *res;
    }

    void MultiGraph::internalRemoveEdge(Edge &edge) {
        edge.getStart().removeOutgoing(edge);
        if(edge != edge.rc()) {
            edge.rc().getStart().removeOutgoing(edge.rc());
            edges_map.erase(edge.rc().getId().innerId());
        }
        edges_map.erase(edge.getId().innerId());
    }

    void MultiGraph::internalRemoveIsolatedVertex(Vertex &vertex) {
        VERIFY(vertex.outDeg() == 0 && vertex.inDeg() == 0);
        if(vertex != vertex.rc())
            vertices_map.erase(vertex.getId().innerId());
    }

    deleted_edges_map MultiGraph::attemptCompressVertex(Vertex &v) {
//Only 1-in-1-out vertices are compressed
        if(v.outDeg()!= 1 || v.inDeg() != 1) {
            return {};
        }
        Edge &edge1 = v.rc().begin()->rc();
        Edge &edge2 = *v.begin();
//Do not compress loops
        if(edge1 == edge2)
            return {};
//Do not compress rc loops
        if(edge1 == edge1.rc() && edge2==edge2.rc()) {
            return {};
        }
//Rc loops are always loops, so this condition should never fail.
        VERIFY_MSG(edge1 != edge2.rc(), "Reverse-complement loop is not a loop. This should never happen");
        VERIFY(edge1.size() >= v.size());
        VERIFY(edge2.size() >= v.size());
        Sequence new_seq = edge1.getSeq() + edge2.getSeq().Subseq(v.size());
        VertexId new_start = edge1.getStart().getId();
        VertexId new_end = edge2.getFinish().getId();
        std::string new_label = ((edge1.isCanonical()?edge1.getLabel(): edge1.getReverseLabel() + "_" +
                                                                        (edge2.isCanonical()?edge2.getLabel(): edge2.getReverseLabel())));
        std::vector<std::string> path = {edge1.getLabel(), edge2.getLabel()};

        if(edge1 == edge1.rc()) {
            new_seq = edge2.rc().getSeq().Prefix(edge2.size() - v.size()) + new_seq;
            new_start = edge2.rc().getStart().getId();
            new_end = edge2.getFinish().getId();
            new_label = edge2.rc().getLabel() + "_" + edge1.getLabel() + "_" + edge2.getLabel();
            path = {edge2.rc().getLabel(), edge1.getLabel(), edge2.getLabel()};
        } else if(edge2 == edge2.rc()) {
            new_seq = new_seq + edge1.rc().getSeq().Subseq(v.size());
            new_start = edge1.getStart().getId();
            new_end = edge1.rc().getFinish().getId();
            new_label = edge1.getLabel() + "_" + edge2.getLabel() + "_" + edge1.rc().getLabel();
            path = {edge1.getLabel(), edge2.getLabel(), edge1.rc().getLabel()};
        }

        deleted_edges_map result_map;
        result_map[new_label] = path;
        addEdge(*new_start, *new_end, new_seq, 0, new_label);
        internalRemoveEdge(edge1);
        internalRemoveEdge(edge2);
        internalRemoveIsolatedVertex(v);
        return result_map;
    }

    std::vector<EdgeId> MultiGraphHelper::uniquePathForward(Edge &edge) {
        std::vector<EdgeId> res = {edge.getId()};
        VertexId cur = edge.getFinish().getId();
        while(cur != edge.getStart().getId() && cur->inDeg() == 1 && cur->outDeg() == 1) {
            res.emplace_back(cur->begin()->getId());
            cur = res.back()->getFinish().getId();
        }
        return std::move(res);
    }

    std::vector<ConstEdgeId> MultiGraphHelper::uniquePathForward(const Edge &edge) {
        std::vector<ConstEdgeId> res = {edge.getId()};
        ConstVertexId cur = edge.getFinish().getId();
        while(cur != edge.getStart().getId() && cur->inDeg() == 1 && cur->outDeg() == 1) {
            res.emplace_back(cur->begin()->getId());
            cur = res.back()->getFinish().getId();
        }
        return std::move(res);
    }

    std::vector<ConstEdgeId> MultiGraphHelper::uniquePath(const Edge &edge) {
        std::vector<ConstEdgeId> path = uniquePathForward(edge.rc());
        return uniquePathForward(path.back()->rc());
    }

    std::vector<EdgeId> MultiGraphHelper::uniquePath(Edge &edge) {
        std::vector<EdgeId> path = uniquePathForward(edge.rc());
        return uniquePathForward(path.back()->rc());
    }

    MultiGraph MultiGraphHelper::MergeAllPaths(const MultiGraph &mg, bool verbose) {
        MultiGraph res;
        std::unordered_set<ConstEdgeId> used;
        std::unordered_map<ConstVertexId, VertexId> old_to_new;
        for(const Edge &edge: mg.edges()) {
            if(used.find(edge.getId()) != used.end())
                continue;
            std::vector<ConstEdgeId> tmp = MultiGraphHelper::uniquePath(edge);
            std::vector<std::vector<ConstEdgeId>> paths_to_add = {{}};
            for(ConstEdgeId e : tmp) {
                paths_to_add.back().emplace_back(e);
                if(e->rc() == edge) {
                    paths_to_add.emplace_back(std::vector<ConstEdgeId>());
                }
                used.emplace(e);
                used.emplace(e->rc().getId());
            }
            VERIFY(paths_to_add.size() <= 2);
            if(paths_to_add.back().empty())
                paths_to_add.pop_back();
            VERIFY(!paths_to_add.empty());
            for(std::vector<ConstEdgeId> &path : paths_to_add) {
                ConstVertexId old_start = path.front()->getStart().getId();
                ConstVertexId old_end = path.back()->getFinish().getId();
                VertexId &new_start = old_to_new[old_start];
                if(!new_start.valid()) {
                    new_start = res.addVertex(old_start->getSeq()).getId();
                    old_to_new[old_start->rc().getId()] = new_start->rc().getId();
                }
                VertexId &new_end = old_to_new[old_end];
                if(!new_end.valid()) {
                    new_end = res.addVertex(old_end->getSeq()).getId();
                    old_to_new[old_end->rc().getId()] = new_end->rc().getId();
                }
                SequenceBuilder sb;
                sb.append(path.front()->getSeq());
                for(size_t i = 1; i < path.size(); i++) {
                    sb.append(path[i]->getSeq().Subseq(path[i]->getStart().size()));
                }
                Edge &new_edge = res.addEdge(*new_start, *new_end, sb.BuildSequence());
                if(verbose) {
                    std::cout << "New getEdge " << new_edge.getId() << " consists of old edges: ";
                    for(auto e : path) {
                        std::cout << e->getId() << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }
        for(const Vertex &vertex : mg.vertices()) {
            if(vertex.inDeg() == 0 && vertex.outDeg() == 0 && vertex.isCanonical()) {
                res.addVertex(vertex.getSeq());
            }
        }
        res.checkConsistency();
        return std::move(res);
    }

    std::vector<Contig> MultiGraphHelper::extractContigs(const MultiGraph &mg, bool cut_overlaps) {
        std::unordered_map<ConstVertexId, size_t> cut;
        for(const Vertex &v : mg.vertices()) {
            if(v.isCanonical()) {
                if(v.outDeg() == 1) {
                    cut[v.getId()] = 0;
                } else {
                    cut[v.getId()] = 1;
                }
                cut[v.rc().getId()] = 1 - cut[v.getId()];
            }
        }
        std::vector<Contig> res;
        size_t cnt = 1;
        for(const Edge &edge: mg.edges()) {
            if(edge.isCanonical()) {
                size_t cut_left = edge.getStart().size() * cut[edge.getStart().getId()];
                size_t cut_right = edge.getFinish().size() * (1 - cut[edge.getFinish().getId()]);
                if(!cut_overlaps) {
                    cut_left = 0;
                    cut_right = 0;
                }
                if(cut_left + cut_right >= edge.size()) {
                    continue;
                }
                res.emplace_back(edge.getSeq().Subseq(cut_left, edge.size() - cut_right), itos(edge.getId().innerId()));
                cnt++;
            }
        }
        return std::move(res);
    }

    void MultiGraphHelper::printExtractedContigs(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool cut_overlaps) {
        std::ofstream os;
        os.open(f);
        for(const Contig &contig : extractContigs(mg, cut_overlaps)) {
            os << ">" << contig.getInnerId() << "\n" << contig.getSeq() << "\n";
        }
        os.close();
    }

    void MultiGraphHelper::printDot(const MultiGraph &mg, const std::experimental::filesystem::path &f) {
        std::ofstream os;
        os.open(f);
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_map<const Vertex *, int> vmap;
        for(const Vertex &vertex : mg.vertices()) {
            os << vertex.getId() << " [label=\"" << vertex.size() << "\" style=filled fillcolor=\"white\"]\n";
        }
        std::unordered_map<EdgeId, std::string> eids;
        for (const Edge &edge : mg.edges()) {
            os << "\"" << edge.getStart().getId() << "\" -> \"" << edge.getFinish().getId() <<
               "\" [label=\"" << edge.getId() << "(" << edge.size() << ")\" color = \"black\"]\n" ;

        }
        os << "}\n";
        os.close();
    }

    void
    MultiGraphHelper::printEdgeGFA(const std::experimental::filesystem::path &f, const std::vector<ConstVertexId> &component,
                             bool labels) {
        std::ofstream os;
        os.open(f);
        os << "H\tVN:Z:1.0" << std::endl;
        std::unordered_map<ConstEdgeId, std::string> eids;
        for(ConstVertexId v : component) {
            for (const Edge &edge : *v) {
                if (edge.isCanonical()) {
                    if (labels) {
                        eids[edge.getId()] = edge.getLabel();
                        eids[edge.rc().getId()] = edge.getLabel();
                    } else {
                        eids[edge.getId()] = itos(edge.getId().innerId());
                        eids[edge.rc().getId()] = itos(edge.getId().innerId());
                    }
                    os << "S\t" << eids[edge.getId()] << "\t" << edge.getSeq() << "\n";
                }
            }
        }
        for (ConstVertexId vertex : component) {
            if(!vertex->isCanonical())
                continue;
            for (const Edge &out_edge : *vertex) {
                std::string outid = eids[out_edge.getId()];
                bool outsign = out_edge.isCanonical();
                for (const Edge &inc_edge : vertex->rc()) {
                    std::string incid = eids[inc_edge.getId()];
                    bool incsign = inc_edge.rc().isCanonical();
                    os << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                       << (outsign ? "+" : "-") << "\t" << vertex->getSeq().size() << "M" << "\n";
                }
            }
        }
        os.close();
    }

    void MultiGraphHelper::printEdgeGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool labels) {
        std::vector<ConstVertexId> component;
        for (const Vertex &vertex: mg.vertices()) {
            component.push_back(vertex.getId());
        }
        printEdgeGFA(f, component, labels);
    }

    void MultiGraphHelper::printVertexGFA(const std::experimental::filesystem::path &f,
                                    const std::vector<ConstVertexId> &component) {
        std::ofstream os;
        os.open(f);
        os << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 1;
        std::unordered_map<ConstVertexId , std::string> vids;
        for(ConstVertexId v : component)
            if(v->isCanonical()) {
                vids[v] = itos(v.innerId());
                vids[v->rc().getId()] = itos(v.innerId());
                os << "S\t" << vids[v] << "\t" << v->getSeq() << "\n";
                cnt++;
            }
        for(ConstVertexId v : component)
            for (const Edge &edge : *v) {
                if (edge.isCanonical()) {
                    VERIFY(edge.size() < edge.getStart().size() + edge.getFinish().size());
                    bool incsign = v->isCanonical();
                    bool outsign = edge.getFinish().isCanonical();
                    os << "L\t" << vids[v] << "\t" << (incsign ? "+" : "-") << "\t"
                       << vids[edge.getFinish().getId()] << "\t" << (outsign ? "+" : "-") << "\t"
                       << (v->size() + edge.getFinish().size() - edge.size()) << "M" << "\n";
                }
            }
        os.close();
    }

    void MultiGraphHelper::printVertexGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f) {
        std::vector<ConstVertexId> component;
        for (const Vertex &vertex: mg.vertices())
            component.push_back(vertex.getId());
        printVertexGFA(f, component);
    }

    std::vector<std::vector<ConstVertexId>> MultiGraphHelper::split(const MultiGraph &mg) {
        std::vector<std::vector<ConstVertexId>> res;
        std::unordered_set<ConstVertexId> visited;
        for(const Vertex &start: mg.vertices()) {
            if(visited.find(start.getId()) != visited.end())
                continue;
            std::vector<ConstVertexId> stack = {start.getId()};
            res.emplace_back(std::vector<ConstVertexId>());
            while(!stack.empty()) {
                ConstVertexId v = stack.back();
                stack.pop_back();
                if(visited.find(v) != visited.end())
                    continue;
                visited.emplace(v);
                visited.emplace(v->rc().getId());
                res.back().emplace_back(v);
                res.back().emplace_back(v->rc().getId());
                for(const Edge &e : *v)
                    stack.emplace_back(e.getFinish().getId());
                for(const Edge &e : v->rc())
                    stack.emplace_back(e.getFinish().getId());
            }
        }
        return std::move(res);
    }

//    deleted_edges_map MultiGraph::deleteAndCompress(Edge &edge) {
//        Vertex &start = edge.start();
//        Vertex &end = edge.end();
//        internalRemoveEdge(edge);
//
//        deleted_edges_map result;
//        if(end != start && end != start.rc()) {
//            result = attemptCompressVertex(end);
//        }
//        deleted_edges_map res2 = attemptCompressVertex(end);
//        result.insert(res2.begin(), res2.end());
//
////this compression may contain edges that result from first one, so additional ugly processing required.
////            auto comp_res = attemprCompressVertex(end_v->id);
////            for (auto p: comp_res) {
////                std::vector<std::string> patched_old;
////                for (auto comp_edge: p.second) {
////                    if (res.find(comp_edge) != res.end()) {
////                        patched_old.insert(patched_old.end(), res[comp_edge].begin(), res[comp_edge].end());
////                    } else
////                        patched_old.push_back(comp_edge);
////                }
////                res[p.first] = patched_old;
////            }
//        return result;
//    }
    MultiGraph MultiGraphHelper::LoadGFA(const std::experimental::filesystem::path &gfa_file, bool int_ids) {
        MultiGraph res;
        std::ifstream is;
        is.open(gfa_file);
        std::unordered_map<std::string, VertexId> vmap;
        for(std::string line; getline(is, line); ) {
            std::vector<std::string> tokens = ::split(line);
            if(tokens[0] == "S") {
                std::string name = tokens[1];
                Vertex &newV = int_ids ? res.addVertex(Sequence(tokens[2]), std::stoi(name), name): res.addVertex(Sequence(tokens[2]));

                VERIFY(vmap.find(name) == vmap.end());
                vmap[name] = newV.getId();
            } else if(tokens[0] == "L") {
                VertexId v1 = vmap[tokens[1]];
                VertexId v2 = vmap[tokens[3]];
                if(tokens[2] == "-")
                    v1 = v1->rc().getId();
                if(tokens[4] == "-")
                    v2 = v2->rc().getId();
                size_t overlap = std::stoull(tokens[5].substr(0, tokens[5].size() - 1));
                if(v1->getSeq().Subseq(v1->getSeq().size() - overlap) != v2->getSeq().Subseq(0, overlap)) {
                    v1 = v1->rc().getId();
                }
                VERIFY(v1->getSeq().Subseq(v1->getSeq().size() - overlap) == v2->getSeq().Subseq(0, overlap));
                res.addEdge(*v1, *v2, v1->getSeq() + v2->getSeq().Subseq(overlap));
            }
        }
        is.close();
        return std::move(res);
    }

    Sequence Edge::getSeq() const {
        if(seq.empty())
            return _start->getSeq() + _end->getSeq().Subseq(sz);
        else
            return seq;
    }

    std::string Edge::getReverseLabel() const {
        if(label.empty())
            return {};
        std::vector<std::string> tokens = ::split(label, "_");
        std::string res;
        for (size_t i = tokens.size() -1; i > 0; i --)
            res += tokens[i] + "_";
        res += tokens[0];
        return res;
    }

    size_t Edge::overlap() const {
        VERIFY(_start->size() + _end->size() > sz);
        return _start->size() + _end->size() - sz;
    }

    bool Edge::isCanonical() const {
        VERIFY(canonical == (id > 0));
        return canonical;
    }

    bool Edge::isTip() const {
        return (_start->inDeg() == 0  || _end->outDeg() == 0);
    }

    bool Edge::isSimpleBridge() {
        if (isTip())
            return false;
        for (Edge &alt_e: getStart()) {
            if (alt_e != *this and !alt_e.isTip())
                return false;
        }
        for (Edge &alt_e: rc().getStart()) {
            if (alt_e != rc() and !alt_e.isTip())
                return false;
        }
        return true;
    }

    bool Edge::operator==(const Edge &other) const {
        if(id == other.id) {
            VERIFY(this == &other);
            return true;
        } else {
            return false;
        }
    }

    bool Vertex::operator==(const Vertex &other) const {
        if(id == other.id) {
            VERIFY_MSG(this == &other, "Different vertices have the same id");
            return true;
        } else {
            return false;
        }
    }

    double EdgeData::getCoverage() const {
        return double(cov) / edge->truncSize();
    }
}