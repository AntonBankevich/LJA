#include <common/disjoint_sets.hpp>
#include "multi_graph.hpp"

namespace multigraph {

    MultiGraph MultiGraphHelper::TransformToEdgeGraph(const MultiGraph &mg, size_t tip_size) {
        MultiGraph dbg;
        std::unordered_map<ConstEdgeId, VertexId> emap;
        for (const MGVertex &v: mg.vertices()) {
            if (v.outDeg() == 0 || emap.find(v.front().getId()) != emap.end()) {
                continue;
            }
            MGVertex &newv = dbg.addVertex(v.getSeq().Subseq(v.size() - v.front().overlapSize()), MGVertexData(""));
            for (const MGEdge &edge: v) {
                const MGVertex &right = edge.getFinish();
                for (const MGEdge &edge1: right.rc()) {
                    emap[edge1.getId()] = newv.rc().getId();
                    emap[edge1.rc().getId()] = newv.getId();
                }
            }
        }
        for (const MGVertex &v: mg.verticesUnique()) {
            VertexId start;
            VertexId end;
            if (v.inDeg() == 0) {
                start = dbg.addVertex(v.getSeq().Subseq(0, std::min(tip_size, v.size() - 1))).getId();
            } else {
                start = emap[v.rc().begin()->getId()]->rc().getId();
            }
            if (v.outDeg() == 0) {
                end = dbg.addVertex(v.getSeq().Subseq(v.size() - std::min(tip_size, v.size() - 1))).getId();
            } else {
                end = emap[v.begin()->getId()];
            }
            start->addEdgeLockFree(*end, v.getSeq(), MGEdgeData(v.getLabel()));
        }
        return std::move(dbg);
    }

    MultiGraph MultiGraphHelper::Delete(const MultiGraph &initial, const std::unordered_set<ConstEdgeId> &to_delete,
                                        const std::unordered_set<ConstVertexId> &to_delete_vertices) {
        MultiGraph res;
        std::unordered_map<ConstVertexId, VertexId> vmap;
        std::unordered_set<ConstEdgeId> visited;
        for (const MGVertex &v: initial.verticesUnique()) {
            if (to_delete_vertices.find(v.getId()) != to_delete_vertices.end())
                continue;
            vmap[v.getId()] = res.addVertex(v).getId();
            vmap[v.rc().getId()] = vmap[v.getId()]->rc().getId();
        }
        for (const MGEdge &edge: initial.edgesUnique()) {
            if (to_delete.find(edge.getId()) != to_delete.end())
                continue;
            if (to_delete_vertices.find(edge.getStart().getId()) == to_delete_vertices.end() ||
                to_delete_vertices.find(edge.getFinish().getId()) == to_delete_vertices.end())
                vmap[edge.getStart().getId()]->addEdgeLockFree(*vmap[edge.getFinish().getId()], edge.getSeq(), edge);
        }
        return std::move(res);
    }

    std::vector<EdgeId> MultiGraphHelper::uniquePathForward(MGEdge &edge) {
        std::vector<EdgeId> res = {edge.getId()};
        VertexId cur = edge.getFinish().getId();
        while (cur != edge.getStart().getId() && cur->inDeg() == 1 && cur->outDeg() == 1) {
            res.emplace_back(cur->begin()->getId());
            cur = res.back()->getFinish().getId();
        }
        return std::move(res);
    }

    std::vector<ConstEdgeId> MultiGraphHelper::uniquePathForward(const MGEdge &edge) {
        std::vector<ConstEdgeId> res = {edge.getId()};
        ConstVertexId cur = edge.getFinish().getId();
        while (cur != edge.getStart().getId() && cur->inDeg() == 1 && cur->outDeg() == 1) {
            res.emplace_back(cur->begin()->getId());
            cur = res.back()->getFinish().getId();
        }
        return std::move(res);
    }

    std::vector<ConstEdgeId> MultiGraphHelper::uniquePath(const MGEdge &edge) {
        std::vector<ConstEdgeId> path = uniquePathForward(edge.rc());
        return uniquePathForward(path.back()->rc());
    }

    std::vector<EdgeId> MultiGraphHelper::uniquePath(MGEdge &edge) {
        std::vector<EdgeId> path = uniquePathForward(edge.rc());
        return uniquePathForward(path.back()->rc());
    }

//    MultiGraph MultiGraphHelper::MergeAllPaths(const MultiGraph &mg, bool verbose) {
//        MultiGraph res;
//        std::unordered_set<ConstEdgeId> used;
//        std::unordered_map<ConstVertexId, VertexId> old_to_new;
//        for(const MGEdge &edge: mg.edges()) {
//            if(used.find(edge.getId()) != used.end())
//                continue;
//            std::vector<ConstEdgeId> tmp = MultiGraphHelper::uniquePath(edge);
//            std::vector<std::vector<ConstEdgeId>> paths_to_add = {{}};
//            for(ConstEdgeId e : tmp) {
//                paths_to_add.back().emplace_back(e);
//                if(e->rc() == edge) {
//                    paths_to_add.emplace_back(std::vector<ConstEdgeId>());
//                }
//                used.emplace(e);
//                used.emplace(e->rc().getId());
//            }
//            VERIFY(paths_to_add.size() <= 2);
//            if(paths_to_add.back().empty())
//                paths_to_add.pop_back();
//            VERIFY(!paths_to_add.empty());
//            for(std::vector<ConstEdgeId> &path : paths_to_add) {
//                ConstVertexId old_start = path.front()->getStart().getId();
//                ConstVertexId old_end = path.back()->getFinish().getId();
//                VertexId &new_start = old_to_new[old_start];
//                if(!new_start.valid()) {
//                    new_start = res.addVertex(old_start->getSeq()).getId();
//                    old_to_new[old_start->rc().getId()] = new_start->rc().getId();
//                }
//                VertexId &new_end = old_to_new[old_end];
//                if(!new_end.valid()) {
//                    new_end = res.addVertex(old_end->getSeq()).getId();
//                    old_to_new[old_end->rc().getId()] = new_end->rc().getId();
//                }
//                SequenceBuilder sb;
//                sb.append(path.front()->getSeq());
//                for(size_t i = 1; i < path.size(); i++) {
//                    sb.append(path[i]->getSeq().Subseq(path[i]->getStart().size()));
//                }
//                MGEdge &new_edge = new_start->addEdge(*new_end, sb.BuildSequence());
//                if(verbose) {
//                    std::cout << "New getEdge " << new_edge.getId() << " consists of old edges: ";
//                    for(auto e : path) {
//                        std::cout << e->getId() << " ";
//                    }
//                    std::cout << std::endl;
//                }
//            }
//        }
//        for(const MGVertex &vertex : mg.vertices()) {
//            if(vertex.inDeg() == 0 && vertex.outDeg() == 0 && vertex.isCanonical()) {
//                res.addVertex(vertex.getSeq());
//            }
//        }
//        MultiGraphHelper::checkConsistency(res);
//        return std::move(res);
//    }

    std::vector<Contig> MultiGraphHelper::extractContigs(const MultiGraph &mg, bool cut_overlaps) {
        std::unordered_map<ConstVertexId, size_t> cut;
        for (const MGVertex &v: mg.vertices()) {
            if (v.isCanonical()) {
                if (v.outDeg() == 1) {
                    cut[v.getId()] = 0;
                } else {
                    cut[v.getId()] = 1;
                }
                cut[v.rc().getId()] = 1 - cut[v.getId()];
            }
        }
        std::vector<Contig> res;
        size_t cnt = 1;
        for (const MGEdge &edge: mg.edges()) {
            if (edge.isCanonical()) {
                size_t cut_left = edge.getStart().size() * cut[edge.getStart().getId()];
                size_t cut_right = edge.getFinish().size() * (1 - cut[edge.getFinish().getId()]);
                if (!cut_overlaps) {
                    cut_left = 0;
                    cut_right = 0;
                }
                if (cut_left + cut_right >= edge.fullSize()) {
                    continue;
                }
                res.emplace_back(edge.getSeq().Subseq(cut_left, edge.fullSize() - cut_right),
                                 "E" + edge.getId().innerId().str());
                cnt++;
            }
        }
        return std::move(res);
    }

    void MultiGraphHelper::printExtractedContigs(const MultiGraph &mg, const std::experimental::filesystem::path &f,
                                                 bool cut_overlaps) {
        std::ofstream os;
        os.open(f);
        for (const Contig &contig: extractContigs(mg, cut_overlaps)) {
            os << ">" << contig.getInnerId() << "\n" << contig.getSeq() << "\n";
        }
        os.close();
    }

    void MultiGraphHelper::printDot(const MultiGraph &mg, const std::experimental::filesystem::path &f) {
        std::ofstream os;
        os.open(f);
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_map<const MGVertex *, int> vmap;
        for (const MGVertex &vertex: mg.vertices()) {
            os << vertex.getId() << " [label=\"" << vertex.size() << "\" style=filled fillcolor=\"white\"]\n";
        }
        std::unordered_map<EdgeId, std::string> eids;
        for (const MGEdge &edge: mg.edges()) {
            os << "\"" << edge.getStart().getId() << "\" -> \"" << edge.getFinish().getId() <<
               "\" [label=\"" << edge.getId() << "(" << edge.fullSize() << ")\" color = \"black\"]\n";

        }
        os << "}\n";
        os.close();
    }

    void
    MultiGraphHelper::printEdgeGFA(const std::experimental::filesystem::path &f,
                                   const std::vector<ConstVertexId> &component,
                                   bool labels) {
        std::ofstream os;
        os.open(f);
        os << "H\tVN:Z:1.0" << std::endl;
        std::unordered_map<ConstEdgeId, std::string> eids;
        for (ConstVertexId v: component) {
            for (const MGEdge &edge: *v) {
                if (edge.isCanonical()) {
                    if (labels) {
                        eids[edge.getId()] = edge.getLabel();
                        eids[edge.rc().getId()] = edge.getLabel();
                    } else {
                        eids[edge.getId()] = edge.getId().innerId().str();
                        eids[edge.rc().getId()] = edge.getId().innerId().str();
                    }
                    os << "S\t" << eids[edge.getId()] << "\t" << edge.getSeq() << "\n";
                }
            }
        }
        for (ConstVertexId vertex: component) {
            if (!vertex->isCanonical())
                continue;
            for (const MGEdge &out_edge: *vertex) {
                std::string outid = eids[out_edge.getId()];
                bool outsign = out_edge.isCanonical();
                for (const MGEdge &inc_edge: vertex->rc()) {
                    std::string incid = eids[inc_edge.getId()];
                    bool incsign = inc_edge.rc().isCanonical();
                    os << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                       << (outsign ? "+" : "-") << "\t" << vertex->getSeq().size() << "M" << "\n";
                }
            }
        }
        os.close();
    }

    void
    MultiGraphHelper::
    printEdgeGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool labels) {
        std::vector<ConstVertexId> component;
        for (const MGVertex &vertex: mg.vertices()) {
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
        std::unordered_map<ConstVertexId, std::string> vids;
        for (ConstVertexId v: component)
            if (v->isCanonical()) {
                vids[v] = itos(v.innerId());
                vids[v->rc().getId()] = itos(v.innerId());
                os << "S\t" << vids[v] << "\t" << v->getSeq() << "\n";
                cnt++;
            }
        for (ConstVertexId v: component)
            for (const MGEdge &edge: *v) {
                if (edge.isCanonical()) {
                    VERIFY(edge.fullSize() < edge.getStart().size() + edge.getFinish().size());
                    bool incsign = v->isCanonical();
                    bool outsign = edge.getFinish().isCanonical();
                    os << "L\t" << vids[v] << "\t" << (incsign ? "+" : "-") << "\t"
                       << vids[edge.getFinish().getId()] << "\t" << (outsign ? "+" : "-") << "\t"
                       << (v->size() + edge.getFinish().size() - edge.fullSize()) << "M" << "\n";
                }
            }
        os.close();
    }

    void MultiGraphHelper::printVertexGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f) {
        std::vector<ConstVertexId> component;
        for (const MGVertex &vertex: mg.vertices())
            component.push_back(vertex.getId());
        printVertexGFA(f, component);
    }

    std::vector<std::vector<ConstVertexId>> MultiGraphHelper::split(const MultiGraph &mg) {
        std::vector<std::vector<ConstVertexId>> res;
        std::unordered_set<ConstVertexId> visited;
        for (const MGVertex &start: mg.vertices()) {
            if (visited.find(start.getId()) != visited.end())
                continue;
            std::vector<ConstVertexId> stack = {start.getId()};
            res.emplace_back(std::vector<ConstVertexId>());
            while (!stack.empty()) {
                ConstVertexId v = stack.back();
                stack.pop_back();
                if (visited.find(v) != visited.end())
                    continue;
                visited.emplace(v);
                visited.emplace(v->rc().getId());
                res.back().emplace_back(v);
                res.back().emplace_back(v->rc().getId());
                for (const MGEdge &e: *v)
                    stack.emplace_back(e.getFinish().getId());
                for (const MGEdge &e: v->rc())
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
        for (std::string line; getline(is, line);) {
            std::vector<std::string> tokens = ::split(line);
            if (tokens[0] == "S") {
                std::string name = tokens[1];
		Sequence seq(tokens[2]);
		MGVertex &newV = int_ids ? res.addVertex(seq, {name}, seq.isCanonical() ? std::stoi(name) : -std::stoi(name)) : res.addVertex(seq);

                vmap[name] = newV.getId();
            } else if (tokens[0] == "L") {
                VertexId v1 = vmap[tokens[1]];
                VertexId v2 = vmap[tokens[3]];
                if (tokens[2] == "-")
                    v1 = v1->rc().getId();
                if (tokens[4] == "-")
                    v2 = v2->rc().getId();
                size_t overlap = std::stoull(tokens[5].substr(0, tokens[5].size() - 1));
                if (v1->getSeq().Subseq(v1->getSeq().size() - overlap) != v2->getSeq().Subseq(0, overlap)) {
                    v1 = v1->rc().getId();
                }
                VERIFY(v1->getSeq().Subseq(v1->getSeq().size() - overlap) == v2->getSeq().Subseq(0, overlap));
                v1->addEdge(*v2, v1->getSeq() + v2->getSeq().Subseq(overlap), MGEdgeData(""));
            }
        }
        is.close();
        return std::move(res);
    }

    struct GFAVertexRecord {
        std::string edgeId;
        bool start;
        bool rc;
        GFAVertexRecord() : start(false), rc(false) {}
        GFAVertexRecord(std::string edgeId, bool start, bool rc) : edgeId(std::move(edgeId)), start(start), rc(rc) {}
        GFAVertexRecord RC() const {
            return {edgeId, start, !rc};
        }
        GFAVertexRecord(const GFAVertexRecord &other) = default;
        GFAVertexRecord(GFAVertexRecord &&other) = default;
        GFAVertexRecord &operator=(const GFAVertexRecord &other) = default;
        GFAVertexRecord &operator=(GFAVertexRecord &&other) = default;
        bool operator==(const GFAVertexRecord &other) const {
            return edgeId == other.edgeId && start == other.start && rc == other.rc;
        }
    };
}
template<>
struct std::hash<multigraph::GFAVertexRecord> {
    size_t operator()(const multigraph::GFAVertexRecord &rec) const {
        return std::hash<std::string>()(rec.edgeId) + size_t(rec.rc) + (size_t(rec.start) << 16);
    }
};
namespace multigraph {
    MultiGraph MultiGraphHelper::LoadEdgeGFA(const std::experimental::filesystem::path &gfa_file, size_t K) {
        std::ifstream is;
        is.open(gfa_file);
        DisjointSet<GFAVertexRecord> vertices;
        std::unordered_map<GFAVertexRecord, VertexId> vertexMap;
        std::unordered_map<GFAVertexRecord, size_t> vertexLengths;
        std::vector<std::tuple<Sequence, std::string, Edge::id_type, Edge::id_type>> edges;
        size_t bad_ids = 0;
        for(std::string line; getline(is, line); ) {
            std::vector<std::string> tokens = ::split(line);
            if(tokens[0] == "S") {
                std::string name = tokens[1];
                Sequence edgeseq = Sequence(tokens[2]);
                ag::EdgeSaveLabel eids = {{}, {}};
                try {
                    eids = Parse<ag::EdgeSaveLabel>(name, 0, name.size());
                } catch (std::invalid_argument &e) {
                    eids = {{}, {}};
                    ++bad_ids;
                }
                Edge::id_type eid = eids.fId;
                Edge::id_type rceid = eids.rcId;
                edges.emplace_back(edgeseq, name, eid, rceid);
                vertices.add({name, true, true});
                vertices.add({name, true, false});
                vertices.add({name, false, true});
                vertices.add({name, false, false});
            } else if(tokens[0] == "L") {
                size_t overlap = std::stoull(tokens[5].substr(0, tokens[5].size() - 1));
                GFAVertexRecord vrec1(tokens[1], tokens[2] == "-", tokens[2] == "-");
                GFAVertexRecord vrec2(tokens[3], tokens[4] != "-", tokens[4] == "-");
                vertices.link(vrec1, vrec2);
                vertices.link(vrec1.RC(), vrec2.RC());
                vertexLengths[vrec1] = overlap;
                vertexLengths[vrec2] = overlap;
                vertexLengths[vrec1.RC()] = overlap;
                vertexLengths[vrec2.RC()] = overlap;
            }
        }
        is.close();
        MultiGraph res;
        for(std::tuple<Sequence, std::string, Edge::id_type, Edge::id_type> edge : edges) {
            Sequence eseq = std::get<0>(edge);
            Edge::id_type eid = bad_ids > 0 ? Edge::id_type() : std::get<2>(edge);
            Edge::id_type rceid = bad_ids > 0 ? Edge::id_type() : std::get<3>(edge);
            GFAVertexRecord srec(std::get<1>(edge), true, false);
            srec = vertices.get(srec);
            size_t slen = vertexLengths.find(srec) == vertexLengths.end() ? K : vertexLengths[srec];
            if(vertexMap.find(srec) == vertexMap.end())
                vertexMap[srec] = res.addVertex(eseq.Subseq(0, slen)).getId();
            VertexId startId = vertexMap[srec];
            GFAVertexRecord erec(std::get<1>(edge), false, false);
            erec = vertices.get(erec);
            size_t elen = vertexLengths.find(erec) == vertexLengths.end() ? K : vertexLengths[erec];
            if(vertexMap.find(erec) == vertexMap.end())
                vertexMap[erec] = res.addVertex(eseq.Subseq(eseq.size() - elen, eseq.size())).getId();
            VertexId endId = vertexMap[erec];
            startId->addEdge(endId->rc(), eseq, {eid.str()}, eid, rceid);
        }
        return std::move(res);
    }

}
