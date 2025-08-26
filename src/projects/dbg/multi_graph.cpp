#include <common/disjoint_sets.hpp>
#include <common/logging.hpp>
#include <assembly_graph/data_structures/component.hpp>
#include <assembly_graph/data_structures/splitters.hpp>
#include <assembly_graph/visualization.hpp>
#include "multi_graph.hpp"

namespace multigraph {

    MultiGraph MultiGraphHelper::TransformToEdgeGraph(logging::Logger &logger, const MultiGraph &mg, size_t tip_size) {
        MultiGraph dbg;
        std::unordered_map<ConstEdgeId, VertexId> emap;
        for (const Vertex &v: mg.vertices()) {
            if (v.outDeg() == 0 || emap.find(v.front().getId()) != emap.end()) {
                continue;
            }
            Vertex &newv = dbg.addVertex(v.getSeq().Subseq(v.size() - v.front().overlapSize()));
            for (const Edge &edge: v) {
                const Vertex &right = edge.getFinish();
                for (const Edge &edge1: right.rc()) {
                    emap[edge1.getId()] = newv.rc().getId();
                    emap[edge1.rc().getId()] = newv.getId();
                }
            }
        }
        for (const Vertex &v: mg.verticesUnique()) {
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
            Edge &edge = dbg.addEdgeLockFree(*start, *end, v.getSeq());
            logger.trace() << "Edge id changed from: " << edge.getId() << " to " << v.getId() << std::endl;
        }
        return std::move(dbg);
    }

    MultiGraph MultiGraphHelper::Delete(const MultiGraph &initial, const std::unordered_set<ConstEdgeId> &to_delete,
                                        const std::unordered_set<ConstVertexId> &to_delete_vertices) {
        MultiGraph res;
        std::unordered_map<ConstVertexId, VertexId> vmap;
        std::unordered_set<ConstEdgeId> visited;
        for (const Vertex &v: initial.verticesUnique()) {
            if (to_delete_vertices.find(v.getId()) != to_delete_vertices.end())
                continue;
            vmap[v.getId()] = res.addVertex(v).getId();
            vmap[v.rc().getId()] = vmap[v.getId()]->rc().getId();
        }
        for (const Edge &edge: initial.edgesUnique()) {
            if (to_delete.find(edge.getId()) != to_delete.end())
                continue;
            if (to_delete_vertices.find(edge.getStart().getId()) == to_delete_vertices.end() ||
                to_delete_vertices.find(edge.getFinish().getId()) == to_delete_vertices.end())
                res.addEdgeLockFree(*vmap[edge.getStart().getId()], *vmap[edge.getFinish().getId()], edge.getSeq());
        }
        return std::move(res);
    }

//    std::vector<EdgeId> MultiGraphHelper::uniquePathForward(Edge &edge) {
//        std::vector<EdgeId> res = {edge.getId()};
//        VertexId cur = edge.getFinish().getId();
//        while (cur != edge.getStart().getId() && cur->inDeg() == 1 && cur->outDeg() == 1) {
//            res.emplace_back(cur->begin()->getId());
//            cur = res.back()->getFinish().getId();
//        }
//        return std::move(res);
//    }
//
//    std::vector<ConstEdgeId> MultiGraphHelper::uniquePathForward(const Edge &edge) {
//        std::vector<ConstEdgeId> res = {edge.getId()};
//        ConstVertexId cur = edge.getFinish().getId();
//        while (cur != edge.getStart().getId() && cur->inDeg() == 1 && cur->outDeg() == 1) {
//            res.emplace_back(cur->begin()->getId());
//            cur = res.back()->getFinish().getId();
//        }
//        return std::move(res);
//    }
//
//    std::vector<ConstEdgeId> MultiGraphHelper::uniquePath(const Edge &edge) {
//        std::vector<ConstEdgeId> path = uniquePathForward(edge.rc());
//        return uniquePathForward(path.back()->rc());
//    }
//
//    std::vector<EdgeId> MultiGraphHelper::uniquePath(Edge &edge) {
//        std::vector<EdgeId> path = uniquePathForward(edge.rc());
//        return uniquePathForward(path.back()->rc());
//    }

//    MultiGraph MultiGraphHelper::MergeAllPaths(const MultiGraph &mg, bool verbose) {
//        MultiGraph res;
//        std::unordered_set<ConstEdgeId> used;
//        std::unordered_map<ConstVertexId, VertexId> old_to_new;
//        for(const Edge &edge: mg.edges()) {
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
//                Edge &new_edge = new_start->addEdge(*new_end, sb.BuildSequence());
//                if(verbose) {
//                    std::cout << "New getEdge " << new_edge.getId() << " consists of old edges: ";
//                    for(auto e : path) {
//                        std::cout << e->getId() << " ";
//                    }
//                    std::cout << std::endl;
//                }
//            }
//        }
//        for(const Vertex &vertex : mg.vertices()) {
//            if(vertex.inDeg() == 0 && vertex.outDeg() == 0 && vertex.isCanonical()) {
//                res.addVertex(vertex.getSeq());
//            }
//        }
//        MultiGraphHelper::checkConsistency(res);
//        return std::move(res);
//    }

    std::vector<Contig> MultiGraphHelper::extractContigs(const MultiGraph &mg, bool cut_overlaps) {
        std::unordered_map<ConstVertexId, size_t> cut;
        for (const Vertex &v: mg.vertices()) {
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
        for (const Edge &edge: mg.edges()) {
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
        std::unordered_map<const Vertex *, int> vmap;
        for (const Vertex &vertex: mg.vertices()) {
            os << vertex.getId() << " [label=\"" << vertex.size() << "\" style=filled fillcolor=\"white\"]\n";
        }
        std::unordered_map<EdgeId, std::string> eids;
        for (const Edge &edge: mg.edges()) {
            os << "\"" << edge.getStart().getId() << "\" -> \"" << edge.getFinish().getId() <<
               "\" [label=\"" << edge.getId() << "(" << edge.fullSize() << ")\" color = \"black\"]\n";

        }
        os << "}\n";
        os.close();
    }

    void MultiGraphHelper::printDot2(const MultiGraph &mg, const std::experimental::filesystem::path &f) {
        std::ofstream os;
        os.open(f);
        os << "digraph {\nnodesep = 0.5;\n";
        std::unordered_map<const Vertex *, int> vmap;
        for (const Vertex &vertex: mg.vertices()) {
            os << vertex.getId() << " [label=\"" << vertex.getId() << " : " << vertex.size() << "\" style=filled fillcolor=\"white\"]\n";
        }
        std::unordered_map<EdgeId, std::string> eids;
        for (const Edge &edge: mg.edges()) {
            os << "\"" << edge.getStart().getId() << "\" -> \"" << edge.getFinish().getId() <<
               "\" [label=\"" << edge.getId() << "(" << edge.fullSize() << ")\" color = \"black\"]\n";

        }
        os << "}\n";
        os.close();
    }

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
		Vertex &newV = int_ids ? res.addVertex(seq, seq.isCanonical() ? std::stoi(name) : -std::stoi(name)) : res.addVertex(seq);

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
                res.addEdge(*v1, *v2, v1->getSeq() + v2->getSeq().Subseq(overlap));
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

    LabelStorage::LabelStorage(ag::AssemblyGraph &fire) : ag::ResolutionListener(fire, "LabelStorage") {
        for(Edge &e : fire.edges())
            labels[e.getId()] = {e.getId()};
    }
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
            res.addEdge(*startId, endId->rc(), eseq, eid, rceid);
        }
        return std::move(res);
    }

}
