#include "path_spy.hpp"

void spg::OldVertexTracker::addEmbedding(VertexId old, Segment<Vertex> embedding, std::vector<ag::VertexId> &subs) {
    subs.emplace_back(old);
    vertex_embedding[old] = embedding;
}

void spg::OldVertexTracker::addEmbedding(VertexId old, Segment<Vertex> embedding) {
    addEmbedding(old, embedding, subvertices[embedding.contig().getId()]);
}

spg::OldVertexTracker::OldVertexTracker(ag::ResolutionFire &fire): ag::ResolutionListener(fire, "OldVertexTracker") {
    ag::AssemblyGraph &ag = getFire<ag::AssemblyGraph>();
    for (Vertex &v : ag.vertices()) {
        fireAddVertex(v);
    }
}

inline Segment<ag::Vertex> spg::OldVertexTracker::getPosition(VertexId vid) const {
    auto it = vertex_embedding.find(vid);
    if (it != vertex_embedding.end())
        return it->second;
    return {};
}

void spg::OldVertexTracker::fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) {
    size_t start = 0;
    size_t end = path.getStart().size();
    std::vector<ag::VertexId> &subs = subvertices.at(new_vertex.getId());
    size_t cnt = 0;
    for (Edge &e : path.edges()) {
        end += e.truncSize();
        start += e.rc().truncSize();
        // addEmbedding(e.getFinish().getId(), {new_vertex, start, end}, subs);
        for (VertexId &vid : subvertices.at(e.getFinish().getId())) {
            auto it = subvertices.find(vid);
            if (it != subvertices.end()) {
                const Segment<ag::Vertex> &seg = vertex_embedding.at(vid);
                VERIFY(seg.contig() == e.getFinish());
                addEmbedding(vid, seg.nest({new_vertex, start, end}), subs);
            }
        }
        cnt++;
        if (cnt == path.size() - 1)
            break;
    }
}

void spg::OldVertexTracker::fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {
    // TODO: replace with correct cyclic coordinates of subsegments
    std::vector<ag::VertexId> &subs = subvertices.at(new_vertex.getId());
    for (Edge &e : path.edges()) {
        // addEmbedding(e.getFinish().getId(), {new_vertex, 0, new_vertex.size()}, subs);
        for (VertexId &vid : subvertices.at(e.getFinish().getId())) {
            auto it = subvertices.find(vid);
            if (it != subvertices.end()) {
                addEmbedding(vid, {new_vertex, 0, new_vertex.size()}, subs);
            }
        }
    }
}

void spg::OldVertexTracker::fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) {
    fireMergePath(path, new_edge.isSuffix() ? new_edge.getStart() : new_edge.getFinish());
}

void spg::OldPathTracker::addPath(const std::string &name, const std::vector<ag::AlignmentChain<Contig, ag::Edge>> &als) {
    std::experimental::filesystem::path p = dir / (itos(paths.size(), 2) + "_" + name);
    ensure_dir_existance(p);
    paths.emplace_back(name, als, p);
    for (VertexId &vid : paths.back().vertices) {
        vertex_watch[vid].emplace_back(paths.size() - 1);
    }
}

void spg::OldPathTracker::fireAddSupreVertex(Vertex &v, Edge &e) {
    auto it = vertex_watch.find(e.getStart().getId());
    if (it != vertex_watch.end())
        for (size_t pid : it->second) {
            for (size_t i = 0; i < paths[pid].vertices.size(); i++) {
                if (paths[pid].vertices[i] == e.getStart().getId() && e.getCode() == paths[pid].edge_codes[i]) {
                    paths[pid].vertices[i] = v.getId();
                    paths[pid].edge_codes[i] = {};
                }
            }
        }
}

void spg::OldPathTracker::fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) {
    std::vector<size_t> pathIds;
    for (Vertex &v : path.vertices()) {
        auto it = vertex_watch.find(v.getId());
        if (it != vertex_watch.end()) {
            pathIds.insert(pathIds.end(), it->second.begin(), it->second.end());
            if (it->first != path.getStart().getId() && it->first != path.getFinish().getId()) {
                vertex_watch.erase(it);
            }
        }
    }
    std::sort(pathIds.begin(), pathIds.end());
    pathIds.erase(std::unique(pathIds.begin(), pathIds.end()), pathIds.end());
    vertex_watch[new_vertex.getId()] = pathIds;
    if (pathIds.size() == 0) return;
    ag::Printer p = *printer + ag::EdgeInfo::Colorer(ag::ConstMapping<Edge>(path.edges().begin(), path.edges().end(), "red")) +
                    ag::VertexInfo::Colorer(ag::ConstMapping(new_vertex, "blue"));
    for (size_t pid : pathIds) {
        printPath(pid, p, "MergePath_" + itos(new_vertex.getInnerId()), CollectIds<Vertex>(path.vertices().begin(), path.vertices().end()));
    }
}

void spg::OldPathTracker::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
    auto it = vertex_watch.find(core.getId());
    if (it == vertex_watch.end()) return;
    std::vector<size_t> pathIds = it->second;
    if (pathIds.size() == 0) return;
    ag::Printer p = *printer + ag::VertexInfo::Colorer(ag::ConstMapping(core, "red")) +
        ag::VertexInfo::Colorer(ag::ConstMapping<Vertex>(resolution.newVertices().begin(), resolution.newVertices().end(), "blue"));
    for (size_t pid : pathIds) {
        printPath(pid, p, "ResolveVertex_" + itos(core.getInnerId()), CollectIds<Vertex>(resolution.newVertices().begin(), resolution.newVertices().end()));
    }
    vertex_watch.erase(it);
}

spg::OldPathTracker::WatchPath::WatchPath(std::string name, const std::vector<ag::AlignmentChain<Contig, Edge>> &chain,
    std::experimental::filesystem::path out_dir) :
    name(std::move(name)), out_dir(out_dir) {
    ag::EdgeId prev_eid = {};
    for (const ag::AlignmentChain<Contig, Edge> &al : chain) {
        if (prev_eid != al.seg_to.contig().getId()) {
            vertices.emplace_back(al.seg_to.contig().getStart().getId());
            edge_codes.emplace_back(al.seg_to.contig().getCode());
            prev_eid = al.seg_to.contig().getId();
        }
    }
}

void spg::OldPathTracker::printPath(size_t id, ag::Printer &printer, const std::string &message, const std::vector<ag::VertexId> &extra_vertices) {
    std::experimental::filesystem::path f = paths[id].out_dir / (itos(paths[id].cnt, 4) + + "_" + message + ".dot");
    paths[id].cnt++;
    std::vector<VertexId> vertices = extra_vertices;
    for (VertexId vid : paths[id].vertices) {
        Segment<Vertex> seg = tracker->getPosition(vid);
        if (seg.valid()) {
            vertices.emplace_back(seg.contig().getId());
            VERIFY_MSG(tracker->checkExists(seg.contig().getId()), seg.contig().getId());
        }
    }
    ag::Printer p = printer + ag::VertexInfo::Colorer(ag::ConstMapping<Vertex>(vertices.begin(), vertices.end(), "orange"));
    p.printDot(f, ag::Component::neighbourhood(getFire<ag::AssemblyGraph>(), vertices, 10000, 20000));
}
