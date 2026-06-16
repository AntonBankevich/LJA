#include "path_spy.hpp"

#include <map>

bool spg::OldVertexTracker::DirectEmbedding::operator<(const DirectEmbedding &other) const {
    return std::tie(inner_vertex, outer_vertex, from, to) <
           std::tie(other.inner_vertex, other.outer_vertex, other.from, other.to);
}

bool spg::OldVertexTracker::DirectEmbedding::operator==(const DirectEmbedding &other) const {
    return inner_vertex == other.inner_vertex && outer_vertex == other.outer_vertex && from == other.from && to == other.to;
}

void spg::OldVertexTracker::addEmbedding(VertexId old, Segment<Vertex> embedding, std::vector<ag::VertexId> &subs) {
    subs.emplace_back(old);
    vertex_embedding[old] = embedding;
}

void spg::OldVertexTracker::addEmbedding(VertexId old, Segment<Vertex> embedding) {
    addEmbedding(old, embedding, subvertices[embedding.contig().getId()]);
}

spg::OldVertexTracker::OldVertexTracker(ag::ResolutionFire &fire, bool debug):
        ag::ResolutionListener(fire, "OldVertexTracker"), debug(debug) {
    ag::AssemblyGraph &ag = getFire<ag::AssemblyGraph>();
    for (Edge &edge : ag.edges()) {
        fireAddEdge(edge);
    }
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

void spg::OldVertexTracker::addMultiEmbedding(Vertex &subvertex, Vertex &supervertex, size_t left, size_t right) {
    VERIFY(right >= left);
    VERIFY(subvertex.size() == right - left);
    VERIFY(subvertex.size() < supervertex.size() || subvertex == supervertex);
    multi_embedding[subvertex.getId()].emplace_back(subvertex, supervertex, left, right);
    multi_subvertices[supervertex.getId()].emplace_back(subvertex, supervertex, left, right);
}

size_t spg::OldVertexTracker::getVertexSize(VertexId vid) const {
    const std::vector<DirectEmbedding> &embeddings = multi_embedding.at(vid);
    if (embeddings.size() > 0) return embeddings.front().to - embeddings.front().from;
    return vid->size();
}

void spg::OldVertexTracker::fireAddEdge(Edge &e) {
    if (e.isPrefix()) {
        addMultiEmbedding(e.getStart(), e.getFinish(), 0, e.getStart().size());
    } else if (e.isSuffix()) {
        addMultiEmbedding(e.getFinish(), e.getStart(), e.rc().truncSize(), e.getStart().size());
    }
}

void spg::OldVertexTracker::fireAddVertex(Vertex &v) {
    subvertices[v.getId()] = {};
    addEmbedding(v.getId(), {v, 0, v.size()});
    multi_embedding[v.getId()] = {};
    multi_subvertices[v.getId()] = {};
    vertex_seq[v.getId()] = v.getSeq();
}

void spg::OldVertexTracker::fireDeleteVertex(Vertex &v) {
    subvertices.erase(v.getId());
    multi_embedding[v.getId()].emplace_back(v, v, 0, v.size());
}

void spg::OldVertexTracker::fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) {
    size_t start = 0;
    size_t end = path.getStart().size();
    std::vector<ag::VertexId> &subs = subvertices.at(new_vertex.getId());
    std::vector<DirectEmbedding> &multi_subs = multi_subvertices[new_vertex.getId()];
    size_t cnt = 0;
    for (Edge &e : path.edges()) {
        end += e.truncSize();
        start += e.rc().truncSize();
        multi_embedding[e.getFinish().getId()] = {{e.getFinish(), new_vertex, start, end}};
        multi_subs.emplace_back(e.getFinish(), new_vertex, start, end);
        VERIFY(e.getFinish().getSeq() == new_vertex.getSeq().Subseq(start, end));
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

void spg::OldVertexTracker::fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) {
    VERIFY(new_edge.isSuffix() || new_edge.isPrefix());
    fireMergePath(path, new_edge.isSuffix() ? new_edge.getStart() : new_edge.getFinish());
}

void spg::OldVertexTracker::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
    if (resolution.newVertices().calculateSize() != 1) {
        for (ag::VertexId vid : subvertices.at(core.getId())) {
            vertex_embedding.erase(vid);
        }
    } else {
        auto new_connection = *resolution.begin();
        Vertex& new_vertex = *new_connection.first;
        Segment<Vertex> seg = {new_vertex, new_connection.second.incoming().rc().truncSize(),
            new_connection.second.incoming().rc().truncSize() + core.size()};
        std::vector<VertexId> & new_subvertices = subvertices[new_vertex.getId()];
        for (ag::VertexId vid : subvertices.at(core.getId())) {
            Segment<Vertex> &embedding = vertex_embedding[vid];
            addEmbedding(vid, embedding.nest(seg), new_subvertices);
        }
    }
    subvertices.erase(core.getId());
}

std::vector<spg::OldVertexTracker::DirectEmbedding> spg::OldVertexTracker::getAllEmbeddings(VertexId vid) const {
    std::vector<DirectEmbedding> res;
    std::set<std::pair<size_t, DirectEmbedding>> recs;
    size_t sz = getVertexSize(vid);
    recs.emplace(sz, DirectEmbedding(*vid, *vid, 0, sz));
    std::vector<DirectEmbedding> result;
    while (!recs.empty()) {
        auto it = recs.begin();
        DirectEmbedding de = it->second;
        VERIFY(de.inner_vertex == vid);
        recs.erase(it);
        bool real = true;
        for (const DirectEmbedding &next : multi_embedding.at(de.outer_vertex)) {
            if (next.inner_vertex == next.outer_vertex) {
                real = false;
                continue;
            }
            recs.emplace(getVertexSize(next.outer_vertex), DirectEmbedding(vid, next.outer_vertex, de.from + next.from, de.to + next.from));
        }
        if (real) {
            res.emplace_back(de);
        }
    }
    return res;
}

std::vector<spg::OldVertexTracker::DirectEmbedding> spg::OldVertexTracker::getAllSubvertices(VertexId vid) const {
    std::vector<DirectEmbedding> result;
    std::set<std::pair<size_t, DirectEmbedding>> recs;
    size_t sz = getVertexSize(vid);
    recs.emplace(sz, DirectEmbedding(*vid, *vid, 0, sz));
    while (!recs.empty()) {
        auto it = recs.end();
        --it;
        DirectEmbedding de = it->second;
        VERIFY(de.outer_vertex == vid);
        recs.erase(it);
        for (const DirectEmbedding &next : multi_subvertices.at(de.inner_vertex)) {
            if (next.inner_vertex == next.outer_vertex) {
                continue;
            }
            recs.emplace(getVertexSize(next.inner_vertex), DirectEmbedding(next.inner_vertex, vid, de.from + next.from, de.from + next.to));
        }
        result.emplace_back(de);
    }
    for (DirectEmbedding &de : result) {
        VERIFY(vertex_seq.at(de.inner_vertex) == vertex_seq.at(de.outer_vertex).Subseq(de.from, de.to));
    }
    return result;
}

void spg::OldPathTracker::addPath(const std::string &name, const std::vector<ag::AlignmentChain<Contig, ag::Edge>> &als) {
    std::experimental::filesystem::path p = dir / (itos(paths.size(), 2) + "_" + name);
    ensure_dir_existance(p);
    paths.emplace_back(name, als, p);
    for (VertexId &vid : paths.back().vertices) {
        vertex_watch[vid].emplace_back(paths.size() - 1);
    }
}

void spg::OldPathTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
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
