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
        SPGListener(fire, "OldVertexTracker"), debug(debug) {
    ag::AssemblyGraph &ag = getFire<ag::AssemblyGraph>();
//    Vertices first: fireAddVertex resets the vertex's records, so it must not run after its edges were recorded.
    for (Vertex &v : ag.vertices()) {
        fireAddVertex(v);
    }
    for (Edge &edge : ag.edges()) {
        fireAddEdge(edge);
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
    direct_supervertices[subvertex.getId()].emplace_back(subvertex, supervertex, left, right);
    direct_subvertices[supervertex.getId()].emplace_back(subvertex, supervertex, left, right);
}

size_t spg::OldVertexTracker::getVertexSize(VertexId vid) const {
    const std::vector<DirectEmbedding> &embeddings = direct_supervertices.at(vid);
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
    direct_supervertices[v.getId()] = {};
    direct_subvertices[v.getId()] = {};
    vertex_seq[v.getId()] = v.getSeq();
}

void spg::OldVertexTracker::fireDeleteVertex(Vertex &v) {
//    Mark v as no longer alive with a self-embedding in its own direct_supervertices. getAllEmbeddings
//    uses this marker (its `real` flag) to skip deleted vertices and climb to the minimal *living*
//    vertices that still contain v's sequence. Without it, getAllEmbeddings can't tell alive from
//    deleted and the "which living vertices is this read a substring of" query is meaningless.
    direct_supervertices[v.getId()].emplace_back(v, v, 0, v.size());
}

void spg::OldVertexTracker::fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) {
//    Every vertex of the path is embedded into new_vertex at its offset; its earlier embeddings stay valid.
    size_t start = 0;
    size_t end = path.getStart().size();
    std::vector<DirectEmbedding> &new_subs = direct_subvertices[new_vertex.getId()];
    direct_supervertices[path.getStart().getId()].emplace_back(path.getStart(), new_vertex, 0, path.getStart().size());
    new_subs.emplace_back(path.getStart(), new_vertex, 0, path.getStart().size());
    for (Edge &e : path.edges()) {
        end += e.truncSize();
        start += e.rc().truncSize();
        VERIFY(e.getFinish().getSeq() == new_vertex.getSeq().Subseq(start, end));
        direct_supervertices[e.getFinish().getId()].emplace_back(e.getFinish(), new_vertex, start, end);
        new_subs.emplace_back(e.getFinish(), new_vertex, start, end);
    }
}

void spg::OldVertexTracker::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
    for (auto p : resolution) {
        addEmbedding(core.getId(), {*p.first, 1, core.size() + 1}, subvertices.at(core.getId()));
//        Also record the embedding in direct_supervertices/direct_subvertices so getAllEmbeddings (which
//        walks those) can climb from the resolved-away core to its living resolution instances. The
//        instance sequence extends core by one nucleotide on each side, so core occupies [1, size+1).
        addMultiEmbedding(core, *p.first, 1, core.size() + 1);
    }
}

bool spg::OldVertexTracker::checkConsistency(std::ostream &os) {
    ag::AssemblyGraph &graph = getGraph();
    bool ok = true;
    std::unordered_set<ag::ConstVertexId> live;
    for (Vertex &v : graph.vertices()) {
        live.insert(v.getId());
//        1.
        auto seq = vertex_seq.find(v.getId());
        if (seq == vertex_seq.end() || seq->second != v.getSeq() || direct_supervertices.count(v.getId()) == 0 ||
                direct_subvertices.count(v.getId()) == 0) {
            os << "OldVertexTracker: live vertex " << v.getId() << " has no records or a wrong sequence" << std::endl;
            ok = false;
        }
    }
    auto contains = [](const std::vector<DirectEmbedding> &list, const DirectEmbedding &de) {
        return std::find(list.begin(), list.end(), de) != list.end();
    };
    for (const auto &p : direct_supervertices) {
        bool dead_marker = false;
        for (const DirectEmbedding &de : p.second) {
            if (de.inner_vertex != p.first) {
                os << "OldVertexTracker: embedding filed under the wrong vertex " << p.first << std::endl;
                ok = false;
            } else if (de.outer_vertex == de.inner_vertex) {
                dead_marker = true;
            } else {
//                2, 3.
                auto inner = vertex_seq.find(de.inner_vertex), outer = vertex_seq.find(de.outer_vertex);
                auto subs = direct_subvertices.find(de.outer_vertex);
                if (inner == vertex_seq.end() || outer == vertex_seq.end() || de.to > outer->second.size() ||
                        de.from > de.to || inner->second != outer->second.Subseq(de.from, de.to)) {
                    os << "OldVertexTracker: embedding of " << de.inner_vertex << " into " << de.outer_vertex
                       << " [" << de.from << "," << de.to << ") does not hold" << std::endl;
                    ok = false;
                }
                if (subs == direct_subvertices.end() || !contains(subs->second, de)) {
                    os << "OldVertexTracker: embedding of " << de.inner_vertex << " into " << de.outer_vertex
                       << " is missing among the subvertices of " << de.outer_vertex << std::endl;
                    ok = false;
                }
            }
        }
//        4.
        if (dead_marker == (live.count(p.first) > 0)) {
            os << "OldVertexTracker: vertex " << p.first << (dead_marker ? " is live but marked dead" : " is dead but not marked") << std::endl;
            ok = false;
        }
    }
    for (const auto &p : direct_subvertices)
        for (const DirectEmbedding &de : p.second) {
            auto sups = direct_supervertices.find(de.inner_vertex);
            if (de.outer_vertex != p.first || sups == direct_supervertices.end() || !contains(sups->second, de)) {
                os << "OldVertexTracker: subvertex record of " << de.inner_vertex << " in " << p.first
                   << " has no matching supervertex record" << std::endl;
                ok = false;
            }
        }
//    5.
    auto recorded = [&](Vertex &inner, const DirectEmbedding &de) {
        auto it = direct_supervertices.find(inner.getId());
        return it != direct_supervertices.end() && contains(it->second, de);
    };
    for (Edge &e : graph.edges()) {
        bool missing = false;
        if (e.isPrefix())
            missing = !recorded(e.getStart(), DirectEmbedding(e.getStart(), e.getFinish(), 0, e.getStart().size()));
        else if (e.isSuffix())
            missing = !recorded(e.getFinish(), DirectEmbedding(e.getFinish(), e.getStart(), e.rc().truncSize(), e.getStart().size()));
        if (missing) {
            os << "OldVertexTracker: the embedding of edge " << e.getId() << " is not recorded" << std::endl;
            ok = false;
        }
    }
    return ok;
}

std::vector<spg::OldVertexTracker::DirectEmbedding> spg::OldVertexTracker::getAllEmbeddings(VertexId vid) const {
    std::vector<DirectEmbedding> res;
    std::set<std::pair<size_t, DirectEmbedding>> recs;
    size_t sz = getVertexSize(vid);
//    Use the VertexId ctor (not *vid): vid may be a deleted/legacy vertex, which is exactly the query
//    this method answers, so it must never be dereferenced.
    recs.emplace(sz, DirectEmbedding(vid, vid, 0, sz));
    std::vector<DirectEmbedding> result;
    while (!recs.empty()) {
        auto it = recs.begin();
        DirectEmbedding de = it->second;
        VERIFY(de.inner_vertex == vid);
        recs.erase(it);
        bool real = true;
        for (const DirectEmbedding &next : direct_supervertices.at(de.outer_vertex)) {
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
        for (const DirectEmbedding &next : direct_subvertices.at(de.inner_vertex)) {
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

std::vector<ag::VertexId> spg::readSubstringVertices(const ag::AlignedRead &read, const OldVertexTracker &tracker) {
    std::set<ag::VertexId> seen;
    if (!read.getPath().valid())
        return {};
    ag::VertexId v = read.getPath().getStartId();
    for (const OldVertexTracker::DirectEmbedding &de : tracker.getAllEmbeddings(v))
        seen.insert(de.outer_vertex);
//    TODO: getAllEmbeddings returns every living vertex that contains v at any level; filtering to the
//    minimal-by-substring containers (and canonicalising orientation) is a further refinement.
    return {seen.begin(), seen.end()};
}
