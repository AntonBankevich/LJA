#include "path_spy.hpp"

#include <map>

bool spg::OldVertexTracker::DirectEmbedding::operator<(const DirectEmbedding &other) const {
    return std::tie(inner_vertex, outer_vertex, from, to) <
           std::tie(other.inner_vertex, other.outer_vertex, other.from, other.to);
}

bool spg::OldVertexTracker::DirectEmbedding::operator==(const DirectEmbedding &other) const {
    return inner_vertex == other.inner_vertex && outer_vertex == other.outer_vertex && from == other.from && to == other.to;
}

void spg::OldVertexTracker::append(EmbeddingMap &map, ag::ConstVertexId key, const DirectEmbedding &embedding) {
    map.upsert(key, [&embedding](std::unique_ptr<EmbeddingRecord> &rec) {rec->emplace_back(embedding);},
               std::make_unique<EmbeddingRecord>(1, embedding));
}

const spg::OldVertexTracker::EmbeddingRecord &spg::OldVertexTracker::getRecord(const EmbeddingMap &map, ag::ConstVertexId key) {
    const EmbeddingRecord *res = nullptr;
    map.find_fn(key, [&res](const std::unique_ptr<EmbeddingRecord> &ptr) {res = ptr.get();});
    VERIFY(res != nullptr);
    return *res;
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

void spg::OldVertexTracker::addMultiEmbedding(Vertex &subvertex, Vertex &supervertex, size_t left, size_t right) {
    VERIFY(right >= left);
    VERIFY(subvertex.size() == right - left);
//    Equal sizes: a vertex into itself (the dead marker) or a resolved core into its copy.
    VERIFY(subvertex.size() <= supervertex.size());
    DirectEmbedding embedding(subvertex, supervertex, left, right);
    append(direct_supervertices, subvertex.getId(), embedding);
    append(direct_subvertices, supervertex.getId(), embedding);
}

size_t spg::OldVertexTracker::getVertexSize(VertexId vid) const {
    const EmbeddingRecord &embeddings = getRecord(direct_supervertices, vid);
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
    ag::ConstVertexId vid = v.getId();
    direct_supervertices.insert_or_assign(vid, std::make_unique<EmbeddingRecord>());
    direct_subvertices.insert_or_assign(vid, std::make_unique<EmbeddingRecord>());
    vertex_seq.insert_or_assign(vid, debug ? v.getSeq() : Sequence());
}

void spg::OldVertexTracker::fireDeleteVertex(Vertex &v) {
//    Mark v as no longer alive with a self-embedding in its own direct_supervertices. getAllEmbeddings
//    uses this marker (its `real` flag) to skip deleted vertices and climb to the minimal *living*
//    vertices that still contain v's sequence. Without it, getAllEmbeddings can't tell alive from
//    deleted and the "which living vertices is this read a substring of" query is meaningless.
    append(direct_supervertices, v.getId(), DirectEmbedding(v, v, 0, v.size()));
}

//    Every vertex of the old path that the dominant vertex contains is embedded into it at its offset; earlier
//    embeddings stay valid.
void spg::OldVertexTracker::fireReplacePath(const ag::PathReplacement &rep) {
    const ag::RAGraphPath &path = rep.oldPath();
    Vertex &dominant = rep.getDominant();
    for (size_t i = 0; i <= path.size(); i++) {
        Vertex &v = path.getVertex(i);
        if (v == dominant || rep.oldPosition(i) < rep.dominantPosition())
            continue;
        size_t start = rep.oldPosition(i) - rep.dominantPosition();
        size_t end = start + v.size();
        if (end > dominant.size())
            continue;
        VERIFY(v.getSeq() == dominant.getSeq().Subseq(start, end));
        DirectEmbedding embedding(v, dominant, start, end);
        append(direct_supervertices, v.getId(), embedding);
        append(direct_subvertices, dominant.getId(), embedding);
    }
}

void spg::OldVertexTracker::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
    for (auto p : resolution) {
//        Record the embedding so getAllEmbeddings can climb from the resolved-away core to its living resolution
//        instances. An instance extends core by one nucleotide on each side, so core occupies [1, size+1), or is a copy
//        of core (1-to-1 pair).
        const size_t d = p.first->size() == core.size() ? 0 : 1;
        addMultiEmbedding(core, *p.first, d, core.size() + d);
    }
}

bool spg::OldVertexTracker::checkConsistency(std::ostream &os) {
    ag::AssemblyGraph &graph = getGraph();
    bool ok = true;
    auto seqs = vertex_seq.lock_table();
    auto sups = direct_supervertices.lock_table();
    auto subs = direct_subvertices.lock_table();
    std::unordered_set<ag::ConstVertexId> live;
    for (Vertex &v : graph.vertices()) {
        live.insert(v.getId());
//        1.
        auto seq = seqs.find(ag::ConstVertexId(v.getId()));
        if (seq == seqs.end() || (debug && seq->second != v.getSeq()) || sups.find(ag::ConstVertexId(v.getId())) == sups.end() ||
                subs.find(ag::ConstVertexId(v.getId())) == subs.end()) {
            os << "OldVertexTracker: live vertex " << v.getId() << " has no records or a wrong sequence" << std::endl;
            ok = false;
        }
    }
    auto contains = [](const std::vector<DirectEmbedding> &list, const DirectEmbedding &de) {
        return std::find(list.begin(), list.end(), de) != list.end();
    };
    for (const auto &p : sups) {
        bool dead_marker = false;
        for (const DirectEmbedding &de : *p.second) {
            if (de.inner_vertex != p.first) {
                os << "OldVertexTracker: embedding filed under the wrong vertex " << p.first << std::endl;
                ok = false;
            } else if (de.outer_vertex == de.inner_vertex) {
                dead_marker = true;
            } else {
//                2, 3.
                auto inner = seqs.find(ag::ConstVertexId(de.inner_vertex));
                auto outer = seqs.find(ag::ConstVertexId(de.outer_vertex));
                auto outer_subs = subs.find(ag::ConstVertexId(de.outer_vertex));
                if (inner == seqs.end() || outer == seqs.end() || de.from > de.to ||
                        (debug && (de.to > outer->second.size() || inner->second != outer->second.Subseq(de.from, de.to)))) {
                    os << "OldVertexTracker: embedding of " << de.inner_vertex << " into " << de.outer_vertex
                       << " [" << de.from << "," << de.to << ") does not hold" << std::endl;
                    ok = false;
                }
                if (outer_subs == subs.end() || !contains(*outer_subs->second, de)) {
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
    for (const auto &p : subs)
        for (const DirectEmbedding &de : *p.second) {
            auto inner_sups = sups.find(ag::ConstVertexId(de.inner_vertex));
            if (de.outer_vertex != p.first || inner_sups == sups.end() || !contains(*inner_sups->second, de)) {
                os << "OldVertexTracker: subvertex record of " << de.inner_vertex << " in " << p.first
                   << " has no matching supervertex record" << std::endl;
                ok = false;
            }
        }
//    5.
    auto recorded = [&](Vertex &inner, const DirectEmbedding &de) {
        auto it = sups.find(ag::ConstVertexId(inner.getId()));
        return it != sups.end() && contains(*it->second, de);
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
        for (const DirectEmbedding &next : getRecord(direct_supervertices, de.outer_vertex)) {
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
        for (const DirectEmbedding &next : getRecord(direct_subvertices, de.inner_vertex)) {
            if (next.inner_vertex == next.outer_vertex) {
                continue;
            }
            recs.emplace(getVertexSize(next.inner_vertex), DirectEmbedding(next.inner_vertex, vid, de.from + next.from, de.from + next.to));
        }
        result.emplace_back(de);
    }
    if (debug) {
        for (DirectEmbedding &de : result) {
            VERIFY(vertex_seq.find(ag::ConstVertexId(de.inner_vertex)) ==
                   vertex_seq.find(ag::ConstVertexId(de.outer_vertex)).Subseq(de.from, de.to));
        }
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

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
spg::PrintStorageAlignments(logging::Logger &logger, size_t threads, const ag::AlignedReadStorage &storage,
                            const OldVertexTracker &tracker, const std::experimental::filesystem::path &dir) {
    logger.info() << "Collecting read alignments from read paths" << std::endl;
    std::vector<ag::VertexId> single_vertices;
    std::unordered_map<ag::ConstVertexId, size_t> single_index;
    size_t unknown_reads = 0;
    for (const ag::AlignedRead &read : storage) {
        const ag::GraphPath &path = read.getPath();
        if (!path.valid() || !path.empty())
            continue;
        ag::VertexId vid = path.getStartId();
        if (!tracker.checkExists(vid)) {
            unknown_reads++;
            continue;
        }
        if (single_index.emplace(vid, single_vertices.size()).second)
            single_vertices.emplace_back(vid);
    }
    std::vector<std::vector<OldVertexTracker::DirectEmbedding>> containers(single_vertices.size());
    omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic, 1) default(none) shared(single_vertices, containers, tracker)
    for (size_t i = 0; i < single_vertices.size(); i++) {
        containers[i] = tracker.getAllEmbeddings(single_vertices[i]);
    }

    ensure_dir_existance(dir);
    std::experimental::filesystem::path good_fname = dir / "good_alignments.txt";
    std::experimental::filesystem::path bad_fname = dir / "partial_alignments.txt";
    logger.info() << "Printing read alignments to " << good_fname << std::endl;
    std::ofstream os(good_fname);
    std::ofstream os_bad(bad_fname);
    size_t good_cnt = 0;
    size_t bad_cnt = 0;
    auto print = [&](const std::string &read_id, size_t read_len, Vertex::id_type contig_id, size_t contig_len,
                     size_t read_from, size_t read_to, size_t contig_from, size_t contig_to) {
        if (read_from >= read_to)
            return;
        VERIFY(read_to - read_from == contig_to - contig_from);
        VERIFY(read_to <= read_len && contig_to <= contig_len);
        bool good = (read_from == 0 || contig_from == 0) && (read_to == read_len || contig_to == contig_len);
        (good ? good_cnt : bad_cnt)++;
        (good ? os : os_bad) << read_id << " " << read_from << " " << read_to << " " << contig_id << " "
                             << contig_from << " " << contig_to << "\n";
    };
    for (const ag::AlignedRead &read : storage) {
        const ag::GraphPath &path = read.getPath();
        if (!path.valid())
            continue;
        if (path.empty()) {
            auto it = single_index.find(path.getStartId());
            if (it == single_index.end())
                continue;
            for (const OldVertexTracker::DirectEmbedding &de : containers[it->second]) {
                size_t vertex_size = de.to - de.from;
                VERIFY(path.leftCut() + path.rightCut() <= vertex_size);
                size_t read_len = vertex_size - path.leftCut() - path.rightCut();
                print(read.getId(), read_len, de.outer_vertex.innerId(), de.outer_vertex->size(),
                      0, read_len, de.from + path.leftCut(), de.to - path.rightCut());
            }
            continue;
        }
//        The read occupies [read_from, read_to) of its path's label; a vertex of the path ends where the label of the
//        path up to it ends.
        size_t read_len = path.len();
        size_t read_from = path.leftCut();
        size_t read_to = read_from + read_len;
        auto print_vertex = [&](const Vertex &v, size_t end) {
            VERIFY(end >= v.size());
            size_t start = end - v.size();
            size_t left = std::max(start, read_from);
            size_t right = std::min(end, read_to);
            if (left < right)
                print(read.getId(), read_len, v.getInnerId(), v.size(),
                      left - read_from, right - read_from, left - start, right - start);
        };
        size_t end = path.getStart().size();
        print_vertex(path.getStart(), end);
        for (Edge &edge : path.edges()) {
            end += edge.truncSize();
            print_vertex(edge.getFinish(), end);
        }
    }
    os.close();
    os_bad.close();
    logger.info() << "Printed " << good_cnt << " alignments and " << bad_cnt << " partial alignments" << std::endl;
    if (unknown_reads > 0)
        logger.info() << unknown_reads << " reads lie inside vertices unknown to the vertex tracker and got no alignments"
                      << std::endl;
    return {good_fname, bad_fname};
}
