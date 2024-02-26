#include "mult_correction.hpp"
using namespace dbg;
void printAl(logging::Logger &logger, std::unordered_map<const dbg::Edge *, CompactPath> &unique_extensions,
             const dbg::GraphPath &al) {
    for(Edge &piece : al.edges()) {
        logger << piece.str() << " ";
        if(unique_extensions.find(&piece) != unique_extensions.end()) {
            logger << "+ ";
        }
    }
    logger << std::endl;
}

struct UEdge {
    UEdge(dbg::Edge *from, dbg::Edge *to, dbg::CompactPath cpath, size_t support) : from(from), to(to), cpath(std::move(cpath)),
                                                                                    support(support) {}

    dbg::Edge *from;
    dbg::Edge *to;
    dbg::CompactPath cpath;
    size_t support;

    UEdge RC() const {
        return {&to->rc(), &from->rc(), cpath.RC(), support};
    }
};
//std::unordered_map<const Edge *, CompactPath> constructUniqueExtensions(logging::Logger &logger, SparseDBG &dbg,
//                                                                         const RecordStorage &reads_storage, const UniqueClassificator &classificator) {
//    std::unordered_map<Edge *, std::vector<UEdge>> bg;
//    for(Edge &edge : dbg.edges()) {
//        if(!classificator.isUnique(edge))
//            continue;
//        Vertex &start = *edge.start();
//        std::vector<Sequence> extensions;
//        for(auto & c : reads_storage.getRecord(start)) {
//            GraphAlignment al = CompactPath(start, c.first).getAlignment();
//            if(al.front().contig() != edge)
//                continue;
//            for(size_t i = 0; i < al.size(); i++) {
//                Segment<Edge> &seg = al[i];
//                if(classificator.isUnique(seg.contig())) {
//                    al = al.subalignment(0, i + 1);
//                    break;
//                }
//            }
//            if(!classificator.isUnique(al.back().contig()))
//                continue;
//            extensions.emplace_back(CompactPath(al).cpath());
//        }
//        std::sort(extensions.begin(), extensions.end());
//        extensions.erase(std::unique(extensions.begin(), extensions.end()), extensions.end());
//        for(Sequence &extension : extensions) {
//            CompactPath new_path(start, extension);
//            bg[&edge].emplace_back(&edge, &new_path.getAlignment().back().contig(), new_path,
//                                   reads_storage.getRecord(start).countStartsWith(extension));
//        }
//    }
//    std::unordered_map<Edge *, UEdge> choice;
//}

inline void findEasyExtensions(const std::vector<Edge *> &uniqueEdges, const RecordStorage &reads_storage,
                               const AbstractUniquenessStorage &classificator,
                               std::unordered_map<const Edge *, CompactPath> &unique_extensions) {
    VERIFY(uniqueEdges.empty() || uniqueEdges.front()->truncSize() >= uniqueEdges.back()->truncSize());
    for(Edge *edgeIt : uniqueEdges) {
        Edge &edge = *edgeIt;
        if (unique_extensions.find(&edge) != unique_extensions.end())
            continue;
        Vertex & start = edge.getStart();
        const VertexRecord &rec = reads_storage.getRecord(start);
        Sequence seq = edge.truncSeq().Subseq(0, 1);
        CompactPath path = rec.getFullUniqueExtension(seq, 1, 0);
        if(path.size() == 1)
            continue;
        dbg::GraphPath al = path.unpack();
        for(size_t i = 1; i < al.size(); i++) {
            Segment<Edge> seg = al[i];
            if(classificator.isUnique(seg.contig())) {
                al = al.subPath(0, i + 1);
                break;
            }
        }
        if(!classificator.isUnique(al.back().contig()) || al.size() == 1)
            continue;
        unique_extensions.emplace(&edge, CompactPath(edge.getFinish(), CompactPath(al).cpath().Subseq(1), 0, 0));
        CompactPath res1(al.RC().subPath(1, al.size()));
        unique_extensions.emplace(&al.back().contig().rc(), res1);
    }
}

dbg::GraphPath greedyExtension(const VertexRecord &rec, const AbstractUniquenessStorage &classificator, Edge &edge) {
    dbg::GraphPath path(edge.getStart());
    path += edge;
    Sequence seq = CompactPath(path).cpath();
    while(true) {
        Sequence best;
        size_t best_val = 0;
        Edge *next_edge = nullptr;
        for(Edge &next_candidate : path.finish()) {
            if(classificator.isError(next_candidate))
                continue;
            Sequence next = seq + next_candidate.truncSeq().Subseq(0, 1);
            size_t val = rec.countStartsWith(next);
            if(val > best_val) {
                best_val = val;
                best = next;
                next_edge = &next_candidate;
            }
        }
        if(best_val == 0)
            break;
        seq = best;
        path += *next_edge;
    }
    return path;
}

inline CompactPath findBulgeExtension(const VertexRecord &rec, Edge &edge, const CompactPath & greedy) {
    if(edge.getFinish().outDeg() != 2)
        return greedy;
    Edge &edge1 = edge.getFinish().front();
    Edge &edge2 = edge.getFinish().back();
    CompactPath cp1 = rec.getFullUniqueExtension(edge.firstNucl() + edge1.firstNucl(), 1, 0);
    CompactPath cp2 = rec.getFullUniqueExtension(edge.firstNucl() + edge2.firstNucl(), 1, 0);
    if(greedy.cpath().startsWith(cp2.cpath())) {
        std::swap(cp1, cp2);
    } else {
        if(!greedy.cpath().startsWith(cp1.cpath()))
            return greedy;
    }
    dbg::GraphPath p1 = cp1.unpack();
    dbg::GraphPath p2 = cp2.unpack();
    size_t b1 = 0;
    size_t b2 = 0;
    Sequence choice;
    if(p2.find(p1.getVertex(2), 2) != size_t (-1)) {
        b1 = 2;
        b2 = p2.find(p1.getVertex(2), 2);
        choice = cp2.cpath().Subseq(0, b2);
        if(p1.find(p2.getVertex(2), 2) != size_t(-1)) {
            return greedy;
        }
    } else if(p1.find(p2.getVertex(2), 2) != size_t(-1)) {
        b1 = p1.find(p2.getVertex(2), 2);
        b2 = 2;
        choice = cp1.cpath().Subseq(0, b1);
    }
    if(!cp1.cpath().Subseq(b1).nonContradicts(cp2.cpath().Subseq(b2)))
        return greedy;
    return CompactPath(edge.getStart(), choice + greedy.cpath().Subseq(b1));
}

inline void findComplexExtensions(const std::vector<Edge *> &uniqueEdges, const RecordStorage &reads_storage,
                                  const AbstractUniquenessStorage &classificator,
                                  std::unordered_map<const Edge *, CompactPath> &unique_extensions) {
    for(Edge *edgeIt : uniqueEdges) {
        Edge &edge = *edgeIt;
        if(unique_extensions.find(&edge) != unique_extensions.end())
            continue;
        const VertexRecord &rec = reads_storage.getRecord(edge.getStart());
        dbg::GraphPath path = greedyExtension(rec, classificator, edge);
        VERIFY(edge == path.frontEdge());
        path = findBulgeExtension(rec, edge, CompactPath(path)).unpack();
        VERIFY(edge == path.frontEdge());
        for(size_t i = 1; i < path.size(); i++) {
            if(classificator.isUnique(path.frontEdge())) {
                path = path.subPath(0, i + 1);
                break;
            }
        }
        if(path.size() == 1) {
            continue;
        }
        VERIFY(edge == path.frontEdge());
        unique_extensions.emplace(&edge, CompactPath::Subpath(path,1, path.size()));
        Edge &last_rc_edge = path.backEdge().rc();
        if(classificator.isUnique(last_rc_edge) && unique_extensions.find(&last_rc_edge) == unique_extensions.end()) {
            unique_extensions.emplace(&last_rc_edge, CompactPath::Subpath(path.RC(), 1, path.size()));
        }
    }
}

inline std::unordered_map<const Edge *, CompactPath> constructUniqueExtensions(logging::Logger &logger,
                                                                               SparseDBG &dbg, const RecordStorage &reads_storage,
                                                                               const AbstractUniquenessStorage &classificator) {
    std::unordered_map<const Edge *, CompactPath> unique_extensions;
    std::vector<Edge*> uniqueEdges;
    for(Edge &edge : dbg.edges()) {
        if (classificator.isUnique(edge))
            uniqueEdges.push_back(&edge);
    }
    struct {
        bool operator()(Edge* a, Edge* b) const {
            if(a == b)
                return false;
            if((a->truncSize() < 10000 && a->getCoverage() < 3) || (b->truncSize() < 10000 && b->getCoverage() < 3)) {
                if(a->intCov() * b->truncSize() != b->intCov() * a->truncSize())
                    return a->intCov() * b->truncSize() > b->intCov() * a->truncSize();
            } else if(a->truncSize() != b->truncSize())
                return a->truncSize() > b->truncSize();
            return *a < *b;
        }
    } customLess;
    std::sort(uniqueEdges.begin(), uniqueEdges.end(), customLess);
    findEasyExtensions(uniqueEdges, reads_storage, classificator, unique_extensions);
    findComplexExtensions(uniqueEdges, reads_storage, classificator, unique_extensions);
    return std::move(unique_extensions);
}

dbg::GraphPath correctRead(std::unordered_map<const Edge *, CompactPath> &unique_extensions,
                         const dbg::GraphPath &initial_al) {
    CompactPath initialCompactPath(initial_al);
    dbg::GraphPath al = initial_al;
    bool bad;
    bool corrected = false;
    for(size_t i = 0; i + 1 < al.size(); i++) {
        if(unique_extensions.find(&al[i].contig()) == unique_extensions.end())
            continue;
        CompactPath &compactPath = unique_extensions.find(&al[i].contig())->second;
        if(compactPath.cpath().nonContradicts(CompactPath::Subpath(al, i + 1, al.size()).cpath()))
            continue;
        corrected = true;
        dbg::GraphPath new_al = al.subPath(0, i + 1);
        size_t corrected_len = al.subPath(i + 1, al.size()).truncLen();
        dbg::GraphPath replacement = compactPath.unpack();
        while(replacement.truncLen() < corrected_len &&
              unique_extensions.find(&replacement.back().contig()) != unique_extensions.end()) {
            replacement += unique_extensions[&replacement.back().contig()].unpack();
        }
        if(replacement.truncLen() < corrected_len) {
            size_t deficite = corrected_len - replacement.truncLen();
//            logger.info() << "Need to correct more than known " << read_id << "\n"
//                          << CompactPath(al.subalignment(i + 1, al.size())) << "\n" << compactPath << std::endl;
            new_al += replacement;
            while(new_al.finish().outDeg() == 1 && deficite > 0) {
                size_t len = std::min(deficite, new_al.finish().front().truncSize());
                new_al += Segment<Edge>(new_al.finish().front(), 0, len);
                deficite -= len;
            }
            bad = true;
        } else {
            for (const Segment<Edge> rep_seg : replacement) {
                if (corrected_len <= rep_seg.size()) {
                    new_al += rep_seg.shrinkRightBy(rep_seg.size() - corrected_len);
                    corrected_len = 0;
                    break;
                } else {
                    new_al += rep_seg;
                    corrected_len -= rep_seg.size();
                }
            }
        }
        al = new_al;
    }
    if(corrected)
        return std::move(al);
    else
        return initial_al;
}

void correctReads(logging::Logger &logger, size_t threads, RecordStorage &reads_storage,
                  std::unordered_map<const Edge *, CompactPath> &unique_extensions) {
    omp_set_num_threads(threads);
    logger.info() << "Correcting reads using unique edge extensions" << std::endl;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(reads_storage, unique_extensions)
    for(size_t i = 0; i < reads_storage.size(); i++) {
        AlignedRead &alignedRead = reads_storage[i];
        if(!alignedRead.valid())
            continue;
        const dbg::GraphPath al = alignedRead.path.unpack();
        if(al.size() > 1) {
            dbg::GraphPath corrected1 = correctRead(unique_extensions, al);
            dbg::GraphPath corrected2 = correctRead(unique_extensions, corrected1.RC()).RC();
            if(al != corrected2) {
                reads_storage.reroute(alignedRead, al, corrected2, "mult correction");
            }
        }
    }
    reads_storage.applyCorrections(logger, threads);
}

void CorrectBasedOnUnique(logging::Logger &logger, size_t threads, SparseDBG &sdbg, RecordStorage &reads_storage,
                          const AbstractUniquenessStorage &classificator, const std::experimental::filesystem::path &ext_file) {
    std::unordered_map<const Edge *, CompactPath> unique_extensions =
            constructUniqueExtensions(logger, sdbg, reads_storage, classificator);
    std::ofstream os;
    os.open(ext_file);
    for(auto &it : unique_extensions) {
        os << it.first->getInnerId() << " " << it.second.cpath() << "\n";
    }
    os.close();
    correctReads(logger, threads, reads_storage, unique_extensions);
    logger.info() << "Collecting bad edges" << std::endl;
    std::unordered_set<Edge const *> bad_edges;
    for(Edge & edge : sdbg.edgesUnique()) {
        if(edge.innerSize() > 5000)
            continue;
        if(reads_storage.getRecord(edge.getStart()).isDisconnected(edge) ||
           reads_storage.getRecord(edge.rc().getStart()).isDisconnected(edge.rc())) {
            bad_edges.emplace(&edge);
            bad_edges.emplace(&edge.rc());
        }
    }
    logger.info() << "Removed " << bad_edges.size() / 2 << " disconnected edges"<< std::endl;
    std::ofstream brs;
    std::function<bool(const Edge&)> is_bad = [&bad_edges](const Edge &edge) {
        return edge.getCoverage() < 2 || bad_edges.find(&edge) != bad_edges.end();
    };
    reads_storage.delayedInvalidateBad(logger, threads, is_bad, "after_mult");
    reads_storage.applyCorrections(logger, threads);
}

SetUniquenessStorage PathUniquenessClassifier(logging::Logger &logger, size_t threads, SparseDBG &dbg, RecordStorage &reads_storage,
                                              const AbstractUniquenessStorage &classificator) {
    logger.info() << "Looking for more unique edges" << std::endl;
    SetUniquenessStorage res;
    for(Edge &edge : dbg.edges()) {
        if(classificator.isUnique(edge)) {
            res.addUnique(edge);
            continue;
        }
        const VertexRecord &rec = reads_storage.getRecord(edge.getStart());
        CompactPath unique_extension = rec.getFullUniqueExtension(edge.truncSeq().Subseq(0, 1), 1, 0);
        dbg::GraphPath path = unique_extension.unpack();
        size_t len = 0;
        for(size_t i = 1; i < path.size(); i++) {
            if(classificator.isUnique(path[i].contig())) {
                if(len < 3000 && rec.countStartsWith(CompactPath::Subpath(path,0, i + 1).cpath()) >= 4) {
                    if(i == 1 && edge.getStart().inDeg() == 2 && edge.getFinish().outDeg() == 2 &&
                            edge.getStart().outDeg() == 1 && edge.getFinish().inDeg() == 1) {
                        if(classificator.isUnique(edge.getStart().rc().front()) && classificator.isUnique(edge.getStart().rc().back()) &&
                           classificator.isUnique(edge.getFinish().front()) && classificator.isUnique(edge.getFinish().back())) {
                            continue;
                        }
                    }
                    res.addUnique(edge);
                    logger.trace() << "Found extra unique edge " << edge.getInnerId() << " " << edge.truncSize() << " " << edge.getCoverage() << std::endl;
                    break;
                }
            }
            len += path[i].size();
        }
    }
    logger.info() << "Finished unique edges search. Found " << res.size() << " unique edges" << std::endl;
    return std::move(res);
}

void DrawMult(const std::experimental::filesystem::path &dir, dbg::SparseDBG &dbg, size_t unique_threshold,
              RecordStorage &reads_storage, AbstractUniquenessStorage &uniquenessStorage) {
    std::vector<Component> split = LengthSplitter(unique_threshold).splitGraph(dbg);
    recreate_dir(dir);
    const std::function<std::string(Edge &)> colorer = [&uniquenessStorage](Edge &edge) {
        if(uniquenessStorage.isUnique(edge))
            return "black";
        if(uniquenessStorage.isError(edge))
            return "red";
        if(!edge.is_reliable)
            return "orange";
        return "blue";
    };
    for(size_t i = 0; i < split.size(); i++) {
        printDot(dir / (itos(i) + ".dot"), split[i], reads_storage.labeler(), colorer);
    }
}

RecordStorage MultCorrect(logging::Logger &logger, size_t threads, SparseDBG &dbg, const std::experimental::filesystem::path &dir,
                          RecordStorage &reads_storage, size_t unique_threshold, double initial_rel_coverage, bool diploid,
                          bool debug) {
    if(debug) {
        recreate_dir(dir);
    }
    UniqueClassificator classificator(dbg, reads_storage, initial_rel_coverage, diploid, debug);
    classificator.classify(logger, unique_threshold, dir/"ongoing");
    if(debug)
        DrawMult(dir / "round1", dbg, unique_threshold, reads_storage, classificator);
    CorrectBasedOnUnique(logger, threads, dbg, reads_storage, classificator, dir/"round1.txt");
    SetUniquenessStorage more_unique = PathUniquenessClassifier(logger, threads, dbg, reads_storage, classificator);
    SetUniquenessStorage more_more_unique = PathUniquenessClassifier(logger, threads, dbg, reads_storage, more_unique);
    if(debug)
        DrawMult(dir / "round2", dbg, unique_threshold, reads_storage, more_more_unique);
    CorrectBasedOnUnique(logger, threads, dbg, reads_storage, more_more_unique, dir/"round2.txt");
    if(debug)
        DrawMult(dir / "final", dbg, unique_threshold, reads_storage, more_unique);
    auto res = std::move(ResolveLoops(logger, threads, dbg, reads_storage, more_unique));
    return std::move(res);
}
