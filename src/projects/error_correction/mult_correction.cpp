#include "mult_correction.hpp"
#include "correction_utils.hpp"
#include "read_cleaning.hpp"
#include "assembly_graph/visualization.hpp"

using namespace dbg;
using namespace ag;
void printAl(logging::Logger &logger, std::unordered_map<const dbg::Edge *, GraphPath> &unique_extensions,
             const ag::GraphPath &al) {
    for(Edge &piece : al.edges()) {
        logger << piece.str() << " ";
        if(unique_extensions.find(&piece) != unique_extensions.end()) {
            logger << "+ ";
        }
    }
    logger << std::endl;
}

struct UEdge {
    UEdge(dbg::Edge *from, dbg::Edge *to, ag::GraphPath cpath, size_t support) : from(from), to(to), cpath(std::move(cpath)),
                                                                                    support(support) {}

    dbg::Edge *from;
    dbg::Edge *to;
    ag::GraphPath cpath;
    size_t support;

    UEdge RC() const {
        return {&to->rc(), &from->rc(), cpath.RC(), support};
    }
};
//std::unordered_map<const Edge *, GraphPath> constructUniqueExtensions(logging::Logger &logger, SparseDBG &dbg,
//                                                                         const dbg::ReadAlignmentStorage &reads_storage, const UniqueClassificator &classificator) {
//    std::unordered_map<Edge *, std::vector<UEdge>> bg;
//    for(Edge &edge : dbg.edges()) {
//        if(!classificator.isUnique(edge))
//            continue;
//        Vertex &start = *edge.getStart();
//        std::vector<Sequence> extensions;
//        for(auto & c : reads_storage.getRecord(start)) {
//            GraphAlignment al = GraphPath(start, c.first).getAlignment();
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
//            extensions.emplace_back(GraphPath(al).cpath());
//        }
//        std::sort(extensions.begin(), extensions.end());
//        extensions.erase(std::unique(extensions.begin(), extensions.end()), extensions.end());
//        for(Sequence &extension : extensions) {
//            GraphPath new_path(start, extension);
//            bg[&edge].emplace_back(&edge, &new_path.getAlignment().back().contig(), new_path,
//                                   reads_storage.getRecord(start).countStartsWith(extension));
//        }
//    }
//    std::unordered_map<Edge *, UEdge> choice;
//}

inline void findEasyExtensions(const std::vector<Edge *> &uniqueEdges, const dbg::DBGAlignedReadStorage &reads_storage,
                               const AbstractUniquenessStorage &classificator,
                               std::unordered_map<const Edge *, GraphPath> &unique_extensions) {
    VERIFY(uniqueEdges.empty() || uniqueEdges.front()->truncSize() >= uniqueEdges.back()->truncSize());
    for(Edge *edgeIt : uniqueEdges) {
        Edge &edge = *edgeIt;
        if (unique_extensions.find(&edge) != unique_extensions.end())
            continue;
        Vertex & start = edge.getStart();
        const ag::SuffixRecord &erec = reads_storage.getSuffixes().getSuffixRecord(edge);
//        Sequence seq = edge.truncSeq().Subseq(0, 1);
        GraphPath al = FullSuffixSupportedExtension(erec, GraphPath(edge.getFinish()), 1, 0);
        if(al.empty())
            continue;
        for(PathPosition pp = al.firstPosition(); pp != al.lastPosition(); ++pp) {
            if(classificator.isUnique(pp.nextEdge())) {
                al = al.subPath(al.firstPosition(), pp + 1);
                break;
            }
        }
        VERIFY(!al.empty());
        if(!classificator.isUnique(al.back().contig()))
            continue;
        al.push_front(edge);
        unique_extensions.emplace(&al.frontEdge(), al.subPath(al.firstPosition() + 1));
        al = al.RC();
        unique_extensions.emplace(&al.frontEdge(), al.subPath(al.firstPosition() + 1));
    }
}

ag::GraphPath greedyExtension(const ag::SuffixRecord &rec, const AbstractUniquenessStorage &classificator, Edge &edge) {
    ag::GraphPath path(edge.getFinish());
    Sequence seq = edge.getCode();
    while(true) {
        size_t best_val = 0;
        Edge *next_edge = nullptr;
        for(Edge &next_candidate : path.getFinish()) {
            if(classificator.isError(next_candidate))
                continue;
            path += next_candidate;
            size_t val = rec.countStartsWith(path);
            if(val > best_val) {
                best_val = val;
                next_edge = &next_candidate;
            }
            path.pop_back();
        }
        if(best_val == 0)
            break;
        path += *next_edge;
    }
    return std::move(path);
}

PathPosition findVertexInPath(const GraphPath &path, Vertex &vertex, PathPosition pp) {
    while(pp != path.endPosition()) {
        if(vertex == pp.getVertex())
            return pp;
        ++pp;
    }
    return pp;
}

//Dealing with the case when unique edge is folowed by a fork. Choosing one direction for this case
inline GraphPath findBulgeExtension(const ag::SuffixRecord &rec, Edge &edge, const GraphPath & greedy) {
    if(edge.getFinish().outDeg() != 2)
        return greedy;
    Edge &edge1 = edge.getFinish().front();
    Edge &edge2 = edge.getFinish().back();
    GraphPath p1 = GraphPath(edge1) + FullSuffixSupportedExtension(rec, GraphPath(edge1), 1, 0);
    GraphPath p2 = GraphPath(edge2) + FullSuffixSupportedExtension(rec, GraphPath(edge2), 1, 0);
//    We make sure that greedy extension corresponds to the first extension
    if(greedy.startsWith(p2)) {
        std::swap(p1, p2);
    } else {
//        If gready extension does not correspond to any of the extensions we drop the analysis
        if(!greedy.startsWith(p1))
            return greedy;
    }
    p1 = greedy;
    PathPosition b1 = p1.firstPosition();
    PathPosition b2 = p2.firstPosition();
    GraphPath choice;
//    Considering the most simple case where one of the edges in the forks forms a bulge with alternative extension
    Vertex &p1v2 = (p1.firstPosition() + 1).getVertex();
    Vertex &p2v2 = (p2.firstPosition() + 1).getVertex();
    if(findVertexInPath(p2, p1v2, p2.firstPosition() + 1) != p2.endPosition()) {
        b1 = p1.firstPosition() + 1;
        b2 = findVertexInPath(p2, p1v2, p2.firstPosition() + 1);
        if(findVertexInPath(p1, p2v2, p1.firstPosition() + 1) != p1.endPosition()) {
            return greedy;
        }
        choice = p2.subPath(p2.firstPosition(), b2);
    } else if(findVertexInPath(p1, p2v2, p1.firstPosition() + 1) != p1.endPosition()) {
        b1 = findVertexInPath(p1, p2v2, p1.firstPosition() + 1);
        b2 = p2.firstPosition() + 1;
        choice = p1.subPath(p1.firstPosition(), b1);
    }
//    choice contains the alternative path in the bulge
//    If after the bulge extentions not contradict each other we drop the analysis.
    if(!p1.subPath(b1).nonContradicts(p2.subPath(b2)))
        return greedy;
//    We correct greedy path to go through the more complex path.
    return choice + p1.subPath(b1);
}

inline void findComplexExtensions(const std::vector<Edge *> &uniqueEdges, const dbg::DBGAlignedReadStorage &reads_storage,
                                  const AbstractUniquenessStorage &classificator,
                                  std::unordered_map<const Edge *, GraphPath> &unique_extensions) {
    for(Edge *edgeIt : uniqueEdges) {
        Edge &edge = *edgeIt;
        if(unique_extensions.find(&edge) != unique_extensions.end())
            continue;
        const ag::SuffixRecord &rec = reads_storage.getSuffixes().getSuffixRecord(edge);
        ag::GraphPath path = greedyExtension(rec, classificator, edge);
        VERIFY(edge.getFinish() == path.getStart());
        path = findBulgeExtension(rec, edge, path);
        VERIFY(edge.getFinish() == path.getStart());
        for(PathPosition pp = path.firstPosition(); pp != path.lastPosition(); ++pp) {
            if(classificator.isUnique(pp.nextEdge())) {
                path = path.subPath(path.firstPosition(), pp + 1);
                break;
            }
        }
        if(path.empty()) {
            continue;
        }
        VERIFY(edge.getFinish() == path.getStart());
        unique_extensions.emplace(&edge, path);
        Edge &last_rc_edge = path.backEdge().rc();
        if(classificator.isUnique(last_rc_edge) && unique_extensions.find(&last_rc_edge) == unique_extensions.end()) {
            path.push_front(edge);
            path.pop_back();
            unique_extensions.emplace(&last_rc_edge, path.RC());
        }
    }
}

inline std::unordered_map<const Edge *, GraphPath> constructUniqueExtensions(logging::Logger &logger,
                                                                               SparseDBG &dbg, const dbg::DBGAlignedReadStorage &reads_storage,
                                                                               const AbstractUniquenessStorage &classificator) {
    std::unordered_map<const Edge *, GraphPath> unique_extensions;
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

//This procedure should not exist int this world
ag::GraphPath correctRead(std::unordered_map<const Edge *, GraphPath> &unique_extensions,
                         const ag::GraphPath &initial_al) {
    GraphPath initialGraphPath(initial_al);
    ag::GraphPath al = initial_al;
    bool bad;
    bool corrected = false;
    for(PathPosition cur = al.firstPosition(); cur + 1 != al.lastPosition(); ++cur) {
        if(unique_extensions.find(&cur.nextEdge()) == unique_extensions.end())
            continue;
        GraphPath replacement = unique_extensions.find(&cur.nextEdge())->second;
        if(replacement.nonContradicts(al.subPath(cur + 1)))
            continue;
        corrected = true;
        ag::GraphPath new_al = al.subPath(al.firstPosition(), cur + 1);
        size_t corrected_len = al.subPath(cur + 1).truncLen();
        while(replacement.truncLen() < corrected_len &&
              unique_extensions.find(&replacement.back().contig()) != unique_extensions.end()) {
            replacement += unique_extensions[&replacement.back().contig()];
        }
        if(replacement.truncLen() < corrected_len) {
            size_t deficite = corrected_len - replacement.truncLen();
            new_al += replacement;
            while(new_al.getFinish().outDeg() == 1 && deficite > 0) {
                size_t len = std::min(deficite, new_al.getFinish().front().truncSize());
                new_al += Segment<Edge>(new_al.getFinish().front(), 0, len);
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
        break;
    }
    if(corrected)
        return std::move(al);
    else
        return initial_al;
}

void correctReads(logging::Logger &logger, size_t threads, dbg::DBGAlignedReadStorage &reads_storage,
                  std::unordered_map<const Edge *, GraphPath> &unique_extensions) {
    omp_set_num_threads(threads);
    logger.info() << "Correcting reads using unique edge extensions" << std::endl;
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(reads_storage, unique_extensions)
    for(size_t i = 0; i < reads_storage.getReads().size(); i++) {
        ag::AlignedRead &alignedRead = reads_storage.getReads()[i];
        if(!alignedRead.valid())
            continue;
        const ag::GraphPath al = alignedRead.getPath();
        if(!al.isSingleton()) {
            ag::GraphPath corrected1 = correctRead(unique_extensions, al);
            ag::GraphPath corrected2 = correctRead(unique_extensions, corrected1.RC()).RC();
            if(al != corrected2) {
                reads_storage.getReads().rerouteRead(alignedRead, corrected2, "mult correction");
            }
        }
    }
    reads_storage.getReads().applyCorrections(logger, threads);
}

void CorrectBasedOnUnique(logging::Logger &logger, size_t threads, SparseDBG &sdbg, DBGAlignedReadStorage &reads_storage,
                          const AbstractUniquenessStorage &classificator, const std::experimental::filesystem::path &ext_file) {
    std::unordered_map<const Edge *, GraphPath> unique_extensions =
            constructUniqueExtensions(logger, sdbg, reads_storage, classificator);
    std::ofstream os;
    os.open(ext_file);
    for(auto &it : unique_extensions) {
        os << it.first->getInnerId() << " " << it.second << "\n";
    }
    os.close();
    correctReads(logger, threads, reads_storage, unique_extensions);
    logger.info() << "Collecting bad edges" << std::endl;
    std::unordered_set<Edge const *> bad_edges;
    for(Edge & edge : sdbg.edgesUnique()) {
        if(edge.innerSize() > 5000)
            continue;
        if(reads_storage.getSuffixes().getSuffixRecord(edge).empty() ||
           reads_storage.getSuffixes().getSuffixRecord(edge.rc()).empty()) {
            bad_edges.emplace(&edge);
            bad_edges.emplace(&edge.rc());
        }
    }
    logger.info() << "Removed " << bad_edges.size() / 2 << " disconnected edges"<< std::endl;
    std::ofstream brs;
    std::function<bool(const Edge&)> is_bad = [&bad_edges](const Edge &edge) {
        return edge.getCoverage() < 2 || bad_edges.find(&edge) != bad_edges.end();
    };
    InvalidateBad(logger, threads, reads_storage.getReads(), 500, is_bad, "after_mult");
    reads_storage.getReads().applyCorrections(logger, threads);
}

SetUniquenessStorage PathUniquenessClassifier(logging::Logger &logger, size_t threads, SparseDBG &dbg, DBGAlignedReadStorage &reads_storage,
                                              const AbstractUniquenessStorage &classificator) {
    logger.info() << "Looking for more unique edges" << std::endl;
    SetUniquenessStorage res;
    for(Edge &edge : dbg.edges()) {
        if(classificator.isUnique(edge)) {
            res.addUnique(edge);
            continue;
        }
        const ag::SuffixRecord &rec = reads_storage.getSuffixes().getSuffixRecord(edge);
        GraphPath path = GraphPath(edge) + FullSuffixSupportedExtension(rec, GraphPath(edge.getFinish()), 1, 0);
        size_t len = 0;
        for(PathPosition pos = path.firstPosition() + 1; pos != path.lastPosition(); ++pos) {
            if(classificator.isUnique(pos.nextEdge())) {
                if(len < 3000 && rec.countStartsWith(path.subPath(path.firstPosition() + 1, pos + 1)) >= 4) {
                    if(pos == path.firstPosition() + 1 && edge.getStart().inDeg() == 2 && edge.getFinish().outDeg() == 2 &&
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
            len += path.getSegment(pos).size();
        }
    }
    logger.info() << "Finished unique edges search. Found " << res.size() << " unique edges" << std::endl;
    return std::move(res);
}

void DrawMult(const std::experimental::filesystem::path &dir, dbg::SparseDBG &dbg, size_t unique_threshold,
              DBGAlignedReadStorage &reads_storage, AbstractUniquenessStorage &uniquenessStorage) {
    std::vector<Component> split = ag::LengthSplitter(unique_threshold).splitGraph(dbg);
    recreate_dir(dir);
    const std::function<std::string(const Edge &)> colorer = [&uniquenessStorage](const Edge &edge) {
        if(uniquenessStorage.isUnique(edge))
            return "black";
        if(uniquenessStorage.isError(edge))
            return "red";
        if(!edge.is_reliable)
            return "orange";
        return "blue";
    };
    Printer printer(ObjInfo<dbg::Edge>({reads_storage.getSuffixes().labeler()}, {colorer}, {}));
    for(size_t i = 0; i < split.size(); i++) {
        // printDot(dir / (itos(i) + ".dot"), split[i], reads_storage.labeler(), colorer); delete if ok
        printer.printDot(dir / (itos(i) + ".dot"), split[i]);
    }
}

ag::AlignedReadStorage MultCorrect(logging::Logger &logger, size_t threads, SparseDBG &dbg, const std::experimental::filesystem::path &dir,
                                       DBGAlignedReadStorage &reads_storage, size_t unique_threshold, double initial_rel_coverage, bool diploid,
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
    for(Edge &edge: dbg.edges()) edge.is_reliable = false;
    return std::move(res);
}
