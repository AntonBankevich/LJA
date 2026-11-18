#include "read_alignment_storage.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include <map>
#include <tuple>
using namespace ag;

namespace {
//    Records are unordered: removes dir (which must be there) by moving the last entry into its place.
    void RemoveDirection(std::vector<AlignedReadDirection> &record, const AlignedReadDirection &dir) {
        auto it = std::find(record.begin(), record.end(), dir);
        VERIFY(it != record.end());
        *it = record.back();
        record.pop_back();
    }
}

void AlignedReadStorageMaintenance::fireResolveVertex(Vertex &core, const VertexResolutionResult  &resolution) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    std::vector<AlignedReadDirection> core_starts = std::move(storage->getOutgoingReadsLockFree(core));
    typedef std::pair<std::vector<AlignedReadDirection> *, std::vector<AlignedReadDirection> *> DestRecs;
    std::unordered_map<VertexId, DestRecs> dest_recs;
    for (Edge &e : core) {
        Vertex &d = e.getFinish();
        dest_recs.emplace(d.getId(), DestRecs{&storage->getOutgoingReadsLockFree(d),
                                              &storage->getSubstringReadsLockFree(d.getId())});
        dest_recs.emplace(d.rc().getId(), DestRecs{&storage->getOutgoingReadsLockFree(d.rc()),
                                                   &storage->getSubstringReadsLockFree(d.rc().getId())});
    }
    for (AlignedReadDirection &dir : core_starts) {
        AlignedRead &read = dir.getRead();
        read.lock();
        Edge &prefix = dir.frontEdge();
        VERIFY(prefix.isPrefix());
        Vertex &dest = prefix.getFinish();
        dir.pop_front(prefix);
        if (dir.empty()) {
//            The read collapses to a substring of dest. Both directions share the same path, so this pop already
//            emptied the opposite view too -- but the opposite direction is still registered in starts[dest.rc()],
//            a bucket that neither this fire nor the core.rc() fire iterates, so move it to the substring index
//            explicitly. dir may be either direction of the read, so the opposite one is dir.RC(), not backward().
            AlignedReadDirection bwd = dir.RC();
            DestRecs &rc_rec = dest_recs.at(dest.rc().getId());
            RemoveDirection(*rc_rec.first, bwd);
            rc_rec.second->emplace_back(bwd);
            dest_recs.at(dest.getId()).second->emplace_back(dir);
        } else {
            dest_recs.at(dest.getId()).first->emplace_back(dir);
        }
        read.unlock();
    }
}

void AlignedReadStorageMaintenance::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    Edge &inc = v.incFront();
    Edge &out = v.front();
    if (!e.getStart().coverage_info.empty()) {return;}
    std::vector<AlignedReadDirection> &old_recs = storage->getOutgoingReads(e.getStart());
    std::vector<AlignedReadDirection> filtered_old_recs;
    std::vector<AlignedReadDirection> &new_recs = storage->getOutgoingReads(v);
    std::vector<AlignedReadDirection> &new_subrecs = storage->getSubstringReads(v.getId());
//    Decide which directions leave on e before rewriting any: both directions of a read share one path, and when e goes
//    from u to u.rc(), e.rc() also starts at u, so the other direction of a read on e is in this bucket too. It
//    belongs to the mirrored call for e.rc() and must be recognised before this call rewrites the shared path.
    std::vector<AlignedReadDirection> on_e;
    for (AlignedReadDirection &dir : old_recs) {
        if (dir.getStart() == e.getStart() && dir.frontEdge() != e)
            filtered_old_recs.emplace_back(dir);
        else
            on_e.emplace_back(dir);
    }
    for (AlignedReadDirection &dir : on_e) {
//        A direction already converted by the call for the other strand: its read is now a subread of v.
        if (dir.getStart() == v) {
            new_subrecs.emplace_back(dir);
        } else {
            dir.pop_front(inc);
            if (dir.getFinish().coverage_info.empty() && dir.getFinish() == e.getFinish() && dir.getFSplits().size() == 0) {
                dir.pop_back(out);
                new_subrecs.emplace_back(dir);
            } else {
                new_recs.emplace_back(dir);
            }
        }
    }
    old_recs = std::move(filtered_old_recs);
    VERIFY(storage->getSubstringReads(e.getStart().getId()).empty());
}


void AlignedReadStorageMaintenance::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    size_t left_skip = 0;
    GraphPath prefix(path.frontEdge().getStart());
    Vertex & start = path.getStart();
    std::vector<AlignedReadDirection> &collected_recs = storage->getOutgoingReads(start);
    for(Edge &edge : path.edges()) {
        if (edge == path.backEdge())
            continue;
        prefix += edge;
        left_skip += edge.rc().truncSeq().size();
//        DBG only: interior k-mers have no subreads (reads are at least k + w long).
        VERIFY(storage->getSubstringReads(edge.getFinish().getId()).empty());
        std::vector<AlignedReadDirection> starting = std::move(storage->getOutgoingReads(edge.getFinish()));
        VERIFY(storage->getOutgoingReads(edge.getFinish()).empty());
        for(const AlignedReadDirection &dir : starting) {
            dir.getRead().lock();
            VERIFY(dir.valid())
            VERIFY(dir.getStart() == edge.getFinish());
            size_t cut_left = left_skip + dir.leftCut();
            dir.setCutLeft(0);
            dir.push_front(prefix);
            dir.setCutLeft(cut_left);
            collected_recs.emplace_back(dir);
            dir.getRead().unlock();
        }
    }
}

//    Only reads starting at a removed vertex (and subreads of one) change; reads entering P from outside keep their codes.
//    A read that ends within P, at a removed vertex or at D, lies within the dominant vertex D and becomes its subread.
//    Any other one goes on past P's finish (or ends at a surviving finish other than D). It is extended to P's start; if
//    P's start is replaced, its P part is exchanged for Q (they differ only by codes that moved into the replacements of
//    the start's in-edges, which these reads do not use), otherwise the read is already its form over Q, as P and Q have
//    the same codes. Then it is shortened to start at D by Q's edge into D. If it wraps around a loop P back to a removed
//    vertex, its end is converted as a start by the call for the other strand, or by this one if that call is skipped (P
//    is then its own rc, so its other direction starts at a removed vertex of P too). Only the read's start is changed,
//    without following codes (until removeMarked() old and new edges share first code letters), so deck positions and
//    suffix ends stay valid.
void AlignedReadStorageMaintenance::fireReplacePath(const PathReplacement &rep) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    const RAGraphPath &old_path = rep.oldPath();
    Vertex &dominant = rep.getDominant();
    const GraphPath whole(old_path);
    GraphPath new_whole(rep.newPath().getStart());
    for (Edge &e : rep.newPath().edges())
        new_whole += e;
    auto convert = [&](const AlignedReadDirection &dir, size_t i, const GraphPath &prefix) {
//        Compared with P's codes from the read's start on, a read with fewer ends inside P, at a removed vertex, and one
//        with exactly as many at P's finish.
        const size_t frest = whole.getFSplits().size() - prefix.getFSplits().size();
        const size_t rrest = whole.getRSplits().size() - prefix.getRSplits().size();
        const bool within = dir.getFSplits().size() <= frest && dir.getRSplits().size() <= rrest;
        const bool at_finish = dir.getFSplits().size() == frest && dir.getRSplits().size() == rrest;
        if (within && (!at_finish || rep.endReplaced() || old_path.getFinish() == dominant)) {
            const size_t from = rep.oldPosition(i) + dir.leftCut() - rep.dominantPosition();
            dir.setPath(GraphPath(dominant, from, dominant.size() - from - dir.getPath().len()));
//            It ended at D, which survives, so its other direction starts there: it becomes a subread too.
            if (at_finish && !rep.endReplaced()) {
                RemoveDirection(storage->getOutgoingReadsLockFree(dominant.rc()), dir.RC());
                storage->getSubstringReadsLockFree(dominant.rc().getId()).emplace_back(dir.RC());
            }
            return;
        }
        const size_t left = dir.leftCut() + rep.oldPosition(i);
        dir.setCutLeft(0);
        dir.push_front(prefix);
//        Q = [D]: P's codes went to the replacements of the edges into a replaced start, and P's rc codes stay on the
//        rc path's last edge if the finish is replaced too.
        if (rep.startReplaced() && rep.newPath().size() == 0)
            dir.absorb_front(whole, dominant, rep.endReplaced());
        else if (rep.startReplaced())
            dir.replace_front(whole, new_whole);
        dir.setCutLeft(left);
        if (dir.getStart() != dominant)
            dir.pop_front(rep.vpEdge());
        VERIFY(dir.getStart() == dominant && !dir.empty());
    };
//    A converted read starts at D; one that does not was converted by the call for the other strand already.
    auto refile = [&](std::vector<AlignedReadDirection> &directions, size_t i, const GraphPath &prefix) {
        for (AlignedReadDirection &dir : directions) {
            dir.getRead().lock();
            if (dir.getStart() != dominant)
                convert(dir, i, prefix);
            if (dir.empty())
                storage->getSubstringReadsLockFree(dir.getStart().getId()).emplace_back(dir);
            else
                storage->getOutgoingReadsLockFree(dir.getStart()).emplace_back(dir);
            dir.getRead().unlock();
        }
    };
//    prefix is P up to its i-th vertex, built by index.
    GraphPath prefix(old_path.getStart());
    for (size_t i = 0; i <= old_path.size(); i++) {
        if (!rep.survives(i)) {
            Vertex &v = old_path.getVertex(i);
            std::vector<AlignedReadDirection> starting = std::move(storage->getOutgoingReadsLockFree(v));
            storage->getOutgoingReadsLockFree(v).clear();
            std::vector<AlignedReadDirection> subreads = std::move(storage->getSubstringReadsLockFree(v.getId()));
            storage->getSubstringReadsLockFree(v.getId()).clear();
            refile(starting, i, prefix);
            refile(subreads, i, prefix);
        }
        if (i < old_path.size())
            prefix += old_path.getEdge(i);
    }
}

void AlignedReadStorageMaintenance::fireAddRead(const AlignedRead &read) {
    VERIFY_MSG(false, "New reads can not be added when maintenance is already activated");
}

void AlignedReadStorageMaintenance::fireRerouteRead(AlignedRead &read) {
    storage->updateStart(read.forward());
    storage->updateStart(read.backward());
}

void AlignedReadStorageMaintenance::fireInvalidateRead(AlignedRead &read) {
    storage->updateStart(read.forward());
    storage->updateStart(read.backward());
}

AlignedReadStorageMaintenance::AlignedReadStorageMaintenance(AssemblyGraph &graph,
                                                             AlignedReadStorage &storage) :
        AlignedReadStorageListener(storage, "AlignedReadStorageMaintenance"), ResolutionListener(graph, "AlignedReadStorageMaintenance"), storage(&storage) {
//        `storage` is not visible to any other thread until this constructor returns, so the whole
//        initial population below is done under a single lock_table() acquisition per map instead of
//        paying cuckoohash_map's per-key lock/unlock on every vertex/edge/read (unlike fireAddVertex/
//        fireAddEdge/getOutgoingReadsLockFree/getSubstringReadsLockFree, which stay per-key-locked
//        since they are also called on `storage` after construction, when concurrent access is real).
    {
        auto lt_inside = storage.reads_inside_vertices.lock_table();
        auto lt_starts = storage.starts.lock_table();
        for(Vertex &vertex : graph.vertices()) {
            lt_inside.insert(vertex.getId(), std::make_unique<std::vector<AlignedReadDirection>>());
            lt_starts.insert(vertex.getId(), std::make_unique<std::vector<AlignedReadDirection>>());
        }
    }
    {
        auto starts_lt = storage.starts.lock_table();
        auto vertices_lt = storage.reads_inside_vertices.lock_table();
        for(AlignedRead &read: storage) {
            if(!read.getPath().empty()) {
                starts_lt.at(read.getPath().getStart().getId())->emplace_back(read.forward());
                starts_lt.at(read.getPath().getFinish().rc().getId())->emplace_back(read.backward());
            } else if(read.valid()) {
                vertices_lt.at(read.getPath().getStart().getId())->emplace_back(read.forward());
                vertices_lt.at(read.getPath().getFinish().rc().getId())->emplace_back(read.backward());
            }
        }
    }
}

void AlignedReadStorageMaintenance::fireAddVertex(Vertex &vertex) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    storage->reads_inside_vertices.insert(vertex.getId(), std::make_unique<std::vector<AlignedReadDirection>>());
    storage->starts.insert(vertex.getId(), std::make_unique<std::vector<AlignedReadDirection>>());
}

void AlignedReadStorageMaintenance::fireDeleteVertex(Vertex &vertex) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    std::vector<AlignedReadDirection> subreads = std::move(storage->getSubstringReadsLockFree(vertex.getId()));
    storage->reads_inside_vertices.erase(vertex.getId());
//    In place: both directions share the path, and the call for vertex.rc() finds it legacy already.
    for (AlignedReadDirection &dir: subreads)
        dir.getPath().makeLegacy();
    if (!subreads.empty())
        storage->legacy_reads.insert(vertex.getId().legacyId(),
                                     std::make_unique<std::vector<AlignedReadDirection>>(std::move(subreads)));
//    Every operation moves the reads starting at a vertex it deletes before the vertex is removed.
    VERIFY(storage->getOutgoingReadsLockFree(vertex).empty());
    storage->starts.erase(vertex.getId());
}

void AlignedReadStorageMaintenance::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    VERIFY(split.size() > 1);
    VERIFY(!split.frontEdge().isSuffix());
    VERIFY(!split.backEdge().isPrefix());
    std::unordered_map<VertexId, std::vector<AlignedReadDirection> *> new_recs;
    for(Vertex &v : split.vertices())
        new_recs[v.getId()] = &storage->getOutgoingReadsRecord(v.getId());
    std::vector<AlignedReadDirection> start_recs = std::move(*new_recs.at(edge.getStart().getId()));
    for(AlignedReadDirection direction : start_recs) {
        VERIFY(direction.valid());
        VERIFY(!direction.empty());
        if (direction.frontEdge() != edge) {
//            This direction does not leave on the edge being split. Usually that means it still starts at
//            edge.getStart() and is simply put back. But the sibling (rc) fireSplitEdge call runs on the
//            same shared read path and may already have truncated this read onto an inner split piece,
//            moving its real start onto a new split-point vertex. Re-file under the direction's current
//            start (always one of split.vertices()) rather than blindly under edge.getStart().
            new_recs.at(direction.getStart().getId())->emplace_back(direction);
            continue;
        }
        // if(!direction.valid())
        //     continue;
        VertexId new_start;
//            This condition takes care of handling singleton paths when this listener is called for rc edge/path
        if(direction.getStart() == edge.getStart() && (direction.getFinish() != split.frontEdge().getFinish() ||
                                                        direction.getFSplits().size() > edge.getCode().size())) {
//                All directions extending beyond edge should be processed carefully to avoid spoiling NuclDeck
//                iterators stored in AlignedRead
            if (direction.getFSplits().size() > edge.getCode().size()) {
                VERIFY(direction.leftCut() >= edge.rc().truncSize() - split.backEdge().rc().truncSize());
                size_t left_cut = direction.leftCut();
                direction.setCutLeft(0);
                direction.replace_front(split.backEdge());
                direction.setCutLeft(left_cut + split.backEdge().rc().truncSize() - edge.rc().truncSize());
            } else {
                size_t to_skip = direction.leftCut();
                size_t start_pos = 0;
                for (Edge &e: split.edges()) {
                    if (to_skip < start_pos + e.rc().truncSize()) {
                        VERIFY(to_skip + direction.getPath().len() <= start_pos + e.fullSize());
                        if(to_skip + direction.getPath().len() <= start_pos + e.getStart().size()) {
                            VERIFY(false);
                            direction.invalidate();
//                            TODO: handle subreads
                        } else {
                            direction.setPath({e, to_skip - start_pos, start_pos + e.fullSize() - (to_skip + direction.getPath().len())});
                        }
                        break;
                    }
                    start_pos += e.rc().truncSize();
                }
            }
        }
        if(direction.valid()) {
            new_start = direction.getStart() == edge.getStart() ? split.getStart().getId() : direction.getStart().getId();
            new_recs[new_start]->emplace_back(direction);
        }
    }
}

//TODO: this operation can break subreads in supregraph. Change to merging to a new Vertex?
void AlignedReadStorageMaintenance::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                                        const AlignmentForm &left_al,
                                                        const AlignmentForm &right_al) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    VERIFY(!left.isSuffix());
    VERIFY(!right.isPrefix());
    VERIFY(left.getFinish().outDeg() == 0 && left.getFinish().inDeg() == 1);
    VERIFY(right.getStart().outDeg() == 1 && right.getStart().inDeg() == 0);
    VERIFY(left_al.targetLength() == right_al.targetLength());
    size_t left_skip = left.fullSize() - left_al.queryLength();
    size_t right_skip = right.fullSize() - right_al.queryLength();
    std::vector<AlignedReadDirection> &new_rec = storage->getOutgoingReadsRecord(new_edge.getStart().getId());
//        Reads on left edge will be handled by rc. The tip vertex is a k-mer, so it has no subreads.
    VERIFY(storage->getSubstringReads(right.getStart().getId()).empty());
    std::vector<AlignedReadDirection> starting = std::move(storage->getOutgoingReadsRecord(right.getStart().getId()));
    VERIFY(storage->getOutgoingReadsRecord(right.getStart().getId()).empty());
    for(AlignedReadDirection dir : starting) {
        dir.getRead().lock();
        size_t new_left = dir.leftCut() >= right_al.queryLength() ?
                          left_skip + right_al.targetLength() + dir.leftCut() - right_al.queryLength() :
                          left_skip + right_al.lastColumnByQpos(dir.leftCut()).getTpos();
//        A read only on right ends where right ends, so its other direction stays filed at right.getFinish().rc().
        if(dir.isSingleton()) {
            size_t old_right = right.fullSize() - dir.rightCut();
            size_t new_right = old_right >= right_al.queryLength() ?
                               left_skip + right_al.targetLength() + old_right - right_al.queryLength() :
                               left_skip + right_al.firstColumnByQpos(old_right).getTpos();
            dir.setCutRight(new_edge.fullSize() - new_right);
        }
        dir.setCutLeft(0);
        dir.forcePushFront(left);// While this operation is underway this path is disconnected.
        dir.setCutLeft(new_left);
        dir.getRead().unlock();
        new_rec.emplace_back(dir);
    }
}

void AlignedReadStorageMaintenance::fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph & graph) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
#pragma omp parallel for default(none) schedule(dynamic, 100)
    for(size_t i = 0; i < storage->size(); i++) {
        AlignedRead &read = (*storage)[i];
        read.resetEdgeCodes();
    }
}

bool AlignedReadStorageMaintenance::checkConsistency(std::ostream &os) {
    const AssemblyGraph &graph = getGraph();
    size_t errors = 0;
    auto report = [&]() -> std::ostream & {
        errors++;
        return os << "Read storage inconsistency: ";
    };
    for (const Vertex &v : graph.vertices()) {
        if (!storage->starts.contains(v.getId()))
            report() << "vertex " << v.getInnerId() << " has no starts record\n";
        if (!storage->reads_inside_vertices.contains(v.getId()))
            report() << "vertex " << v.getInnerId() << " has no reads_inside_vertices record\n";
    }
//    Given every live vertex has its records, extra records can only belong to deleted vertices.
    if (storage->starts.size() > graph.size())
        report() << storage->starts.size() - graph.size() << " starts records belong to deleted vertices\n";
    if (storage->reads_inside_vertices.size() > graph.size())
        report() << storage->reads_inside_vertices.size() - graph.size()
                 << " reads_inside_vertices records belong to deleted vertices\n";
    bool ok = storage->checkRecords(os);
    return errors == 0 && ok;
}

void AlignedReadStorageMaintenance::fireMergeLoop(const GraphPath &path, Vertex &new_vertex) {
    VERIFY_MSG(!storage->inReadCorrection(), "The graph can not change while read corrections are pending");
    VERIFY(false);
    // for (Edge &edge : path.edges()) {
    //     if (!edge.isPrefix()) {
    //         for (AlignedReadDirection &dir: this->storage->getOutgoingReads(edge)) {
    //             storage->delayedInvalidateRead(dir.getRead(), "Merge loop");
    //             storage->apply(dir.getRead());
    //         }
    //     }
    // }
//            TODO: handle subreads here
}

void AlignedReadStorage::delayedInvalidateRead(AlignedRead &read, const string &message) { // NOLINT(readability-convert-member-functions-to-static)
    VERIFY(read_correction);
    read.delayedInvalidate();
    this->fireDelayedInvalidateRead(read, message);
}

void AlignedReadStorage::rerouteRead(AlignedRead &alignedRead, GraphPath corrected,
                                     const string &message) {
//    VERIFY(corrected.truncLen() >= 500);
    VERIFY(read_correction);
    alignedRead.correct(Normalize(std::move(corrected)));
    this->fireDelayedRerouteRead(alignedRead, message);
}

bool AlignedReadStorage::apply(AlignedRead &alignedRead) {
    VERIFY(read_correction);
    if (!alignedRead.checkCorrected())
        return false;
    if(alignedRead.getCorrected().valid())
        this->fireRerouteRead(alignedRead);
    else
        this->fireInvalidateRead(alignedRead);
    alignedRead.applyCorrection();
    return true;
}

void
AlignedReadStorage::printReadFasta(logging::Logger &logger, const std::experimental::filesystem::path &path) const {
    logger.info() << "Printing reads to fasta file " << path << std::endl;
    std::string acgt = "ACGT";
    std::ofstream os;
    os.open(path);
    for (const AlignedRead &read: reads) {
        const GraphPath &al = read.getPath();
        if (!al.valid() || al.isLegacy())
            continue;
        os << ">" << read.getId() << "\n" << read.getPath().Seq() << "\n";
    }
    os.close();
}

void AlignedReadStorage::printReadPaths(logging::Logger &logger,
                                        const std::experimental::filesystem::path &aln_path,
                                        const std::experimental::filesystem::path &gfa_path,
                                        const std::experimental::filesystem::path &rp_path,
                                        size_t k) const {
    logger.info() << "Printing reads paths to file " << aln_path << std::endl;
    std::ofstream os;
    os.open(aln_path);
    Save(os);
    os.close();
    os.open(rp_path);
    os << gfa_path.c_str() << std::endl;
    os << aln_path.c_str() << std::endl;
    os << k << std::endl;
    os.close();
}

void AlignedReadStorage::printFullAlignments(logging::Logger &logger,
                                             const std::experimental::filesystem::path &path) const {
    logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
    std::ofstream os;
    os.open(path);
    for (const AlignedRead &read: reads) {
        const GraphPath &al = read.getPath();
        if (!al.valid())
            continue;
        os << read.getId() << " " << read.getPath().str() << "\n";
        os << "-" << read.getId() << " " << read.getPath().RC().str() << "\n";
    }
    os.close();
}

void AlignedReadStorage::printSequences(const std::experimental::filesystem::path &path) const {
    std::ofstream os;
    os.open(path);
    for (const AlignedRead &read: reads) {
        const GraphPath &al = read.getPath();
        if (!al.valid())
            continue;
        os << ">" << read.getId() << "\n" << al.Seq() << "\n";
    }
    os.close();
}

void AlignedReadStorage::applyCorrections(logging::Logger &logger, size_t threads) {
    VERIFY(read_correction);
    if (size() > 10000)
        logger.info() << "Applying corrections to reads" << std::endl;
    omp_set_num_threads(int(threads));
    ParallelCounter cnt(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(cnt)
    for (size_t i = 0; i < reads.size(); i++) { // NOLINT(modernize-loop-convert)
        if (apply(reads[i]))
            cnt += 1;
    }
    this->fireAppliedCorrections(cnt.get());
    read_correction = false;
    if (size() > 10000)
        logger.info() << "Applied correction to " << cnt.get() << " reads" << std::endl;
}

void AlignedReadStorage::Save(std::ostream &os) const {
    os << size() << "\n";
    for (const AlignedRead &alignedRead: *this) {
        os << alignedRead << "\n";
    }
}

AlignedReadStorage AlignedReadStorage::Load(logging::Logger &logger, size_t threads,
                                            std::istream &is, AssemblyGraph &graph) {
    IdIndex<Vertex> index(graph.vertices().begin(), graph.vertices().end());
    return {logger, threads, graph, LoadReadAlignments(is, index)};
}

std::vector<AlignedRead> AlignedReadStorage::LoadReadAlignments(std::istream &is, IdIndex<Vertex> &index) {
    size_t sz;
    is >> sz;
    std::vector<AlignedRead> reads;
    for (size_t i = 0; i < sz; i++) {
        reads.emplace_back(AlignedRead::Load(is, index));
    }
    return reads;
}

//    Moves dir from the record of its current path to the record of its pending correction: a non-empty path is filed
//    under its start in starts, an empty one (a subread) in reads_inside_vertices. Reads are corrected only while
//    their vertices are alive, so neither path is legacy.
void AlignedReadStorage::updateStart(AlignedReadDirection dir) {
    PathDirection corrected = dir.getCorrected();
    VERIFY(!dir.valid() || !dir.getPath().isLegacy());
    VERIFY(!corrected.valid() || !corrected.getPath().isLegacy());
    std::function<DirectionRecord &(bool, Vertex &)> record = [this](bool empty, Vertex &v) -> DirectionRecord & {
        return empty ? this->getSubstringReadsRecord(v.getId()) : this->getOutgoingReadsRecord(v.getId());
    };
    bool unchanged = dir.valid() && corrected.valid() && dir.empty() == corrected.empty() &&
                     dir.getStart() == corrected.getStart();
    if (!unchanged) {
        if (dir.valid()) {
            Vertex &v = dir.getStart();
            DirectionRecord &old = record(dir.empty(), v);
            v.lock();
            RemoveDirection(old, dir);
            v.unlock();
        }
        if (corrected.valid()) {
            Vertex &v = corrected.getStart();
            DirectionRecord &rec = record(corrected.empty(), v);
            v.lock();
            rec.emplace_back(dir);
            v.unlock();
        }
    }
}

//    Whether some sample segment of v lies within [lo, hi) of v's own coordinates. A kpomer chunk stands for each of its
//    (k+1)-mers, so it has one within [lo, hi) iff their overlap with the chunk is at least k + 1 long.
static bool HasSegmentWithin(Vertex &v, __int64_t lo, __int64_t hi, size_t k) {
    for (CoverageSamples::SampleView seg : v.coverage_info.getSamples()) {
        __int64_t start = seg.start();
        __int64_t finish = seg.finish();
        if (seg.sample.type == CoverageSamples::SampleType::kpomer) {
            VERIFY(k > 0);
            if (std::max(start, lo) + __int64_t(k) + 1 <= std::min(finish, hi))
                return true;
        } else if (start >= lo && finish <= hi) {
            return true;
        }
    }
    return false;
}

//    End of a read path: a k-mer end (DBG) must be needed by the read, any other end must have a sample segment inside
//    the read. edges are the path's edges, F the length the path spells and the read is [L, F - R) of it. Returns false
//    on a violation.
static bool CheckFinishSide(Vertex &finish, const std::vector<Edge *> &edges, __int64_t F, __int64_t L, __int64_t R,
                            size_t k, const char *side, const std::string &name,
                            const std::function<std::ostream &()> &report) {
    if (finish.size() == k) {
        if (edges.empty() || R < __int64_t(edges.back()->truncSize()))
            return true;
        report() << name << ": the read does not need the " << side << " k-mer of its path\n";
    } else if (HasSegmentWithin(finish, L - (F - __int64_t(finish.size())), __int64_t(finish.size()) - R, k)) {
        return true;
    } else {
        report() << name << ": the " << side << " vertex of its path has no sample segment inside the read\n";
    }
    return false;
}

bool AlignedReadStorage::checkRecords(std::ostream &os) {
    size_t errors = 0;
    auto report = [&]() -> std::ostream & {
        errors++;
        return os << "Read storage inconsistency: ";
    };
//    Each direction of a valid read must be filed exactly once, under the vertex it starts at, in the map matching its
//    path type. balance counts expected minus filed entries per (map, vertex id, read, direction).
    enum RecordKind {START_RECORD = 0, SUBREAD_RECORD = 1, LEGACY_RECORD = 2};
    const char *record_names[] = {"starts", "reads_inside_vertices", "legacy_reads"};
    std::map<std::tuple<int, Vertex::id_type, const AlignedRead *, bool>, int> balance;
    for (const AlignedRead &read : reads) {
        const GraphPath &path = read.getPath();
        if (!path.valid())
            continue;
        int kind = path.isLegacy() ? LEGACY_RECORD : path.empty() ? SUBREAD_RECORD : START_RECORD;
        balance[{kind, path.getStartId().innerId(), &read, false}]++;
        balance[{kind, path.RC().getStartId().innerId(), &read, true}]++;
    }
    auto scan = [&](libcuckoo::cuckoohash_map<ConstVertexId, std::unique_ptr<DirectionRecord>> &map, int kind) {
        auto lt = map.lock_table();
        for (auto &rec : lt)
            for (const AlignedReadDirection &dir : *rec.second)
                balance[{kind, rec.first.innerId(), &dir.getRead(), dir.isRC()}]--;
    };
    scan(starts, START_RECORD);
    scan(reads_inside_vertices, SUBREAD_RECORD);
    scan(legacy_reads, LEGACY_RECORD);
    for (const auto &entry : balance) {
        if (entry.second == 0)
            continue;
        report() << "read " << std::get<2>(entry.first)->getId()
                 << (std::get<3>(entry.first) ? " backward" : " forward") << " direction is "
                 << (entry.second > 0 ? "missing from" : "stale or duplicated in") << " the "
                 << record_names[std::get<0>(entry.first)] << " record of vertex " << std::get<1>(entry.first) << "\n";
    }

    for (const AlignedRead &read : reads) {
        const GraphPath &path = read.getPath();
        if (!path.valid() || path.isLegacy())
            continue;
        const std::string name = "read " + read.getId();
        if (!starts.contains(path.getStartId()) || !starts.contains(path.RC().getStartId())) {
            report() << name << ": path starts or ends at a deleted vertex\n";
            continue;
        }
//        Follow the path the way PathPosition does: at each vertex take the unique outgoing edge that is a suffix edge or
//        whose code starts with the next fsplits character, consuming its code from the front of fsplits and its rc's
//        code from the back of rsplits.
        const NuclDeck &fs = path.getFSplits();
        const NuclDeck &rs = path.getRSplits();
        Vertex *v = &path.getStart();
        std::vector<Edge *> edges;
        size_t fi = 0;
        size_t ri = rs.size();
        size_t full = v->size();
        bool broken = false;
        while (fi < fs.size() || ri > 0) {
            Edge *next = nullptr;
            size_t matches = 0;
            for (Edge &e : *v) {
                if (e.isCodeEmpty() || (fi < fs.size() && e.getCode()[0] == fs[fi])) {
                    next = &e;
                    matches++;
                }
            }
            if (matches != 1) {
                report() << name << ": path can not be followed from vertex " << v->getInnerId() << " ("
                         << matches << " matching outgoing edges)\n";
                broken = true;
                break;
            }
            const Sequence &code = next->getCode();
            const Sequence &rc_code = next->rc().getCode();
            bool match = (!code.empty() || !rc_code.empty()) && fi + code.size() <= fs.size() && rc_code.size() <= ri;
            for (size_t i = 0; match && i < code.size(); i++)
                match = fs[fi + i] == code[i];
            for (size_t i = 0; match && i < rc_code.size(); i++)
                match = rs[ri - rc_code.size() + i] == rc_code[i];
            if (!match) {
                report() << name << ": path codes do not match edge " << next->getInnerId() << "\n";
                broken = true;
                break;
            }
            fi += code.size();
            ri -= rc_code.size();
            full += next->truncSize();
            edges.push_back(next);
            v = &next->getFinish();
        }
        if (broken)
            continue;
        if (v != &path.getFinish()) {
            report() << name << ": path codes lead to vertex " << v->getInnerId() << " instead of its finish "
                     << path.getFinish().getInnerId() << "\n";
            continue;
        }
        if (path.leftCut() + path.rightCut() >= full) {
            report() << name << ": cuts " << path.leftCut() << " and " << path.rightCut() << " leave nothing of length "
                     << full << "\n";
            continue;
        }
//        Read path normal form: both ends (the start as the finish of the RC), then the path itself. A path may have
//        extra vertices without sample segments at either end (e.g. supre-vertices made by ExtendAdjacent during a
//        resolution); they are skipped, and the rest of the path must be in normal form.
        std::function<std::ostream &()> report_fn = report;
        std::vector<Edge *> rc_edges;
        for (size_t i = edges.size(); i > 0; i--)
            rc_edges.push_back(&edges[i - 1]->rc());
        bool extra_finish = path.getFinish().size() != k && path.getFinish().coverage_info.empty();
        bool extra_start = path.getStart().size() != k && path.getStart().coverage_info.empty();
        bool ends_ok = extra_finish || CheckFinishSide(path.getFinish(), edges, full, path.leftCut(), path.rightCut(), k,
                                                       "finish", name, report_fn);
        ends_ok = (extra_start || CheckFinishSide(path.getStart().rc(), rc_edges, full, path.rightCut(), path.leftCut(),
                                                  k, "start", name, report_fn)) && ends_ok;
        GraphPath inner = path;
        if (extra_finish) {
            inner.setCutRight(0);
            while (!inner.empty() && inner.getFinish().coverage_info.empty())
                inner.pop_back();
        }
        if (extra_start) {
            inner.setCutLeft(0);
            while (!inner.empty() && inner.getStart().coverage_info.empty())
                inner.pop_front();
        }
        bool inner_segmented = !inner.getStart().coverage_info.empty() && !inner.getFinish().coverage_info.empty();
        if (ends_ok && (inner_segmented || (!extra_start && !extra_finish)) && Normalize(inner) != inner)
            report() << name << ": path " << path.str() << " is not in normal form " << Normalize(path).str() << "\n";
    }
    if (errors > 0)
        os << "Read storage inconsistency: " << errors << " violation(s) found" << std::endl;
    return errors == 0;
}

bool AlignedReadStorage::checkConsistency(std::ostream &os) {
    bool ok = checkRecords(os);
    bool listeners_ok = this->fireCheckConsistency(os);
    return ok && listeners_ok;
}

AlignedReadStorage::AlignedReadStorage(logging::Logger &logger, size_t threads,
                                       AssemblyGraph &graph, std::vector<AlignedRead> reads)
        : reads(std::move(reads)), k(graph.getK()) {
    omp_set_num_threads(int(threads));
//    this-> matters: the parameter `reads` was moved from and is empty.
#pragma omp parallel for default(none) schedule(dynamic, 100)
    for(size_t i = 0; i < this->reads.size(); i++) {
        this->reads[i].setPath(Normalize(this->reads[i].getPath()));
        this->fireAddRead(this->reads[i]);
    }
    maintenance = new ag::AlignedReadStorageMaintenance(graph, *this);
}

AlignedReadStorage::AlignedReadStorage(AlignedReadStorage &&other) noexcept {
    *this = std::move(other);
}

AlignedReadStorage &AlignedReadStorage::operator=(AlignedReadStorage &&other) noexcept {
    AlignedReadStorageFire::operator=(std::move(other));
    std::swap(reads, other.reads);
    std::swap(k, other.k);
    std::swap(starts, other.starts);
    std::swap(reads_inside_vertices, other.reads_inside_vertices);
    std::swap(legacy_reads, other.legacy_reads);
    std::swap(maintenance, other.maintenance);
    std::swap(read_correction, other.read_correction);
    if(maintenance != nullptr)
        maintenance->storage = this;
    if(other.maintenance != nullptr)
        other.maintenance->storage = &other;
    return *this;
}

AlignedReadStorage::AlignedReadStorage(AssemblyGraph &graph, std::vector<AlignedRead> reads)
        : reads(std::move(reads)), k(graph.getK()) {
//    this-> matters: the parameter `reads` was moved from and is empty.
    for(size_t i = 0; i < this->reads.size(); i++) {
        this->reads[i].setPath(Normalize(this->reads[i].getPath()));
        this->fireAddRead(this->reads[i]);
    }
    maintenance = new ag::AlignedReadStorageMaintenance(graph, *this);
}

void AlignedReadStorage::Save(const std::experimental::filesystem::path &path) {
    std::ofstream os;
    os.open(path);
    Save(os);
    os.close();
}

AlignedReadStorage AlignedReadStorage::Load(logging::Logger &logger, size_t threads,
                                            const std::experimental::filesystem::path &path,
                                            AssemblyGraph &graph) {
    std::ifstream is;
    is.open(path);
    AlignedReadStorage result(Load(logger, threads, is, graph));
    is.close();
    return std::move(result);
}

AlignedReadStorage::~AlignedReadStorage() {delete maintenance;}

const std::vector<AlignedReadDirection> & AlignedReadStorage::getOutgoingReadsLockFree(const Vertex &vertex) const {
    return getOutgoingReadsRecord(vertex.getId());
}

const std::vector<AlignedReadDirection> & AlignedReadStorage::getOutgoingReads(const Vertex &vertex) const {
    return getOutgoingReadsLockFree(vertex);
}

std::vector<AlignedReadDirection> & AlignedReadStorage::getOutgoingReadsLockFree(const Vertex &vertex) {
    return getOutgoingReadsRecord(vertex.getId());
}

std::vector<AlignedReadDirection> & AlignedReadStorage::getOutgoingReads(const Vertex &vertex) {
    return getOutgoingReadsLockFree(vertex);
}

//AITODO: simplify this outdated system of calls for vertex records
std::vector<AlignedReadDirection> & AlignedReadStorage::getSubstringReadsLockFree(VertexId vertex) {
    return getSubstringReadsRecord(vertex);
}

const std::vector<AlignedReadDirection> & AlignedReadStorage::getSubstringReads(VertexId vertex) const {
    return getSubstringReadsLockFree(vertex);
}

const std::vector<AlignedReadDirection> & AlignedReadStorage::getSubstringReadsLockFree(VertexId vertex) const {
    return getSubstringReadsRecord(vertex);
}

std::vector<AlignedReadDirection> & AlignedReadStorage::getSubstringReads(VertexId vertex) {
    return getSubstringReadsLockFree(vertex);
}
