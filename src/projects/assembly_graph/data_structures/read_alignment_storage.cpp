#include "read_alignment_storage.hpp"
#include "assembly_graph/assembly_graph.hpp"
using namespace ag;

void AlignedReadStorageMaintenance::migrateSubstringReads(Vertex &from, Vertex &to, size_t shift) {
    std::vector<AlignedReadDirection> &src = storage->getSubstringReadsLockFree(from.getId());
    std::vector<AlignedReadDirection> &dst = storage->getSubstringReadsLockFree(to.getId());
    for (AlignedReadDirection &dir : src) {
        dir.getRead().lock();
        size_t len = dir.getPath().len();
        size_t new_left = shift + dir.leftCut();
        dir.setPath(GraphPath(to, new_left, to.size() - new_left - len));
        dir.getRead().unlock();
        dst.emplace_back(dir);
    }
    src.clear();
}

void AlignedReadStorageMaintenance::fireResolveVertex(Vertex &core, const VertexResolutionResult  &resolution) {
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
//            The read collapses to a substring of dest. forward and backward share the same path, so this
//            pop already emptied the backward view too -- but the backward direction is still registered in
//            starts[dest.rc()], a bucket that neither this fire nor the core.rc() fire iterates, so move it
//            to the substring index explicitly (operator== matches by read, so find locates the stale copy).
            AlignedReadDirection bwd = read.backward();
            DestRecs &rc_rec = dest_recs.at(dest.rc().getId());
            rc_rec.first->erase(std::find(rc_rec.first->begin(), rc_rec.first->end(), bwd));
            rc_rec.second->emplace_back(bwd);
            dest_recs.at(dest.getId()).second->emplace_back(dir);
        } else {
            dest_recs.at(dest.getId()).first->emplace_back(dir);
        }
        read.unlock();
    }
}

void AlignedReadStorageMaintenance::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
// Nothing to do here. All read paths stay as they are.

    // std::vector<AlignedReadDirection> &recs = storage->getOutgoingReads(e);
    // std::vector<AlignedReadDirection> &new_recs = storage->getOutgoingReads(v.front());
    // std::vector<AlignedReadDirection> &subread_recs = storage->getSubstringReads(v.getId());
    // for(AlignedReadDirection &dir : recs) {
    //     if(dir.getStart() == v) {
    //         //In case this path consists only of e, it could have already been processed by rc call
    //         VERIFY(dir.empty());
    //         subread_recs.emplace_back(dir);
    //     } else if(dir.isSingleton()) {
    //         dir.setPath(GraphPath(v, dir.leftCut(), dir.rightCut()));
    //         subread_recs.emplace_back(dir);
    //     } else {
    //         dir.pop_front(v.incFront());
    //         new_recs.emplace_back(dir);
    //     }
    // }
}


void AlignedReadStorageMaintenance::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    size_t left_skip = 0;
    GraphPath prefix(path.frontEdge().getStart());
    Vertex & start = path.getStart();
    std::vector<AlignedReadDirection> &collected_recs = storage->getOutgoingReads(start);
    for(Edge &edge : path.edges()) {
        if (edge == path.backEdge())
            continue;
        prefix += edge;
        left_skip += edge.rc().truncSeq().size();
        for(const AlignedReadDirection &dir : storage->getOutgoingReads(edge.getFinish())) {
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

void AlignedReadStorageMaintenance::fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
    VERIFY(new_vertex != path.getStart());//All-suffix and all-prefix paths should be handled by mergePathToEdge
    VERIFY(new_vertex != path.getFinish());
//    path.getStart()/getFinish() survive this merge (only inner vertices are isolated below) and keep
//    their own reads -- they are not "the same entity" as new_vertex just because one side is
//    unbranching, same as SampleSegmentTracker no longer folds their samples into new_vertex either.
    Edge &forwardEdge = new_vertex.front();//Last edge is now a suffix edge from new_vertex to end of the path
    Edge &backwardEdge = new_vertex.incFront();
    size_t left_skip = 0;
    size_t right_skip = new_vertex.size() - path.getStart().size();
    std::vector<AlignedReadDirection> &recs = storage->getOutgoingReads(new_vertex);
    std::vector<AlignedReadDirection> &vrecs = storage->getSubstringReads(new_vertex.getId());
    GraphPath suffix(path);
    GraphPath prefix(path.frontEdge().getStart());
    GraphPath out_edge(new_vertex.front());
    for(Vertex &path_vertex : path.vertices()) {
        if (path_vertex != new_vertex && path_vertex != new_vertex.rc()) {
            for(AlignedReadDirection dir : storage->getOutgoingReads(path_vertex)) {
                dir.getRead().lock();
                if(dir.getStart() == new_vertex) {
                    // Here we check the case that this read has already been handled by rc call
                    VERIFY(dir.getFinish() == new_vertex);
                    vrecs.emplace_back(dir);
                } else {
                    size_t cut_left = left_skip + dir.leftCut();
                    if (dir.getFSplits().size() + prefix.getFSplits().size() > backwardEdge.getCode().size()) {
                        //When the read path goes beyond merged path, we
                        dir.setCutLeft(0);
                        dir.replace_back(suffix, out_edge);
                        dir.setCutLeft(cut_left);
                        VERIFY(dir.getStart() == new_vertex);
                        recs.emplace_back(dir);
                    } else {
                        size_t len = dir.getPath().len();
                        dir.setPath(GraphPath(new_vertex, cut_left, new_vertex.size() - cut_left - len));
                        vrecs.emplace_back(dir);
                    }
                }
                dir.getRead().unlock();
            }
            for(AlignedReadDirection dir : storage->getSubstringReads(path_vertex.getId())) {
                if(dir.getStart() == new_vertex) {
                    // Here we check the case that this read has already been handled by rc call
                    VERIFY(dir.getFinish() == new_vertex);
                    vrecs.emplace_back(dir);
                } else {
                    size_t cut_left = left_skip + dir.leftCut();
                    size_t len = dir.getPath().len();
                    dir.setPath(GraphPath(new_vertex, cut_left, new_vertex.size() - cut_left - len));
                    vrecs.emplace_back(dir);
                }
                dir.getRead().unlock();
            }
        }
        left_skip += suffix.frontEdge().rc().truncSize();
        right_skip -= suffix.frontEdge().truncSize();
        prefix += suffix.frontEdge();
        suffix.pop_front();
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
    storage->reads_inside_vertices.insert(vertex.getId(), std::make_unique<std::vector<AlignedReadDirection>>());
    storage->starts.insert(vertex.getId(), std::make_unique<std::vector<AlignedReadDirection>>());
}

void AlignedReadStorageMaintenance::fireDeleteVertex(Vertex &vertex) {
    std::vector<AlignedReadDirection> tmp = std::move(storage->getSubstringReadsLockFree(vertex.getId()));
    storage->reads_inside_vertices.erase(vertex.getId());
    storage->reads_inside_vertices.insert(vertex.getId().legacyId(),
                                           std::make_unique<std::vector<AlignedReadDirection>>(std::move(tmp)));
    std::vector<AlignedReadDirection> &recs = storage->getSubstringReadsLockFree(vertex.getId().legacyId());
    for (AlignedReadDirection &dir: recs) {
        GraphPath path = dir.getPath();
        dir.setPath(GraphPath::LegacyPath(vertex.getId(), vertex.rc().getId(), path.leftCut(), path.rightCut()));
    }
    storage->starts.erase(vertex.getId());
}

void AlignedReadStorageMaintenance::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
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
    VERIFY(!left.isSuffix());
    VERIFY(!right.isPrefix());
    VERIFY(left.getFinish().outDeg() == 0 && left.getFinish().inDeg() == 1);
    VERIFY(right.getStart().outDeg() == 1 && right.getStart().outDeg() == 1);
    VERIFY(left_al.targetLength() == right_al.targetLength());
    size_t left_skip = left.fullSize() - left_al.queryLength();
    size_t right_skip = right.fullSize() - right_al.queryLength();
    std::vector<AlignedReadDirection> &new_rec = storage->getOutgoingReadsRecord(new_edge.getStart().getId());
    std::vector<AlignedReadDirection> &new_rec_rc = storage->getOutgoingReadsRecord(new_edge.rc().getStart().getId());
//        Reads on left edge will be handled by rc
    for(AlignedReadDirection dir : storage->getOutgoingReadsRecord(right.getStart().getId())) {
        dir.getRead().lock();
        size_t new_left = dir.leftCut() >= right_al.queryLength() ?
                          left_skip + right_al.targetLength() + dir.leftCut() - right_al.queryLength() :
                          left_skip + right_al.lastColumnByQpos(dir.leftCut()).getTpos();
        if(dir.isSingleton()) {
            size_t old_right = right.fullSize() - dir.rightCut();
            size_t new_right = old_right >= right_al.queryLength() ?
                               left_skip + right_al.targetLength() + old_right - right_al.queryLength() :
                               left_skip + right_al.firstColumnByQpos(old_right).getTpos();
            dir.setCutRight(new_edge.fullSize() - new_right);
            new_rec_rc.emplace_back(dir.RC());
        }
        dir.setCutLeft(0);
        dir.forcePushFront(left);// While this operation is underway this path is disconnected.
        dir.setCutLeft(new_left);
        dir.getRead().unlock();
        new_rec.emplace_back(dir);
    }
}

void AlignedReadStorageMaintenance::fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph & graph) {
#pragma omp parallel for default(none) schedule(dynamic, 100)
    for(size_t i = 0; i < storage->size(); i++) {
        AlignedRead &read = (*storage)[i];
        read.resetEdgeCodes();
    }
}

bool AlignedReadStorageMaintenance::fireCheckConsistency() {
    size_t cnt_starts = 0;
//        fireCheckConsistency is only ever called single-threaded, never concurrently with graph
//        modifications; lock_table() here is just libcuckoo's only full-table iteration API, not
//        synchronization against any other running thread.
    {
        auto lt = storage->starts.lock_table();
        for (auto &rec : lt) {
            ConstVertexId vid = rec.first;
            std::vector<AlignedReadDirection> &dirs = *rec.second;
            for(AlignedReadDirection dir : dirs) {
                VERIFY(!dir.empty() && dir.getStart() == *vid);
            }
            cnt_starts += dirs.size();
        }
    }
    size_t cnt_reads = 0;
    for(AlignedRead & read: storage->reads) {
        if(read.getPath().empty())
            continue;
        cnt_reads++;
    }
    VERIFY(cnt_starts == cnt_reads * 2);
    return true;
}

void AlignedReadStorageMaintenance::fireMergeLoop(const GraphPath &path, Vertex &new_vertex) {
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
    read.delayedInvalidate();
    this->fireDelayedInvalidateRead(read, message);
}

void AlignedReadStorage::rerouteRead(AlignedRead &alignedRead, GraphPath corrected,
                                     const string &message) {
//    VERIFY(corrected.truncLen() >= 500);
    Normalize(corrected);
    alignedRead.correct(std::move(corrected));
    this->fireDelayedRerouteRead(alignedRead, message);
}

bool AlignedReadStorage::apply(AlignedRead &alignedRead) {
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

void AlignedReadStorage::updateStart(AlignedReadDirection dir) {
    if (dir.valid() && dir.getCorrected().valid() && !dir.empty() && !dir.getCorrected().empty() && dir.getStart() == dir.getCorrected().getStart())
        return;
    if (dir.valid() && dir.getCorrected().valid() && dir.empty() && dir.getCorrected().empty() && dir.getStart() == dir.getCorrected().getStart()) {
        return;
    }
    if(dir.valid()) {
        if(dir.empty()) {
            std::vector<AlignedReadDirection> &old = this->getSubstringReads(dir.getStart().getId());
            dir.getStart().lock();
            // TODO: Check performance. This is deletion from vector.
            old.erase(std::find(old.begin(), old.end(), dir));
            dir.getStart().unlock();
        } else {
            std::vector<AlignedReadDirection> &old = this->getOutgoingReadsRecord(dir.getStart().getId());
            dir.getStart().lock();
            // TODO: Check performance. This is deletion from vector.
            old.erase(std::find(old.begin(), old.end(), dir));
            dir.getStart().unlock();
        }
    }
    if(dir.getCorrected().valid()) {
        if (dir.empty()) {
            std::vector<AlignedReadDirection> &old = this->getSubstringReads(dir.getStart().getId());
            dir.getStart().lock();
            // TODO: Check performance. This is deletion from vector.
            old.emplace_back(dir);
            dir.getStart().unlock();
        } else {
            auto &rec = this->getOutgoingReadsRecord(dir.getCorrected().getStart().getId());
            Vertex &v = dir.getCorrected().getStart();
            v.lock();
            rec.emplace_back(dir);
            v.unlock();
        }
    }
}

bool AlignedReadStorage::checkConsistency() {
    for(AlignedRead &al : reads) {
        if(!al.valid())
            continue;
        std::vector<AlignedReadDirection> &start = getOutgoingReadsRecord(al.getPath().getStart().getId());
        bool f1 = false;
        for(AlignedReadDirection &dir : start) {
            if(dir.getRead().getId() == al.getId()) {
                f1 = true;
                break;
            }
        }
        std::vector<AlignedReadDirection> &rcstart = getOutgoingReadsRecord(al.getPath().getFinish().rc().getId());
        bool f2 = false;
        for(AlignedReadDirection &dir : rcstart) {
            if(dir.getRead().getId() == al.getId()) {
                f2 = true;
                break;
            }
        }
        if(!f1 || !f2) {
            VERIFY_MSG(false, al.getId() << " " << al.getPath().getStart().getId() << " " << al.getPath().getFinish().getId() << " " << al.getPath().getFSplits().str() << " " << al.getPath().getRSplits().str());
            return false;
        }
    }
    return this->fireCheckConsistency();
}

AlignedReadStorage::AlignedReadStorage(logging::Logger &logger, size_t threads,
                                       AssemblyGraph &graph, std::vector<AlignedRead> reads)
        : reads(std::move(reads)) {
    omp_set_num_threads(int(threads));
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(reads)
    for(size_t i = 0; i < reads.size(); i++) {
        this->fireAddRead(reads[i]);
    }
    maintenance = new ag::AlignedReadStorageMaintenance(graph, *this);
}

AlignedReadStorage::AlignedReadStorage(AlignedReadStorage &&other) noexcept {
    *this = std::move(other);
}

AlignedReadStorage &AlignedReadStorage::operator=(AlignedReadStorage &&other) noexcept {
    AlignedReadStorageFire::operator=(std::move(other));
    std::swap(reads, other.reads);
    std::swap(starts, other.starts);
    std::swap(reads_inside_vertices, other.reads_inside_vertices);
    std::swap(maintenance, other.maintenance);
    if(maintenance != nullptr)
        maintenance->storage = this;
    if(other.maintenance != nullptr)
        other.maintenance->storage = &other;
    return *this;
}

AlignedReadStorage::AlignedReadStorage(AssemblyGraph &graph, std::vector<AlignedRead> reads)
        : reads(std::move(reads)) {
    for(size_t i = 0; i < reads.size(); i++) {
        Normalize(reads[i].getPath());
        this->fireAddRead(reads[i]);
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
