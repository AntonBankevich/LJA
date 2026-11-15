#include "suffix_tracker.hpp"
using namespace ag;

size_t SuffixRecord::getNumberOfPaths() const {
    size_t res = 0;
    for(const auto &p : paths)
        res += p.second;
    return res;
}

size_t SuffixRecord::countStartsWith(const GraphPath &path) const {
    VERIFY(path.empty() || path.getStart() == eid->getFinish());
    if(path.empty()) return countStartsWith(Sequence());
    return countStartsWith(Sequence(path.getFSplits().begin(), path.getFSplits().end() - path.backEdge().getCode().size() + 1));
}

bool SuffixRecord::empty() const {
    if(paths.empty())
        return true;
    for(const auto &path: paths) {
        if(path.second > 0)
            return false;
    }
    return true;
}

size_t SuffixRecord::countStartsWith(const Sequence &seq) const {
//        lock();
    int cnt = 0;
    for (const std::pair<Sequence, int> &rec: paths) {
        if (rec.first.startsWith(seq)) {
            cnt += rec.second;
        }
    }
//        unlock();
    VERIFY(cnt >= 0);
    return cnt;
}

std::string SuffixRecord::str() const {
    std::stringstream ss;
    for (const auto &path: paths) {
        ss << path.first << " " << path.second << std::endl;
    }
    return ss.str();
}

void SuffixRecord::updateZero(size_t old_val, size_t new_val) {
    if(new_val == 0)
        zero_cnt++;
    if(old_val == 0)
        zero_cnt--;
}

void SuffixRecord::removeZero() {
    std::vector<std::pair<Sequence, int>> new_paths;
    for (std::pair<Sequence, int> &rec: this->paths) {
        if (rec.second != 0) {
            new_paths.emplace_back(std::move(rec.first), rec.second);
        }
    }
    std::swap(this->paths, new_paths);
    this->zero_cnt = 0;
}

void SuffixRecord::lockFreeChangePathCnt(const Sequence &min_seq, const Sequence &max_seq, int diff) {
    if(min_seq.empty()) {
        return;
    }
    if (diff == 0) return;
    if(diff > 0) {
        addPath(min_seq, diff);
        return;
    }
    for (auto it = paths.begin(); diff != 0 && it != paths.end(); ++it) {
        if((it->first == min_seq || it->first == max_seq) && it->second + diff >= 0) {
            it->second += diff;
            updateZero(it->second - diff, it->second);
            diff = 0;
        }
    }
    for (auto it = paths.begin(); diff != 0 && it != paths.end(); ++it) {
        if(it->second != 0 && max_seq.startsWith(it->first) && it->first.startsWith(min_seq)) {
            VERIFY(it->second > 0);
            int taken = std::min(it->second, -diff);
            it->second -= taken;
            updateZero(it->second + taken, it->second);
            diff += taken;
        }
    }
    VERIFY_MSG(diff == 0, "Attempting to remove path that is not present in the record.");
//    Entries dropping to zero are kept for a while (addPath revives one if the same suffix comes back); clean up once
//    they make up a third of the record.
    if (zero_cnt > paths.size() / 3) {
        removeZero();
    }
}

void SuffixRecord::changePathCnt(const Sequence &min_seq, const Sequence &max_seq, int diff) {
    lock();
    lockFreeChangePathCnt(min_seq, max_seq, diff);
    unlock();
}

void SuffixRecord::addPath(const Sequence &seq, int diff) {
    VERIFY(diff > 0);
    VERIFY(!seq.empty());
    for (auto & path : paths) {
        if(seq == path.first) {
            updateZero(path.second, 1);
            path.second += diff;
            return;
        }
    }
    paths.emplace_back(seq, diff);
}

void SuffixRecord::directAddPath(const Sequence &seq, size_t cnt) {
    if (!seq.empty()) {
        paths.emplace_back(seq, cnt);
    }
}

void SuffixRecord::clear() {
    paths.clear();
}

void SuffixRecord::resetCodes(Vertex &start) {
    std::vector<std::pair<Sequence, int>> old = std::move(paths);
    for(const std::pair<Sequence, int> &rec : old) {
        if(rec.second == 0)
            continue;
        GraphPath path(start, rec.first);
        path.resetEdgeCodes();
        changePathCnt(Sequence(path.getFSplits()), Sequence(path.getFSplits()), rec.second);
    }
}

std::function<std::string(const Edge & )> SuffixTracker::labeler() const {
    return [this](const Edge &edge) {
        if (!hasSuffixRecord(edge))
            return std::string();
        const SuffixRecord &rec = getSuffixRecord(edge);
        std::stringstream ss;
        size_t cnt = 0;
        for (const auto &ext: rec) {
            if (ext.second > 0) {
                if (cnt < 30)
                    ss << "\\n" << ext.first << "(" << ext.second << ")";
                cnt++;
            }
        }
        if (cnt > 30) {
            ss << "\\n" << "and another " << (cnt - 30) << " records";
        }
        return ss.str();
    };
}

void SuffixTracker::processPath(PathPosition left, PathPosition right, int diff) {
    while (right != left && right.prevEdge().getCode().empty()) {
        --right;
    }
    //A path with no edge with a code (e.g. a read that is just res->end after a merge) visits no edge
    //that has an edge with a code after it, so it contributes to no suffix record.
    if (right == left)
        return;
    PathPosition from_pos = left;
    PathPosition to_pos = from_pos + 1;
    size_t clen = left.nextEdge().truncSize();
    VERIFY(to_pos <= right);
    Sequence read_seq(left.getFPos(), right.getFPos());
    while (from_pos != right) {
        Edge &edge = from_pos.nextEdge();
        clen -= edge.truncSize();
        ++from_pos;
        SuffixRecord *erec = findRecord(edge);
        const size_t need = erec == nullptr ? 0 : erec->getMaxSuffixLength();
//        An edge without a record stores nothing, but to_pos must stay ahead of from_pos for the edges after it. A
//        stored suffix never ends with an edge without a code.
        while (to_pos < right && (clen < need || (erec == nullptr && to_pos == from_pos) ||
                                  (to_pos != from_pos && to_pos.prevEdge().getCode().empty()))) {
            Edge &new_edge = to_pos.nextEdge();
            clen += new_edge.truncSize();
            to_pos += new_edge;
        }
        VERIFY(to_pos <= right);
        VERIFY(to_pos == from_pos || !to_pos.prevEdge().getCode().empty());
        if (erec != nullptr && clen >= min_len) {
//            Trailing virtually suffix edges (and edges without a code) do not count for the canonical form: an operation
//            may give a read a coded virtually suffix edge in place of a code-less one at its end (flip), so the form
//            must not depend on them. The canonical form then ends with the first letter of the last remaining edge.
            PathPosition min_pos = to_pos;
            while (min_pos != from_pos && (min_pos.prevEdge().isVirtualSuffix() || min_pos.prevEdge().getCode().empty()))
                --min_pos;
            if (min_pos == from_pos) {
                erec->changePathCnt(Sequence(), Sequence(), diff);
            } else {
                Sequence subseq_max_code = read_seq.Subseq(from_pos.getFPos() - left.getFPos(),
                        to_pos.getFPos() - left.getFPos());
                size_t min_size = min_pos.getFPos() - from_pos.getFPos() - min_pos.prevEdge().getCode().size() + 1;
                erec->changePathCnt(subseq_max_code.Subseq(0, min_size), subseq_max_code, diff);
            }
        }
    }
}

void SuffixTracker::fillFromStorage(logging::Logger &logger, size_t threads) {
    logger.info() << "Collecting and storing read suffixes" << std::endl;
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(storage) schedule(dynamic, 100)
    for (size_t i = 0; i < storage->size(); i++) {
        if(!(*storage)[i].getPath().empty())
            fireAddRead((*storage)[i]);
    }
    logger.info() << "Finished collecting and storing read suffixes" << std::endl;
}

//    Legacy paths (subreads of deleted vertices) are empty, so they hold no suffixes; their vertex is gone, so they must
//    not be walked.
void SuffixTracker::fireAddRead(const AlignedRead &read) {
    if(read.valid() && !read.getPath().isLegacy()) {
        addSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
        addSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
    }
}

void SuffixTracker::fireRerouteRead(AlignedRead &read) {
    if(read.valid() && !read.getPath().isLegacy()) {
        removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
        removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
    }
    if(read.getCorrected().valid() && !read.getCorrected().isLegacy()) {
        addSubpath(read.getCorrected().firstPosition(), read.getCorrected().lastPosition());
        addSubpath(read.getCorrected().lastPosition().RC(), read.getCorrected().firstPosition().RC());
    }
}

void SuffixTracker::fireInvalidateRead(AlignedRead &read) {
    if(read.valid() && !read.getPath().isLegacy()) {
        removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
        removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
    }
}

SuffixTracker::SuffixTracker(AlignedReadStorage &storage, AssemblyGraph &graph,
                             size_t _min_len, size_t _max_len) :
        AlignedReadStorageListener(storage, "SuffixTracker"), ResolutionListener(graph, "SuffixTracker"), storage(&storage),
        min_len(_min_len), max_len(_max_len) {
//    Records on every edge that is not virtually prefix (on a DBG, where no edge has a virtual type, on every edge).
//    Operations then give a new edge a record iff the edge it takes over from had one.
    for(Edge &e : graph.edges()) {
        if (!e.isVirtualPrefix())
            edge_data.insert(e.getId(), std::make_unique<SuffixRecord>(e, max_len));
    }
}

void SuffixTracker::moveRecord(const Edge &from, Edge &to, const Edge *through) {
    SuffixRecord *src = findRecord(from);
    if (src == nullptr)
        return;
    size_t max_len = src->getMaxSuffixLength();
    if (through != nullptr)
        max_len -= std::min(max_len, through->truncSize());
    std::unique_ptr<SuffixRecord> rec = std::make_unique<SuffixRecord>(to, max_len);
    rec->paths = std::move(src->paths);
    VERIFY(src->paths.empty());
    rec->zero_cnt = src->zero_cnt;
    src->zero_cnt = 0;
    if (through != nullptr && !through->getCode().empty()) {
        const Sequence &cut = through->getCode();
        size_t kept = 0;
        for (std::pair<Sequence, int> &p : rec->paths)
            if (p.second > 0 && p.first.size() > cut.size() && p.first.startsWith(cut))
                rec->paths[kept++] = {p.first.Subseq(cut.size()), p.second};
        rec->paths.resize(kept);
        rec->zero_cnt = 0;
    }
    edge_data.insert_or_assign(to.getId(), std::move(rec));
}

void SuffixTracker::copyRecord(const Edge &from, Edge &to, const Edge *through) {
    const SuffixRecord *src = findRecord(from);
    if (src == nullptr)
        return;
    size_t max_len = src->getMaxSuffixLength();
    if (through != nullptr)
        max_len -= std::min(max_len, through->truncSize());
    std::unique_ptr<SuffixRecord> rec = std::make_unique<SuffixRecord>(to, max_len);
    if (through == nullptr || through->getCode().empty()) {
        rec->paths = src->paths;
        rec->zero_cnt = src->zero_cnt;
    } else {
        const Sequence &cut = through->getCode();
        for (const std::pair<Sequence, int> &p : *src) {
            if (p.second > 0 && p.first.size() > cut.size() && p.first.startsWith(cut))
                rec->directAddPath(p.first.Subseq(cut.size()), p.second);
        }
    }
    edge_data.insert_or_assign(to.getId(), std::move(rec));
}

//    Every new edge takes over from one old edge, after which reads continue the same way, and gets its record. The
//    source is deleted by the operation, so it is moved, unless the mirrored call needs it too: a path whose last edge
//    is the rc of its first one (forward source path.backEdge(), mirrored source path.frontEdge().rc()).
void SuffixTracker::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    if (path.backEdge() == path.frontEdge().rc())
        copyRecord(path.backEdge(), new_edge);
    else
        moveRecord(path.backEdge(), new_edge);
}

//    Vertex sizes strictly change along the run, so its last edge is never the rc of its first one.
void SuffixTracker::fireMergePrefixPath(const RAGraphPath &path, Edge &new_edge) {
    moveRecord(path.backEdge(), new_edge);
}

void SuffixTracker::fireMergePath(const RAGraphPath &path, Vertex &vertex) {
    if (path.backEdge() == path.frontEdge().rc())
        copyRecord(path.backEdge(), vertex.front());
    else
        moveRecord(path.backEdge(), vertex.front());
}

//    The mirrored call's source is left.rc(), which is right itself if the tips are each other's rc.
void
SuffixTracker::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &,
                                   const AlignmentForm &) {
    if (right == left.rc())
        copyRecord(right, new_edge);
    else
        moveRecord(right, new_edge);
}

//    No read crosses a split vertex (the split contract), so the pieces before the last have empty records. A self-rc
//    edge fires only once.
void SuffixTracker::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
    for (Edge &piece : split.edges())
        if (piece != split.backEdge())
            edge_data.insert_or_assign(piece.getId(), std::make_unique<SuffixRecord>(piece, max_len));
    moveRecord(edge, split.backEdge());
}

//    Forward call: u -> w_i replaces core -> w_i, which is moved unless it is its own rc (then the mirrored call uses it
//    too). Mirrored call: edge = out is the single outgoing edge of the rc side, and w_i.rc() -> u.rc() replaces the
//    path e -> out for each incoming e, so it takes e's record continued through out (each e is used once).
void SuffixTracker::fireFlipVertex(Edge &edge, const std::unordered_map<EdgeId, EdgeId> &extensions) {
    if (edge.isSuffix()) {
        for (Edge &e : edge.getFinish()) {
            if (e == e.rc())
                copyRecord(e, *extensions.at(e.getId()));
            else
                moveRecord(e, *extensions.at(e.getId()));
        }
    } else {
        for (Edge &e : edge.getStart().incoming())
            moveRecord(e, *extensions.at(e.getId()), &edge);
    }
}

void SuffixTracker::fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
    std::function<void(size_t, Edge &)> task = [this](size_t, Edge &e) {
        SuffixRecord *rec = findRecord(e);
        if (rec != nullptr)
            rec->resetCodes(e.getFinish());
    };
    processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, task);
}

bool SuffixTracker::fireCheckConsistency(std::ostream &os) {
//        Single-threaded only; lock_table() is just libcuckoo's only full-table iteration API.
    bool ok = true;
    for(AlignedRead &read : *storage) {
        if(read.valid() && !read.getPath().isLegacy()) {
            removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
            removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
        }
    }
    {
        auto lt = edge_data.lock_table();
        for (auto &p : lt) {
            SuffixRecord &rec = *p.second;
            rec.removeZero();
            if(rec.begin() != rec.end()) {
                os << "SuffixTracker: record of edge " << p.first << " keeps suffixes not given by the stored reads:\n"
                   << rec.str();
                ok = false;
            }
        }
    }
    for(AlignedRead &read : *storage) {
        if(read.valid() && !read.getPath().isLegacy()) {
            addSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
            addSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
        }
    }
    return ok;
}

//e = u -> w becomes u -> v -> w; reads continue after v -> w as they did after e. A self-rc e fires only once.
void SuffixTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    VERIFY_MSG(e.getVirtualType() != VirtualType::unset,
               "A DBG can not be converted to a supregraph with a SuffixTracker attached: attach it after conversion");
    moveRecord(e, v.front());
}

void SuffixTracker::fireDeleteEdge(Edge &edge) {
    edge_data.erase(edge.getId());
}

//The pair in -> out is replaced by an edge into out's finish w (new_vertex -> w or a direct edge), after which reads
//continue as they did after the path in -> out: in's suffixes continuing through out, out's code cut off. in always has
//a record (core has several incoming edges), and w may now have several incoming edges, so every such edge gets one.
void SuffixTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    for (const auto &rec : resolution)
        copyRecord(rec.second.incoming(), rec.first->front(), &rec.second.outgoing());
    for (const auto &rec : resolution.newEdges())
        copyRecord(rec.second.incoming(), *rec.first, &rec.second.outgoing());
}

bool SuffixTracker::checkConsistency(std::ostream &os) {
    size_t errors = 0;
    for (Edge &e : getGraph().edges()) {
        if (!e.isVirtualPrefix() && findRecord(e) == nullptr) {
            os << "Suffix tracker inconsistency: edge " << e.getInnerId() << " " << e.getStart().getInnerId() << " -> "
               << e.getFinish().getInnerId() << " is not virtually prefix but has no record\n";
            errors++;
        }
    }
    return errors == 0;
}
