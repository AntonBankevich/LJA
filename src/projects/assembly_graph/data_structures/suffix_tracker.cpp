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

void SuffixRecord::lockFreeChangePathCnt(const Sequence &seq, int diff) {
    if (seq.empty() || diff == 0)
        return;
    if (diff > 0) {
        addPath(seq, diff);
        return;
    }
    auto it = std::find_if(paths.begin(), paths.end(),
                           [&seq](const std::pair<Sequence, int> &rec) { return rec.first == seq; });
    VERIFY_MSG(it != paths.end() && it->second + diff >= 0,
               "Attempting to remove path that is not present in the record.");
    it->second += diff;
    updateZero(it->second - diff, it->second);
//    Entries dropping to zero are kept for a while (addPath revives one if the same suffix comes back); clean up once
//    they make up a third of the record.
    if (zero_cnt > paths.size() / 3) {
        removeZero();
    }
}

void SuffixRecord::changePathCnt(const Sequence &seq, int diff) {
    lock();
    lockFreeChangePathCnt(seq, diff);
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
    paths.clear();
    zero_cnt = 0;
    for(const std::pair<Sequence, int> &rec : old) {
        if(rec.second == 0)
            continue;
        GraphPath path(start, rec.first);
        path.resetEdgeCodes();
        changePathCnt(Sequence(path.getFSplits()), rec.second);
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

void SuffixTracker::processPath(PathPosition left, PathPosition right, NuclDeck::IteratorInfo cut, int diff) {
    VERIFY_MSG(cut.getCur() <= right.getFPos().getDeck().end().getCur(), "Codes were removed past a read's suffix cut");
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
        const int begin = from_pos.getFPos().getCur();
        const int end = std::min(to_pos.getFPos().getCur(), cut.getCur());
        if (erec != nullptr && clen >= min_len && end > begin) {
            const int offset = left.getFPos().getCur();
            erec->changePathCnt(read_seq.Subseq(begin - offset, end - offset), diff);
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
void SuffixTracker::processPath(const GraphPath &path, NuclDeck::IteratorInfo fend, NuclDeck::IteratorInfo rend,
                                int diff) {
    if (path.valid() && !path.isLegacy()) {
        processPath(path.firstPosition(), path.lastPosition(), fend, diff);
        processPath(path.lastPosition().RC(), path.firstPosition().RC(), rend, diff);
    }
}

void SuffixTracker::fireAddRead(const AlignedRead &read) {
    processRead(read, 1);
}

//    The corrected path becomes the read's path as it is now (no graph operation can come in between), and its suffix
//    ends then are its current ends (AlignedRead::setPath).
void SuffixTracker::fireRerouteRead(AlignedRead &read) {
    processRead(read, -1);
    const GraphPath &corrected = read.getCorrected();
    processPath(corrected, corrected.getFSplits().end(), corrected.getRSplits().end(), 1);
}

void SuffixTracker::fireInvalidateRead(AlignedRead &read) {
    processRead(read, -1);
}

SuffixTracker::SuffixTracker(AlignedReadStorage &storage, AssemblyGraph &graph,
                             size_t _min_len, size_t _max_len) :
        AlignedReadStorageListener(storage, "SuffixTracker"), ResolutionListener(graph, "SuffixTracker"), storage(&storage),
        min_len(_min_len), max_len(_max_len) {
    for(Edge &e : graph.edges()) {
        if (holdsRecord(e))
            edge_data.insert(e.getId(), std::make_unique<SuffixRecord>(e, max_len));
    }
}

void SuffixTracker::transferRecord(const Edge &from, Edge &to, bool copy, const Sequence &cut, size_t cut_len) {
    SuffixRecord &src = getSuffixRecord(from);
    size_t max_len = src.getMaxSuffixLength();
    max_len -= std::min(max_len, cut_len);
    std::unique_ptr<SuffixRecord> rec = std::make_unique<SuffixRecord>(to, max_len);
    if (cut.empty()) {
        rec->paths = copy ? src.paths : std::move(src.paths);
        rec->zero_cnt = src.zero_cnt;
    } else {
        for (const std::pair<Sequence, int> &p : src)
            if (p.second > 0 && p.first.size() > cut.size() && p.first.startsWith(cut))
                rec->directAddPath(p.first.Subseq(cut.size()), p.second);
    }
    edge_data.insert_or_assign(to.getId(), std::move(rec));
}

void SuffixTracker::moveRecord(const Edge &from, Edge &to, const Edge *through) {
    transferRecord(from, to, false, through == nullptr ? Sequence() : through->getCode(),
                   through == nullptr ? 0 : through->truncSize());
}

void SuffixTracker::copyRecord(const Edge &from, Edge &to, const Edge *through) {
    transferRecord(from, to, true, through == nullptr ? Sequence() : through->getCode(),
                   through == nullptr ? 0 : through->truncSize());
}

//    The new edge takes over from the path's last edge, after which reads continue the same way.
void SuffixTracker::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    moveRecord(path.backEdge(), new_edge);
}

//    Reads continue after Q's last edge as after P's last edge, and after a replacement as after the edge it replaces,
//    except that P's codes are cut off when they moved into the replacements of the edges into P's start.
void SuffixTracker::fireReplacePath(const PathReplacement &rep) {
    const RAGraphPath &old_path = rep.oldPath();
    if (rep.newPath().size() > 0 && holdsRecord(rep.newPath().backEdge()))
        moveRecord(old_path.backEdge(), rep.newPath().backEdge());
    if (rep.startReplaced()) {
        const bool cut = rep.codesToReplacedStart();
        const Sequence codes = cut ? rep.oldCodes() : Sequence();
        const size_t cut_len = cut ? old_path.truncLen() : 0;
        for (Edge &in : old_path.getStart().incoming())
            if (holdsRecord(rep.replacement(in)))
                transferRecord(in, rep.replacement(in), false, codes, cut_len);
    }
    if (rep.endReplaced())
        for (Edge &out : old_path.getFinish())
            if (!rep.closesLoop(out) && holdsRecord(rep.replacement(out)))
                moveRecord(out, rep.replacement(out));
}

void
SuffixTracker::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &,
                                   const AlignmentForm &) {
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
    for(AlignedRead &read : *storage)
        processRead(read, -1);
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
    for(AlignedRead &read : *storage)
        processRead(read, 1);
    return ok;
}

//e = u -> w becomes u -> v -> w; reads continue after v -> w as they did after e. A self-rc e fires only once.
void SuffixTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    VERIFY_MSG(e.getVirtualType() != VirtualType::unset,
               "A DBG can not be converted to a supregraph with a SuffixTracker attached: attach it after conversion");
    if (holdsRecord(v.front()))
        moveRecord(e, v.front());
}

void SuffixTracker::fireDeleteEdge(Edge &edge) {
    edge_data.erase(edge.getId());
}

//The pair in -> out is replaced by an edge into out's finish w (new_vertex -> w or a direct edge), after which reads
//continue as they did after the path in -> out: in's suffixes continuing through out, out's code cut off. in always has
//a record (it is virtually suffix: core has several incoming edges); new_vertex -> w always needs one, a direct edge
//only if virtually suffix. A copy of core (1-to-1 pair) stands for core itself: u -> copy, virtually suffix, takes in's
//whole record (all of in's reads continue through out); copy -> w is a prefix edge and holds none.
void SuffixTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    for (const auto &rec : resolution) {
        if (rec.first->size() == core.size())
            copyRecord(rec.second.incoming(), rec.first->incFront());
        else
            copyRecord(rec.second.incoming(), rec.first->front(), &rec.second.outgoing());
    }
    for (const auto &rec : resolution.newEdges())
        if (holdsRecord(*rec.first))
            copyRecord(rec.second.incoming(), *rec.first, &rec.second.outgoing());
}

bool SuffixTracker::checkConsistency(std::ostream &os) {
    size_t errors = 0;
    for (Edge &e : getGraph().edges()) {
        if (holdsRecord(e) && findRecord(e) == nullptr) {
            os << "Suffix tracker inconsistency: edge " << e.getInnerId() << " " << e.getStart().getInnerId() << " -> "
               << e.getFinish().getInnerId() << " must hold a record but has none\n";
            errors++;
        }
    }
    return errors == 0;
}
