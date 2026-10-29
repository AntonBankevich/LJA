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
    for (auto & path : paths) {
        if((path.first == min_seq || path.first == max_seq) && path.second + diff >= 0) {
            path.second += diff;
            updateZero(path.second - diff, path.second);
            diff = 0;
            return;
        }
    }
    std::vector<size_t> subseqs;
    for (size_t i = 0; i < paths.size(); i++) {
        std::pair<Sequence, int> &path = paths[i];
        if(path.second != 0 && max_seq.startsWith(path.first) && path.first.startsWith(min_seq) ) {
            VERIFY(path.second > 0);
            if (path.second + diff >= 0) {
                path.second += diff;
                updateZero(path.second - diff, path.second);
                diff = 0;
                return;
            } else {
                diff += path.second;
                path.second = 0;
                updateZero(1, 0);
            }
            subseqs.push_back(i);
        }
    }
    VERIFY_MSG(diff == 0, "Attempting to remove path that is not present in the record.");
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
    while (right != left && right.prevEdge().isSuffix()) {
        --right;
    }
    PathPosition from_pos = left;
    PathPosition to_pos = from_pos + 1;
    size_t clen = left.nextEdge().truncSize();
    VERIFY(to_pos <= right);
    Sequence read_seq(left.getFPos(), right.getFPos());
    while (from_pos != right) {
        Edge &edge = from_pos.nextEdge();
        clen -= edge.truncSize();
        ++from_pos;
        SuffixRecord &erec = getSuffixRecord(edge);
        while (to_pos < right && (clen < erec.getMaxSuffixLength())) {
            Edge &new_edge = to_pos.nextEdge();
            clen += new_edge.truncSize();
            to_pos += new_edge;
        }
        VERIFY(to_pos <= right);
        VERIFY(!to_pos.prevEdge().isSuffix());
        if (clen >= min_len) {
            if (from_pos == to_pos) {
                erec.changePathCnt(Sequence(), Sequence(), diff);
            } else {
                Sequence subseq_max_code = read_seq.Subseq(from_pos.getFPos() - left.getFPos(),
                        to_pos.getFPos() - left.getFPos());
                erec.changePathCnt(subseq_max_code.Subseq(0, subseq_max_code.size() - to_pos.prevEdge().getCode().size() + 1),
                    subseq_max_code, diff);
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

void SuffixTracker::fireAddRead(const AlignedRead &read) {
    if(read.valid()) {
        addSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
        addSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
    }
}

void SuffixTracker::fireRerouteRead(AlignedRead &read) {
    if(read.valid()) {
        removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
        removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
    }
    if(read.getCorrected().valid()) {
        addSubpath(read.getCorrected().firstPosition(), read.getCorrected().lastPosition());
        addSubpath(read.getCorrected().lastPosition().RC(), read.getCorrected().firstPosition().RC());
    }
}

void SuffixTracker::fireInvalidateRead(AlignedRead &read) {
    if(read.valid()) {
        removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
        removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
    }
}

SuffixTracker::SuffixTracker(AlignedReadStorage &storage, AssemblyGraph &graph,
                             size_t _min_len, size_t _max_len) :
        AlignedReadStorageListener(storage, "SuffixTracker"), ResolutionListener(graph, "SuffixTracker"), storage(&storage),
        min_len(_min_len), max_len(_max_len) {
    for(Edge &e : graph.edges()) {
        edge_data.insert(e.getId(), std::make_unique<SuffixRecord>(e, max_len));
    }
}

void SuffixTracker::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    SuffixRecord &erec = getSuffixRecord(new_edge);
    SuffixRecord &lasterec = getSuffixRecord(path.backEdge());
    erec.paths = std::move(lasterec.paths);
}

void SuffixTracker::fireMergePath(const RAGraphPath &path, Vertex &vertex) {
    VERIFY(vertex != path.getFinish());
    VERIFY(vertex != path.getStart());
    //Unlike a resolved vertex, a merged vertex can have reads starting inside it (at one of the path's
    //former inner vertices), so incoming can not be derived from outgoing. Instead each boundary edge of
    //vertex inherits from the matching boundary edge of path: path.frontEdge() occupied the same entrance
    //position before the collapse, path.backEdge() the same exit position.
    for (Edge &edge: path.getFinish()) {
        Edge &new_edge = vertex.getNonCorporealOutgoing(edge.firstNucl());
        getSuffixRecord(new_edge).paths = std::move(getSuffixRecord(edge).paths);
        getSuffixRecord(new_edge).zero_cnt = getSuffixRecord(edge).zero_cnt;
    }
    for (Edge &edge: path.getStart().incoming()) {
        Edge & new_edge = vertex.rc().getNonCorporealOutgoing(edge.rc().firstNucl()).rc();
        SuffixRecord &new_rec = getSuffixRecord(new_edge);
        for (std::pair<Sequence, int> &rec: getSuffixRecord(edge).paths) {
            if (rec.first.size() > new_edge.getCode().size() - edge.getCode().size()) {
                new_rec.directAddPath(rec.first.Subseq(new_edge.getCode().size() - edge.getCode().size()), rec.second);
            }
        }
    }
}

void
SuffixTracker::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &,
                                   const AlignmentForm &) {
    SuffixRecord &erec = getSuffixRecord(new_edge);
    SuffixRecord &rrec = getSuffixRecord(right);
    SuffixRecord &lrec = getSuffixRecord(left);
    erec.paths = std::move(rrec.paths);
}

void SuffixTracker::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
    SuffixRecord &erec = getSuffixRecord(edge);
    SuffixRecord &last_rec = getSuffixRecord(split.backEdge());
    last_rec.paths = std::move(erec.paths);
}

void SuffixTracker::fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
    std::function<void(size_t, Edge &)> task = [this](size_t, Edge &e) {
        getSuffixRecord(e).resetCodes(e.getFinish());
    };
    processObjects(graph.edges().begin(), graph.edges().end(), logger, threads, task);
}

bool SuffixTracker::fireCheckConsistency() {
//        fireCheckConsistency is only ever called single-threaded, never concurrently with graph
//        modifications; lock_table() here is just libcuckoo's only full-table iteration API, not
//        synchronization against any other running thread.
    {
        auto lt = edge_data.lock_table();
        for (auto &p : lt) {
            SuffixRecord &rec = *p.second;
        }
    }
    for(AlignedRead &read : *storage) {
        if(read.valid()) {
            removeSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
            removeSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
        }
    }
    {
        auto lt = edge_data.lock_table();
        for (auto &p : lt) {
            SuffixRecord &rec = *p.second;
            rec.removeZero();
            VERIFY(rec.begin() == rec.end());
            if(rec.begin() != rec.end())
                return false;
        }
    }
    for(AlignedRead &read : *storage) {
        if(read.valid()) {
            addSubpath(read.getPath().firstPosition(), read.getPath().lastPosition());
            addSubpath(read.getPath().lastPosition().RC(), read.getPath().firstPosition().RC());
        }
    }
    return true;
}

//e (X->Y) is replaced by X->v (prefix edge = v.incFront()) and v->Y (suffix edge = v.front()), which
//together carry exactly e's code. Every read suffix stored on e continues identically after both of v's
//boundary edges - a read entering v via the prefix edge, and (for reads continuing past Y) the suffix
//edge, pick up the same continuation e had. Copy e's record onto both so a tracker carried through the
//conversion matches a from-scratch fill on the supregraph. The rc-mirrored fire handles v.rc()'s pair.
void SuffixTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    SuffixRecord &erec = getSuffixRecord(e);
    for (Edge *dst : {&v.front(), &v.incFront()}) {
        SuffixRecord &drec = getSuffixRecord(*dst);
        drec.paths = erec.paths;
        drec.zero_cnt = erec.zero_cnt;
    }
}

void SuffixTracker::fireAddEdge(Edge &edge) {
    edge_data.insert(edge.getId(), std::make_unique<SuffixRecord>(edge, max_len));
}

void SuffixTracker::fireDeleteEdge(Edge &edge) {
    edge_data.erase(edge.getId());
}

void SuffixTracker::fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
    for (Edge &inc: core.incoming()) {
        const SuffixRecord &rec = getSuffixRecord(inc);
        std::unordered_map<unsigned char, std::pair<size_t, SuffixRecord *>> rec_map;
        for(Edge &out : core) {
            if(!resolution.contains(inc, out))
                continue;
            rec_map[out.firstNucl()] = {out.getCode().size(), &getSuffixRecord(resolution.get(inc, out).front())};
            size_t len = getSuffixRecord(inc).max_suffix_len;
            len -= std::min(len, out.truncSize());
            rec_map[out.firstNucl()].second->max_suffix_len = len;
        }
        for (const std::pair<Sequence, int> &info: rec) {
            if (info.second == 0)
                continue;
            VERIFY(!info.first.empty());
            Edge &out = core.getOutgoing(info.first[0]);
            if(out.getCode().size() < info.first.size()) {
                rec_map.at(info.first[0]).second->directAddPath(info.first.Subseq(rec_map.at(info.first[0]).first), info.second);
            }
        }
    }
    //There can be no read that starts in any of the new vertices. Thus incoming record is the same as outgoing record
    for (auto &rec : resolution) {
        SuffixRecord &inc_rec = getSuffixRecord(rec.first->incFront());
        SuffixRecord &out_rec = getSuffixRecord(rec.first->front());
        inc_rec.paths = out_rec.paths;
    }
}
