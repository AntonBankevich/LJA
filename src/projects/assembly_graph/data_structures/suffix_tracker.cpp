#include "suffix_tracker.hpp"
using namespace ag;

size_t SuffixRecord::getNumberOfEnds() const {
    VERIFY(!getEdge().isSuffix() || num_of_ends == 0);
    return num_of_ends;
}

size_t SuffixRecord::getNumberOfPaths() const {
    size_t res = 0;
    for(const auto &p : paths)
        res += p.second;
    VERIFY(res == num_of_paths);
    return res + num_of_ends;
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
    if(seq.empty())
        cnt += num_of_ends;
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
        VERIFY(num_of_ends + diff >= 0);
        num_of_ends += diff;
        return;
    }
    num_of_paths += diff;
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
    if (seq.empty()) {
        num_of_ends += cnt;
    } else {
        num_of_paths += cnt;
        paths.emplace_back(seq, cnt);
    }
}

void SuffixRecord::clear() {
    num_of_paths = 0;
    num_of_ends = 0;
    paths.clear();
}

void SuffixRecord::resetCodes(Vertex &start) {
    std::vector<std::pair<Sequence, int>> old = std::move(paths);
    size_t old_num_of_paths = num_of_paths;
    num_of_paths = 0;
    for(const std::pair<Sequence, int> &rec : old) {
        if(rec.second == 0)
            continue;
        GraphPath path(start, rec.first);
        path.resetEdgeCodes();
        changePathCnt(Sequence(path.getFSplits()), Sequence(path.getFSplits()), rec.second);
    }
    VERIFY(old_num_of_paths == num_of_paths);
}

std::function<std::string(const Edge & )> SuffixTracker::labeler() const {
    return [this](const Edge &edge) {
        if (edge.isPrefix()) return std::string();
        const SuffixRecord &rec = getSuffixRecord(edge);
        std::stringstream ss;
        size_t cnt = 0;
        ss << "_(" << rec.num_of_ends << ")";
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
    PathPosition from_pos = left;
    PathPosition to_pos = from_pos + 1;
    size_t clen = left.nextEdge().truncSize();
    VERIFY(to_pos <= right);
    Sequence read_seq(left.getFPos(), right.getFPos());
    while (from_pos != right) {
        Edge &edge = from_pos.nextEdge();
        clen -= edge.truncSize();
        ++from_pos;
        if(!edge.isPrefix()) {
            SuffixRecord &erec = getSuffixRecord(edge);
            while (to_pos < right && (clen < erec.getMaxSuffixLength())) {
                Edge &new_edge = to_pos.nextEdge();
                clen += new_edge.truncSize();
                to_pos += new_edge;
            }
            VERIFY(to_pos <= right);
            VERIFY(!to_pos.prevEdge().isSuffix());
//            if(to_pos > right)
//                to_pos = right;
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
        if(!e.isPrefix())
            edge_data.insert(e.getId(), std::make_unique<SuffixRecord>(e, max_len));
    }
}

void SuffixTracker::fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
    if (new_edge.isPrefix())
        return;
    int ends = 0;
    for(Edge &e : path.edges()) {
        if (!e.isPrefix())
            ends += getSuffixRecord(e).num_of_ends;
    }
    SuffixRecord &erec = getSuffixRecord(new_edge);
    SuffixRecord &lasterec = getSuffixRecord(path.backEdge());
    erec.paths = std::move(lasterec.paths);
    erec.num_of_ends = ends;
    erec.num_of_paths = lasterec.num_of_paths;
}

void SuffixTracker::fireMergePath(const RAGraphPath &path, Vertex &vertex) {
    VERIFY(vertex != path.getFinish());
    VERIFY(vertex != path.getStart());
    VERIFY(!path.backEdge().isPrefix());
    SuffixRecord &old_rec = getSuffixRecord(path.backEdge());
    SuffixRecord &new_rec = getSuffixRecord(vertex.front());
    new_rec.paths = std::move(old_rec.paths);
    new_rec.num_of_ends = old_rec.num_of_ends;
    new_rec.num_of_paths = old_rec.num_of_paths;

    // RAGraphPath rc_path = path.RC();
    // for(Edge &e : rc_path.edges()) {
    //     if (!e.rc().isPrefix()) {
    //         edge_data.at(vertex.front().getId()).paths = std::move(edge_data.at(e.rc().getId()).paths);
    //         return;
    //     }
    // }
}

void
SuffixTracker::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &,
                                   const AlignmentForm &) {
    SuffixRecord &erec = getSuffixRecord(new_edge);
    SuffixRecord &rrec = getSuffixRecord(right);
    SuffixRecord &lrec = getSuffixRecord(left);
    VERIFY(lrec.num_of_paths == 0);//Since lrec is a tip, it can have no recorded continuations
    erec.paths = std::move(rrec.paths);
    erec.num_of_ends = lrec.num_of_ends + rrec.num_of_ends;
    erec.num_of_paths = rrec.num_of_paths;
}

void SuffixTracker::fireSplitEdge(Edge &edge, const RAGraphPath &split) {
    SuffixRecord &erec = getSuffixRecord(edge);
    SuffixRecord &last_rec = getSuffixRecord(split.backEdge());
    last_rec.paths = std::move(erec.paths);
    for(Edge &e : split.edges()) {
        getSuffixRecord(e.rc()).lockFreeChangePathCnt(Sequence(), Sequence(), storage->startCnt(e));
    }
    last_rec.num_of_paths = erec.num_of_paths;
}

void SuffixTracker::fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
    std::function<void(size_t, Edge &)> task = [this](size_t, Edge &e) {
        if (!e.isPrefix())
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
            VERIFY(!p.first->isPrefix());
            SuffixRecord &rec = *p.second;
            VERIFY(p.first->isSuffix() || rec.num_of_ends == storage->startCnt(p.first->rc()));
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

void SuffixTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    SuffixRecord &erec = getSuffixRecord(e);
    SuffixRecord & new_rec = getSuffixRecord(v.front());
    new_rec.paths = std::move(erec.paths);
    new_rec.num_of_paths = erec.num_of_paths;
}

void SuffixTracker::fireAddEdge(Edge &edge) {
    if (edge.isPrefix())
        return;
    edge_data.insert(edge.getId(), std::make_unique<SuffixRecord>(edge, max_len));
}

void SuffixTracker::fireDeleteEdge(Edge &edge) {
    if (edge.isPrefix())
        return;
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
}
