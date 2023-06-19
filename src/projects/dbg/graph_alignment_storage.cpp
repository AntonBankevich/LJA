#include "graph_alignment_storage.hpp"

using namespace dbg;
void AlignedRead::correct(CompactPath &&cpath) {
    VERIFY_MSG(!corrected, "Attempt to correct path while previous correction was not yet applied");
    corrected_path = std::move(cpath);
    corrected = true;
}

void AlignedRead::applyCorrection() {
    if(corrected)
        path = std::move(corrected_path);
    corrected_path = {};
    corrected = false;
}

void AlignedRead::delayedInvalidate() {
    VERIFY(!corrected);
    corrected_path = {};
    corrected = true;
}

void VertexRecord::addPath(const Sequence &seq) {
    lock();
    cov += 1;
    for(std::pair<Sequence, size_t> &path : paths) {
        if(path.first == seq) {
            path.second += 1;
            if(path.second == 1)
                zero_cnt -= 1;
            unlock();
            return;
        }
    }
    paths.emplace_back(seq, 1);
    unlock();
}

void VertexRecord::removePath(const Sequence &seq) {
    lock();
    bool found = false;
    for(std::pair<Sequence, size_t> &path : paths) {
        if(path.first == seq) {
            found = true;
            VERIFY(path.second > 0);
            path.second -= 1;
            cov -= 1;
            if(path.second == 0)
                zero_cnt += 1;
            break;
        }
    }
    if(!found) {
        std::cout << "Error" << std::endl;
        unlock();
        std::cout << seq << std::endl;
        std::cout << this->str() << std::endl;
    }
    VERIFY(found);
    if(zero_cnt > paths.size() / 3) {
        std::vector<std::pair<Sequence, size_t>> new_paths;
        for(std::pair<Sequence, size_t> &rec : paths) {
            if(rec.second != 0) {
                new_paths.emplace_back(std::move(rec.first), rec.second);
            }
        }
        std::swap(paths, new_paths);
        zero_cnt = 0;
    }
    unlock();
}

bool VertexRecord::isDisconnected(const Edge &edge) const {
    if(edge.end()->outDeg() == 0)
        return false;
    for(const std::pair<Sequence, size_t> &rec : paths) { // NOLINT(readability-use-anyofallof)
        if(rec.second > 0 && rec.first[0] == edge.seq[0] && rec.first.size() > 1)
            return false;
    }
    return true;
}

size_t VertexRecord::countStartsWith(const Sequence &seq) const {
//        lock();
    size_t cnt = 0;
    for(const std::pair<Sequence, size_t> &rec : paths) {
        if(rec.first.startsWith(seq)) {
            cnt += rec.second;
        }
    }
//        unlock();
    return cnt;
}

std::vector<GraphAlignment> VertexRecord::getBulgeAlternatives(const Vertex &end, double threshold) const {
//        lock();
    std::vector<std::pair<Sequence, size_t>> candidates;
    for(const auto & extension : paths) {
        Path unpacked = CompactPath(v, extension.first).getPath();
        for(size_t i = 1; i <= unpacked.size(); i++) {
            if(end == unpacked.getVertex(i)){
                candidates.emplace_back(extension.first.Subseq(0, i), extension.second);
            }
        }
    }
//        unlock();
    if(candidates.empty())
        return {};
    std::sort(candidates.begin(), candidates.end());
    std::vector<GraphAlignment> res;
    size_t cnt = 0;
    for(size_t i = 0; i < candidates.size(); i++) {
        if(i > 0 && candidates[i-1].first != candidates[i].first) {
            if(cnt >= threshold)
                res.emplace_back(CompactPath(v, candidates[i - 1].first).getAlignment());
            cnt = 0;
        }
        cnt += candidates[i].second;
    }
    if(cnt >= threshold)
        res.emplace_back(CompactPath(v, candidates.back().first).getAlignment());
    return std::move(res);
}

CompactPath VertexRecord::getFullUniqueExtension(const Sequence &start, size_t min_good_cov, size_t max_bad_cov, size_t max_size) const {
    Sequence res = start;
    while(res.size() < max_size) {
        unsigned char next = getUniqueExtension(res, min_good_cov, max_bad_cov);
        if(next == (unsigned char)-1)
            break;
        std::vector<unsigned char> ext = {next};
        res = res + Sequence(ext);
    }
    return {v, res};
}

unsigned char VertexRecord::getUniqueExtension(const Sequence &start, size_t min_good, size_t max_bad) const {
    std::vector<size_t> counts(4);
    for(const auto & extension : paths) {
        if(extension.first.size() > start.size() && extension.first.startsWith(start)) {
            counts[extension.first[start.size()]] += extension.second;
        }
    }
    size_t bad = 0;
    size_t good = 0;
    size_t res = 0;
    for(size_t c = 0; c < 4; c++) {
        if(counts[c] <= max_bad)
            bad++;
        if(counts[c] >= min_good) {
            good++;
            res = c;
        }
    }
    if(bad != 3 || good != 1)
        return (unsigned char)(-1);
    return (unsigned char)(res);
}

std::vector<GraphAlignment> VertexRecord::getTipAlternatives(size_t len, double threshold) const {
    len += std::max<size_t>(30, len / 20);
//        lock();
    std::vector<std::pair<Sequence, size_t>> candidates;
    for(const auto & extension : paths) {
        GraphAlignment unpacked = CompactPath(v, extension.first).getAlignment();
        if(unpacked.len() >= len) {
            unpacked.cutBack(unpacked.len() - len);
            candidates.emplace_back(CompactPath(unpacked).cpath(), extension.second);
        }
    }
//        unlock();
    if(candidates.empty())
        return {};
    std::sort(candidates.begin(), candidates.end());
    std::vector<GraphAlignment> res;
    size_t cnt = 0;
    for(size_t i = 0; i < candidates.size(); i++) {
        if(i > 0 && candidates[i - 1].first != candidates[i].first) {
            if(cnt > threshold) {
                GraphAlignment cp = CompactPath(v, candidates[i - 1].first).getAlignment();
                cp.cutBack(cp.len() - len);
                res.emplace_back(cp);
            }
            cnt = 0;
        }
        cnt += candidates[i].second;
    }
    if(cnt > threshold) {
        GraphAlignment cp = CompactPath(v, candidates.back().first).getAlignment();
        cp.cutBack(cp.len() - len);
        res.emplace_back(cp);
    }
    return std::move(res);
}

std::string VertexRecord::str() const {
    std::stringstream ss;
    lock();
    for(const auto & path : paths) {
        ss << path.first << " " << path.second << std::endl;
    }
    unlock();
    return ss.str();
}

ReadLogger::CountingSS &ReadLogger::CountingSS::operator<<(const string &s) {
    log << s;
    len += s.size();
    return *this;
}

ReadLogger::CountingSS &ReadLogger::CountingSS::operator<<(const size_t &s) {
    log << s;
    len += 5;
    return *this;
}

void ReadLogger::CountingSS::clear() {
    log = std::stringstream();
    len = 0;
}

void ReadLogger::dump(ReadLogger::CountingSS &sublog) {
#pragma omp critical
    {
        os << sublog.str();
    };
    sublog.clear();
}

void ReadLogger::flush() {
    for(CountingSS &sublog : logs) {
        dump(sublog);
    }
}

ReadLogger::~ReadLogger() {
    flush();
    os.close();
}

void ReadLogger::logRead(AlignedRead &alignedRead) {
    CountingSS &ss = logs[omp_get_thread_num()];
    ss << alignedRead.id << " initial " << alignedRead.path.getAlignment().str(true) << "\n";
    if(ss.size() > 100000) {
        dump(ss);
    }
}

void ReadLogger::logRerouting(AlignedRead &alignedRead, const GraphAlignment &initial, const GraphAlignment &corrected,
                              const string &message) {
    CountingSS &ss = logs[omp_get_thread_num()];
    size_t left = 0;
    size_t right = 0;
    size_t left_len = 0;
    size_t right_len = 0;
    while(left < corrected.size() && left < initial.size()) {
        if(initial[left] != corrected[left])
            break;
        left_len += initial[left].size();
        left++;
    }
    while(left + right < corrected.size() && left + right < initial.size()) {
        if(initial[initial.size() - right - 1] != corrected[corrected.size() - right - 1])
            break;
        right_len += initial[initial.size() - right - 1].size();
        right++;
    }
    ss << alignedRead.id << " " << message  << " " << left << "(" << left_len << ") " << right << "(" << right_len << ")\n";
    ss << alignedRead.id << "  initial  " << initial.subalignment(left, initial.size() - right).str(true) << "\n";
    ss << alignedRead.id << " corrected " << corrected.subalignment(left, corrected.size() - right).str(true) << "\n";
//        ss << alignedRead.id << " rc  initial  " << initial.RC().str(true) << "\n";
//        ss << alignedRead.id << " rc corrected " << corrected.RC().str(true) << "\n";
    if(ss.size() > 100000) {
        dump(ss);
    }
}

void RecordStorage::processPath(const CompactPath &cpath, const std::function<void(Vertex &, const Sequence &)> &task,
                                const std::function<void(Segment<Edge>)> &edge_task) const {
    GraphAlignment al = cpath.getAlignment();
    for(size_t i = 0; i < al.size(); i++) {
        Edge &edge = al[i].contig();
        size_t seg_left = i == 0 ? cpath.leftSkip() : 0;
        size_t seg_right = i == cpath.size() - 1 ? edge.size() - cpath.rightSkip() : edge.size();
        edge_task(Segment<Edge>(edge, seg_left, seg_right));
    }
    size_t j = 1;
    size_t clen = al[0].contig().size();
    for (size_t i = 1; i <= al.size(); i++) {
        clen -= al[i - 1].contig().size();
        while (j < al.size() && clen < max_len) {
            clen += al[j].contig().size();
            j++;
        }
        if (clen >= min_len)
            task(al.getVertex(i - 1), cpath.cpath().Subseq(i - 1, j));
    }
}

//TODO Remove threads parameter
RecordStorage::RecordStorage(SparseDBG &dbg, size_t _min_len, size_t _max_len, size_t threads,
                             ReadLogger &readLogger, bool _track_cov, bool log_changes, bool track_suffixes) :
        min_len(_min_len), max_len(_max_len), track_cov(_track_cov), readLogger(&readLogger), log_changes(log_changes), track_suffixes(track_suffixes) {
    for(auto &it : dbg) {
        data.emplace(&it.second, VertexRecord(it.second));
        data.emplace(&it.second.rc(), VertexRecord(it.second.rc()));
    }
}

std::function<std::string(Edge &)> RecordStorage::labeler() const {
    if(track_suffixes)
        return [this](Edge &edge) {
            const VertexRecord &rec = getRecord(*edge.start());
            std::stringstream ss;
            size_t cnt = 0;
            for (const auto &ext : rec) {
                if (ext.first[0] == edge.seq[0]) {
                    if(cnt < 30)
                        ss << ext.first << "(" << ext.second << ")\\n";
                    cnt++;
                }
            }
            if(cnt > 30) {
                ss << "and another " << (cnt - 30) << " records\\n";
            }
            return ss.str();
        };
    else return [](Edge &) {
        return "";
    };
}

void RecordStorage::addRead(AlignedRead &&read) {
    reads.emplace_back(std::move(read));
    addSubpath(read.path);
    addSubpath(read.path.RC());
}

void RecordStorage::delayedInvalidateRead(AlignedRead &read, const std::string &message) { // NOLINT(readability-convert-member-functions-to-static)
    if(log_changes)
        readLogger->logInvalidate(read, message);
    read.delayedInvalidate();
}

void RecordStorage::delayedInvalidateBad(logging::Logger &logger, size_t threads, double threshold, const std::string &message) {
    const std::function<bool(const Edge &)> &is_bad = [threshold](const Edge &edge) {
        return edge.getCoverage() < threshold;
    };
    delayedInvalidateBad(logger, threads, is_bad, message);
}

void RecordStorage::delayedInvalidateBad(logging::Logger &logger, size_t threads, const std::function<bool(
        const Edge &)> &is_bad, const std::string &message) {
    ParallelCounter cnt(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(cnt, message, is_bad)
    for(size_t i = 0; i < reads.size(); i++) {
        AlignedRead &alignedRead = reads[i];
        if(!alignedRead.valid())
            continue;
        GraphAlignment al = alignedRead.path.getAlignment();
        size_t l = 0;
        size_t r = al.size();
        while(l < al.size() && is_bad(al[l].contig())) {
            l++;
        }
        while(r > 0 && is_bad(al[r - 1].contig())) {
            r--;
        }
        bool middle_bad = false;
        size_t len = 0;
        for(size_t j = l; j < r; j++) {
            if(is_bad(al[j].contig())) {
                middle_bad = true;
                break;
            }
            len += al[j].size();
        }
        if(middle_bad || l >= r) {
            delayedInvalidateRead(alignedRead, message);
            cnt += 1;
        } else if(r - l != al.size()) {
            GraphAlignment sub = al.subalignment(l, r);
            if(sub.len() < 1000)
                delayedInvalidateRead(alignedRead, message);
            else {
                std::string extra_message = message + "_EndsClipped_" + itos(al.subalignment(0, l).len()) + "_" + itos(al.subalignment(r, al.size()).len());
                reroute(alignedRead, al, al.subalignment(l, r), extra_message);
            }
            cnt += 1;
        }
    }
    logger.info() << "Could not correct " << cnt.get() << " reads. They were removed or truncated." << std::endl;
}

void RecordStorage::invalidateSubreads(logging::Logger &logger, size_t threads) {
    for(AlignedRead &alignedRead : reads) {
        const VertexRecord &rec = this->getRecord(alignedRead.path.start());
        size_t cnt = rec.countStartsWith(alignedRead.path.cpath());
        VERIFY_MSG(cnt >= 1, "This function assumes that suffixes are stored for complete reads")
        if(rec.countStartsWith(alignedRead.path.cpath()) >= 2) {
            delayedInvalidateRead(alignedRead, "Subread");
        }
    }
    logger.info() << "Uncorrected reads were removed." << std::endl;
}

void RecordStorage::addSubpath(const CompactPath &cpath) {
    if(!cpath.valid())
        return;
    std::function<void(Vertex &, const Sequence &)> vertex_task = [](Vertex &v, const Sequence &s) {};
    if(track_suffixes)
        vertex_task = [this](Vertex &v, const Sequence &s) {
            data.find(&v)->second.addPath(s);
        };
    std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg){};
    if(track_cov)
        edge_task = [](Segment<Edge> seg){
            seg.contig().incCov(seg.size());
        };
    processPath(cpath, vertex_task, edge_task);
}

void RecordStorage::removeSubpath(const CompactPath &cpath) {
    if(!cpath.valid())
        return;
    std::function<void(Vertex &, const Sequence &)> vertex_task = [](Vertex &v, const Sequence &s) {};
    if(track_suffixes)
        vertex_task = [this](Vertex &v, const Sequence &s) {
            data.find(&v)->second.removePath(s);
        };
    std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg){};
    if(track_cov)
        edge_task = [](Segment<Edge> seg) {
            seg.contig().incCov(size_t(-seg.size()));
        };
    processPath(cpath, vertex_task, edge_task);
}

bool RecordStorage::apply(AlignedRead &alignedRead) {
    if(!alignedRead.checkCorrected())
        return false;
    this->removeSubpath(alignedRead.path);
    this->removeSubpath(alignedRead.path.RC());
    alignedRead.applyCorrection();
    this->addSubpath(alignedRead.path);
    this->addSubpath(alignedRead.path.RC());
    return true;
}

void RecordStorage::reroute(AlignedRead &alignedRead, const GraphAlignment &initial, const GraphAlignment &corrected,
                            const string &message) {
    if (log_changes)
        readLogger->logRerouting(alignedRead, initial, corrected, message);
    if(corrected.len() < 500)
        delayedInvalidateRead(alignedRead, "Deleted_after_" + message);
    else {
        alignedRead.correct(CompactPath(corrected));
    }
}

void RecordStorage::reroute(AlignedRead &alignedRead, const GraphAlignment &corrected, const string &message) {
    reroute(alignedRead, alignedRead.path.getAlignment(), corrected, message);
}

void RecordStorage::applyCorrections(logging::Logger &logger, size_t threads) {
    if(size() > 10000)
        logger.info() << "Applying corrections to reads" << std::endl;
    omp_set_num_threads(threads);
    ParallelCounter cnt(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(cnt)
    for(size_t i = 0; i < reads.size(); i++) { // NOLINT(modernize-loop-convert)
        if(apply(reads[i]))
            cnt += 1;
    }
    flush();
    if(size() > 10000)
        logger.info() << "Applied correction to " << cnt.get() << " reads" << std::endl;
}

void RecordStorage::printReadAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const {
    logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
    std::string acgt = "ACGT";
    std::ofstream os;
    os.open(path);
    for(const AlignedRead &read : reads) {
        const CompactPath& al = read.path;
        if(!al.valid())
            continue;
        os  << read.id << " " << al.start().hash() << int(al.start().isCanonical())
            << " " << al.cpath().str() << "\n";
        CompactPath rc_al = al.RC();
        os  << "-" << read.id << " " << rc_al.start().hash() << int(rc_al.start().isCanonical())
            << " " << rc_al.cpath().str() << "\n";
    }
    os.close();
}

void RecordStorage::printReadFasta(logging::Logger &logger, const std::experimental::filesystem::path &path) const {
    logger.info() << "Printing reads to fasta file " << path << std::endl;
    std::string acgt = "ACGT";
    std::ofstream os;
    os.open(path);
    for(const AlignedRead &read : reads) {
        const CompactPath &al = read.path;
        if(!al.valid())
            continue;
        os  << ">" << read.id << "\n" << read.path.getAlignment().Seq() << "\n";
    }
    os.close();
}

void RecordStorage::printFullAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const {
    logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
    std::string acgt = "ACGT";
    std::ofstream os;
    os.open(path);
    for(const AlignedRead &read : reads) {
        const CompactPath &al = read.path;
        if(!al.valid())
            continue;
        os << read.id << " " << read.path.getAlignment().str(true) << "\n";
        os << "-" << read.id << " " << read.path.getAlignment().RC().str(true) << "\n";
    }
    os.close();
}

//void RecordStorage::updateExtensionSize(logging::Logger &logger, size_t threads, size_t new_max_extension) {
//    logger << "Updating stored read subpaths size" << std::endl;
//    for(auto &it : this->data) {
//        it.second.clear();
//    }
//    this->max_len = new_max_extension;
//    bool tmp_track_cov = track_cov;
//    track_cov = false;
//#pragma omp parallel for default(none) schedule(dynamic, 100)
//    for(size_t i = 0; i < size(); i++) {
//        addSubpath(reads[i].path);
//        addSubpath(reads[i].path.RC());
//    }
//    track_cov = tmp_track_cov;
//}

const VertexRecord &RecordStorage::getRecord(const Vertex &v) const {
    VERIFY(track_suffixes);
    return data.find(&v)->second;
}

void RecordStorage::trackSuffixes(logging::Logger &logger, size_t threads) {
    logger.info() << "Collecting and storing read suffixes" << std::endl;
    VERIFY(!track_suffixes);
    track_suffixes = true;
    omp_set_num_threads(threads);
    std::function<void(Vertex &, const Sequence &)> vertex_task  = [this](Vertex &v, const Sequence &s) {
        data.find(&v)->second.addPath(s);
    };
    std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg){};
#pragma omp parallel for default(none) shared(vertex_task, edge_task)
    for(size_t i = 0; i < reads.size(); i++) {
        if(reads[i].valid()) {
            processPath(reads[i].path, vertex_task, edge_task);
            processPath(reads[i].path.RC(), vertex_task, edge_task);
        }
    }
}

void RecordStorage::untrackSuffixes() {
    if(track_suffixes) {
        track_suffixes = false;
        for (auto &it: this->data) {
            it.second.clear();
        }
    }
}

void RecordStorage::Save(std::ostream &os) const {
    os << size() << "\n";
    for(const AlignedRead &alignedRead : *this) {
        os << alignedRead << "\n";
    }
}

void RecordStorage::Load(std::istream &is, SparseDBG &dbg) {
    size_t sz;
    is >> sz;
    for(size_t i = 0; i < sz; i++) {
        addRead(AlignedRead::Load(is, dbg));
    }
}

void SaveAllReads(const std::experimental::filesystem::path &fname, const std::vector<RecordStorage *> &recs) {
    std::ofstream os;
    os.open(fname);
    os << recs.size() << "\n";
    for(RecordStorage *rs : recs) {
        rs->Save(os);
    }
    os.close();
}

void LoadAllReads(const std::experimental::filesystem::path &fname, const std::vector<RecordStorage *> &recs,
                  dbg::SparseDBG &dbg) {
    std::ifstream is;
    is.open(fname);
    size_t sz;
    is >> sz;
    VERIFY(sz == recs.size())
    for(RecordStorage *recordStorage : recs) {
        recordStorage->Load(is, dbg);
    }
    is.close();
}
