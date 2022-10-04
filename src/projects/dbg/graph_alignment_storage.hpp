#pragma once

#include "compact_path.hpp"

class AlignedRead {
private:
    dbg::CompactPath corrected_path;
    bool corrected = false;
public:
    std::string id;
    dbg::CompactPath path;

    AlignedRead() = default;
    AlignedRead(AlignedRead &&other) = default;
    AlignedRead &operator=(AlignedRead &&other) = default;
    explicit AlignedRead(std::string readId) : id(std::move(readId)), corrected(false) {}
    AlignedRead(std::string readId, dbg::GraphAlignment &_path) : id(std::move(readId)), path(_path), corrected(false) {}
    AlignedRead(std::string readId, dbg::CompactPath _path) : id(std::move(readId)), path(std::move(_path)), corrected(false) {}

    bool operator<(const AlignedRead& other) const {return id < other.id;}

    void delayedInvalidate();
    bool checkCorrected() const {return corrected;}
    bool valid() const {return path.valid();}

    void correct(dbg::CompactPath &&cpath);
    void applyCorrection();

    static AlignedRead Load(std::istream &is, dbg::SparseDBG &dbg) {
        std::string id;
        is >> id;
        return {id, dbg::CompactPath::Load(is, dbg)};
    }
};

inline std::ostream& operator<<(std::ostream  &os, const AlignedRead &alignedRead) {
    return os << alignedRead.id << " " << alignedRead.path;
}

class RecordStorage;
struct VertexRecord {
    friend RecordStorage;
private:
    typedef std::vector<std::pair<Sequence, size_t>> Storage;
    typedef Storage::const_iterator const_iterator;
    dbg::Vertex &v;
    Storage paths;
    size_t zero_cnt = 0;
    size_t cov = 0;

    void lock() const {v.lock();}
    void unlock() const {v.unlock();}

    void addPath(const Sequence &seq);
    void removePath(const Sequence &seq);
    void clear() {paths.clear();}
public:
    explicit VertexRecord(dbg::Vertex &_v) : v(_v) {}
    VertexRecord(const VertexRecord &) = delete;
    VertexRecord(VertexRecord &&other)  noexcept : v(other.v), paths(std::move(other.paths)),
                                                   zero_cnt(other.zero_cnt), cov(other.cov) {}

    VertexRecord & operator=(const VertexRecord &) = delete;

    size_t coverage() const {return cov;}
    std::string str() const;
    const_iterator begin() const {return paths.begin();}
    const_iterator end() const {return paths.end();}

    size_t countStartsWith(const Sequence &seq) const;

    bool isDisconnected(const dbg::Edge &edge) const;
    std::vector<dbg::GraphAlignment> getBulgeAlternatives(const dbg::Vertex &end, double threshold) const;
    std::vector<dbg::GraphAlignment> getTipAlternatives(size_t len, double threshold) const;
    unsigned char getUniqueExtension(const Sequence &start, size_t min_good_cov, size_t max_bad_cov) const;
    dbg::CompactPath getFullUniqueExtension(const Sequence &start, size_t min_good_cov, size_t max_bad_cov, size_t max_size = size_t(-1)) const;
};

inline std::ostream& operator<<(std::ostream  &os, const VertexRecord &rec) {return os << rec.str();}

class ReadLogger {
private:
    class CountingSS {
    private:
        std::stringstream log;
        size_t len;
    public:
        CountingSS() : log(), len(0) {}

        std::string str() {return log.str();}
        size_t size() const {return len;}

        CountingSS &operator<<(const std::string &s);
        CountingSS &operator<<(const size_t &s);

        void clear();
    };

    std::vector<CountingSS> logs;
    std::ofstream os;

    void dump(CountingSS &sublog);
public:
    ReadLogger(size_t threads, const std::experimental::filesystem::path &out_file) : logs(threads), os() {os.open(out_file);}
    ~ReadLogger();

    ReadLogger(ReadLogger &&other)  = default;
    ReadLogger &operator=(ReadLogger &&other) = default;
    ReadLogger(const ReadLogger &other) = delete;
    ReadLogger &operator=(const ReadLogger &other) = delete;

    void flush();
    void logRead(AlignedRead &alignedRead);
    void logRerouting(AlignedRead &alignedRead, const dbg::GraphAlignment &initial, const dbg::GraphAlignment &corrected, const std::string &message);
    void logInvalidate(AlignedRead &alignedRead, const std::string &message) {
        CountingSS &ss = logs[omp_get_thread_num()];
        ss << alignedRead.id << " invalidated " << message << ")\n";
        ss << alignedRead.id << "    final    " << alignedRead.path.getAlignment().str(true) << "\n";
        if(ss.size() > 100000) {
            dump(ss);
        }
    }
};


class RecordStorage {
private:
    std::vector<AlignedRead> reads;
    std::unordered_map<const dbg::Vertex *, VertexRecord> data;
    ReadLogger *readLogger;
public:
    size_t min_len;
    size_t max_len;
    bool track_suffixes;
    bool track_cov;
    bool log_changes;

private:
    void processPath(const dbg::CompactPath &cpath, const std::function<void(dbg::Vertex &, const Sequence &)> &task,
                            const std::function<void(Segment<dbg::Edge>)> &edge_task = [](Segment<dbg::Edge>){}) const;
public:
    RecordStorage(dbg::SparseDBG &dbg, size_t _min_len, size_t _max_len, size_t threads,
                  ReadLogger &readLogger, bool _track_cov = false, bool log_changes = false, bool track_suffixes = true);

    typedef typename std::vector<AlignedRead>::iterator iterator;
    typedef typename std::vector<AlignedRead>::const_iterator const_iterator;

    RecordStorage &operator=(const RecordStorage &other) = delete;
    RecordStorage(const RecordStorage &other) = delete;
    RecordStorage &operator=(RecordStorage &&other) = default;
    RecordStorage(RecordStorage &&other) = default;

    const VertexRecord &getRecord(const dbg::Vertex &v) const;
    iterator begin() {return reads.begin();}
    iterator end() {return reads.end();}
    const_iterator begin() const {return reads.begin();}
    const_iterator end() const {return reads.end();}
    AlignedRead &operator[](size_t ind) {return reads[ind];}
    const AlignedRead &operator[](size_t ind) const {return reads[ind];}
    size_t getMinLen() const {return min_len;}
    size_t getMaxLen() const {return max_len;}
    bool isTrackingCov() const {return track_cov;}
    bool isTrackingSuffixes() const {return track_suffixes;}
    size_t size() const {return reads.size();}

    std::function<std::string(dbg::Edge &)> labeler() const;

    void addSubpath(const dbg::CompactPath &cpath);
    void removeSubpath(const dbg::CompactPath &cpath);
    void addRead(AlignedRead &&read);
    void delayedInvalidateRead(AlignedRead &read, const std::string &message);
    void reroute(AlignedRead &alignedRead, const dbg::GraphAlignment &initial, const dbg::GraphAlignment &corrected, const std::string &message);
    void reroute(AlignedRead &alignedRead, const dbg::GraphAlignment &corrected, const std::string &message);
    bool apply(AlignedRead &alignedRead);

    void delayedInvalidateBad(logging::Logger &logger, size_t threads, double threshold, const std::string &message);
    void delayedInvalidateBad(logging::Logger &logger, size_t threads, const std::function<bool(const dbg::Edge &)> &is_bad, const std::string &message);
    void invalidateSubreads(logging::Logger &logger, size_t threads);

    template<class I>
    void fill(I begin, I end, dbg::SparseDBG &dbg, size_t min_read_size, logging::Logger &logger, size_t threads);
    void trackSuffixes(logging::Logger &logger, size_t threads);
    void untrackSuffixes();

    //    void updateExtensionSize(logging::Logger &logger, size_t threads, size_t new_max_extension);
    void applyCorrections(logging::Logger &logger, size_t threads);
    void printReadAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const;
    void printReadFasta(logging::Logger &logger, const std::experimental::filesystem::path &path) const;
    void printFullAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const;
    ReadLogger &getLogger() {return *readLogger;}
    void flush() {readLogger->flush();}

    void Save(std::ostream &os) const;

    void Load(std::istream &is, dbg::SparseDBG &dbg);
};

void SaveAllReads(const std::experimental::filesystem::path &fname, const std::vector<RecordStorage *> &recs);

void LoadAllReads(const std::experimental::filesystem::path &fname, const std::vector<RecordStorage *> &recs, dbg::SparseDBG &dbg);

template<class I>
void RecordStorage::fill(I begin, I end, dbg::SparseDBG &dbg, size_t min_read_size, logging::Logger &logger, size_t threads) {
    if (track_cov) {
        logger.info() << "Cleaning edge coverages" << std::endl;
        for(dbg::Edge & edge: dbg.edges()) {
            edge.incCov(-edge.intCov());
        }
    }
    logger.info() << "Collecting alignments of sequences to the graph" << std::endl;
    if(track_suffixes) {
        logger.info() << "Storing suffixes of read paths of length up to " << this->max_len << std::endl;
    }
    ParallelRecordCollector<std::tuple<size_t, std::string, dbg::CompactPath>> tmpReads(threads);
    ParallelCounter cnt(threads);
    std::function<void(size_t, StringContig &)> read_task = [this, min_read_size, &tmpReads, &cnt, &dbg](size_t pos, StringContig & scontig) {
        Contig contig = scontig.makeContig();
        if(contig.size() < min_read_size) {
            tmpReads.emplace_back(pos, contig.id, dbg::CompactPath());
            return;
        }
        dbg::GraphAlignment path = dbg::GraphAligner(dbg).align(contig.seq);
        dbg::CompactPath cpath(path);
        dbg::GraphAlignment rcPath = path.RC();
        dbg::CompactPath crcPath(rcPath);
        addSubpath(cpath);
        addSubpath(crcPath);
        cnt += cpath.size();
        tmpReads.emplace_back(pos, contig.id, cpath);
    };
    processRecords(begin, end, logger, threads, read_task);
    reads.resize(tmpReads.size());
    for(auto &rec : tmpReads) {
        VERIFY(std::get<0>(rec) < reads.size());
        reads[std::get<0>(rec)] = AlignedRead(std::get<1>(rec), std::move(std::get<2>(rec)));
    }
    logger.info() << "Alignment collection finished. Total length of alignments is " << cnt.get() << std::endl;
}

//struct GraphError {
//    AlignedRead *read;
//    size_t from;
//    size_t to;
//    size_t _size;
//    dbg::GraphAlignment correction;
//
//    GraphError(AlignedRead *read, size_t from, size_t to, size_t size) :
//            read(read), from(from), to(to), _size(size) {}
//
//    size_t size() const {
//        return _size;
//    }
//
//    bool operator==(const GraphError &other) const {
//        return read == other.read && from == other.from && to == other.to;
//    }
//
//};

//class ErrorStorage {
//private:
//    std::vector<GraphError> errors;
//    RecordStorage &recs;
//public:
//
//    explicit ErrorStorage(RecordStorage &_recs): recs(_recs){
//    };
//
//    void fill(logging::Logger &logger, size_t threads, double error_threshold, double extend_threshold) {
//        logger.info() << "Correcting low covered regions in reads" << std::endl;
//        ParallelRecordCollector<GraphError> result(threads);
//#pragma omp parallel for default(none) shared(error_threshold, extend_threshold, result, logger)
//        for(size_t read_ind = 0; read_ind < recs.size(); read_ind++) { // NOLINT(modernize-loop-convert)
//            AlignedRead &alignedRead = recs[read_ind];
//            std::stringstream ss;
//            CompactPath &initial_cpath = alignedRead.path;
//            dbg::GraphAlignment path = initial_cpath.getAlignment();
//            dbg::GraphAlignment corrected_path(path.start());
//            for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
//                VERIFY_OMP(corrected_path.finish() == path.getVertex(path_pos));
//                dbg::Edge &edge = path[path_pos].contig();
//                if (edge.getCoverage() >= error_threshold) {
//                    continue;
//                }
//                size_t step_back = 0;
//                size_t step_front = 0;
//                size_t size = edge.size();
//                while (step_back < corrected_path.size() &&
//                       corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage() <
//                       extend_threshold) {
//                    size += corrected_path[corrected_path.size() - step_back - 1].size();
//                    step_back += 1;
//                }
//                while (step_front + path_pos + 1 < path.size() &&
//                       path[step_front + path_pos + 1].contig().getCoverage() < extend_threshold) {
//                    size += path[step_front + path_pos + 1].size();
//                    step_front += 1;
//                }
//                result.emplace_back(&alignedRead, corrected_path.size() - step_back, step_front + path_pos + 1, size);
//                path_pos += step_front;
//            }
//        }
//        for(GraphError &error : result)
//            errors.emplace_back(error);
//    }
//
//    void apply() {
//        std::function<bool(const GraphError &, const GraphError &)> compare =
//                [](const GraphError &a, const GraphError &b) {
//            if (a.read != b.read)
//                return a.read < b.read;
//            if(a.from != b.from)
//                return a.from < b.from;
//            return a.to < b.to;
//        };
//        __gnu_parallel::sort(errors.begin(), errors.end(), compare);
////#pragma omp parallel for default(none) shared(errors, recs)
//        for(size_t i = 0; i < errors.size(); i++) {
//            if(i > 0 && errors[i].read == errors[i - 1].read) {
//                continue;
//            }
//            std::vector<GraphError> corrected_errors;
//            AlignedRead &read = *errors[i].read;
//            for(size_t j = i; j < errors.size(); j++) {
//                if(errors[j].read != errors[i].read)
//                    break;
//                if(errors[j].correction.valid()) {
//                    corrected_errors.push_back(errors[j]);
//                }
//            }
//            if(corrected_errors.empty())
//                continue;
//            dbg::GraphAlignment initial = read.path.getAlignment();
//            std::vector<Segment<Edge>> corrected;
//            dbg::Vertex &start = corrected_errors[0].from == 0 ? corrected_errors[0].correction.start() : initial.start();
//            size_t prev = 0;
//            for(GraphError & error : corrected_errors) {
//                for(size_t j = prev; j < error.from; j++) {
//                    corrected.emplace_back(initial[j]);
//                }
//                for(size_t j = 0; j < error.size(); j++) {
//                    corrected.emplace_back(error.correction[j]);
//                }
//                prev = error.to;
//            }
//            for(size_t j = prev; j < initial.size(); j++) {
//                corrected.emplace_back(initial[j]);
//            }
//            dbg::GraphAlignment corrected_alignment(&start, std::move(corrected));
//            recs.reroute(read, initial, corrected_alignment, "unknown");
//        }
//    }
//
//    void sortByLength() {
//        std::function<bool(const GraphError &, const GraphError &)> compare =
//                [](const GraphError &a, const GraphError &b) {
//                    return a._size < b._size;
//                };
//        __gnu_parallel::sort(errors.begin(), errors.end(), compare);
//    }
//};
