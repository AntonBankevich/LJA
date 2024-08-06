#pragma once
#include "compact_path.hpp"
#include "assembly_graph.hpp"
#include "fstream"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include <experimental/filesystem>

//TODO: get rid of read names entirely and replace with integer ids
inline std::string encodeReadId(const std::string &str) {
    std::string res = str;
    std::replace(res.begin(), res.end(), ' ', '^');
    std::replace(res.begin(), res.end(), '\t', '$');
    return std::move(res);
}

inline std::string decodeReadId(const std::string &str) {
    std::string res = str;
    std::replace(res.begin(), res.end(), '^', ' ');
    std::replace(res.begin(), res.end(), '$', '\t');
    return std::move(res);
}

namespace ag {
    template<class Traits>
    class AlignedRead {
    private:
        CompactPath<Traits> corrected_path;
        bool corrected = false;
    public:
        std::string id;
        CompactPath<Traits> path;

        AlignedRead() = default;

        AlignedRead(AlignedRead &&other) noexcept = default;

        AlignedRead &operator=(AlignedRead &&other) = default;

        explicit AlignedRead(std::string readId) : id(std::move(readId)), corrected(false) {}

        AlignedRead(std::string readId, const ag::GraphPath<Traits> &_path) : id(std::move(readId)), path(_path),
                                                                        corrected(false) {}

        AlignedRead(std::string readId, CompactPath<Traits> _path) : id(std::move(readId)), path(std::move(_path)),
                                                             corrected(false) {}

        bool operator<(const AlignedRead &other) const { return id < other.id; }

        void delayedInvalidate();

        bool checkCorrected() const { return corrected; }

        bool valid() const { return path.valid(); }

        void correct(CompactPath<Traits> &&cpath);

        void applyCorrection();

        static AlignedRead Load(std::istream &is, const IdIndex<typename Traits::Vertex> &index) {
            std::string id;
            is >> id;
            id = decodeReadId(id);
            return {id, CompactPath<Traits>::Load(is, index)};
        }
    };

    template<class Traits>
    inline std::ostream &operator<<(std::ostream &os, const AlignedRead<Traits> &alignedRead) {
        return os << encodeReadId(alignedRead.id) << " " << alignedRead.path;
    }

    template<class Traits>
    class RecordStorage;

    template<class Traits>
    struct VertexRecord {
        friend RecordStorage<Traits>;
    public:
        typedef typename Traits::Vertex Vertex;
    private:
        typedef std::vector<std::pair<Sequence, size_t>> Storage;
        typedef Storage::const_iterator const_iterator;
        Vertex &v;
        Storage paths;
        size_t zero_cnt = 0;
        size_t cov = 0;

        void lock() const { v.lock(); }

        void unlock() const { v.unlock(); }

        void addPath(const Sequence &seq);

        void removePath(const Sequence &seq);

        void clear() { paths.clear(); }

    public:
        explicit VertexRecord(Vertex &_v) : v(_v) {}

        VertexRecord(const VertexRecord &) = delete;

        VertexRecord(VertexRecord &&other) noexcept: v(other.v), paths(std::move(other.paths)),
                                                     zero_cnt(other.zero_cnt), cov(other.cov) {}

        VertexRecord &operator=(const VertexRecord &) = delete;

        size_t coverage() const { return cov; }

        std::string str() const;

        const_iterator begin() const { return paths.begin(); }

        const_iterator end() const { return paths.end(); }

        size_t countStartsWith(const Sequence &seq) const;

        bool isDisconnected(const typename Traits::Edge &edge) const;

        std::vector<ag::GraphPath<Traits>> getBulgeAlternatives(const Vertex &end, double threshold) const;

        std::vector<ag::GraphPath<Traits>> getTipAlternatives(size_t len, double threshold) const;

        unsigned char getUniqueExtension(const Sequence &start, size_t min_good_cov, size_t max_bad_cov) const;

        ag::CompactPath<Traits> getFullUniqueExtension(const Sequence &start, size_t min_good_cov, size_t max_bad_cov,
                                                       size_t max_size = size_t(-1)) const;
    };

    template<class Traits>
    inline std::ostream &operator<<(std::ostream &os, const VertexRecord<Traits> &rec) { return os << rec.str(); }

    class ReadLogger {
    private:
        class CountingSS {
        private:
            std::stringstream log;
            size_t len;
        public:
            CountingSS() : log(), len(0) {}

            std::string str() { return log.str(); }

            size_t size() const { return len; }

            CountingSS &operator<<(const std::string &s);

            CountingSS &operator<<(const size_t &s);

            void clear();
        };

        std::vector<CountingSS> logs;
        std::ofstream os;

        void dump(CountingSS &sublog);

    public:
        ReadLogger(size_t threads, const std::experimental::filesystem::path &out_file) : logs(threads), os() {
            os.open(out_file);
        }

        ~ReadLogger();

        ReadLogger(ReadLogger &&other) = default;

        ReadLogger &operator=(ReadLogger &&other) = default;

        ReadLogger(const ReadLogger &other) = delete;

        ReadLogger &operator=(const ReadLogger &other) = delete;

        void flush();

        template<class Traits>
        void logRead(AlignedRead<Traits> &alignedRead) {
            CountingSS &ss = logs[omp_get_thread_num()];
            ss << alignedRead.id << " initial " << alignedRead.path.unpack().covStr(true) << "\n";
            if (ss.size() > 100000) {
                dump(ss);
            }
        }

        template<class Traits>
        void logRerouting(AlignedRead<Traits> &alignedRead, const ag::GraphPath<Traits> &initial,
                          const ag::GraphPath<Traits> &corrected, const std::string &message) {
            CountingSS &ss = logs[omp_get_thread_num()];
            size_t left = 0;
            size_t right = 0;
            size_t left_len = 0;
            size_t right_len = 0;
            while (left < corrected.size() && left < initial.size()) {
                if (initial[left] != corrected[left])
                    break;
                left_len += initial[left].size();
                left++;
            }
            while (left + right < corrected.size() && left + right < initial.size()) {
                if (initial[initial.size() - right - 1] != corrected[corrected.size() - right - 1])
                    break;
                right_len += initial[initial.size() - right - 1].size();
                right++;
            }
            ss << alignedRead.id << " " << message << " " << left << "(" << left_len << ") " << right << "("
               << right_len << ")\n";
            ss << alignedRead.id << "  initial  " << initial.subPath(left, initial.size() - right).covStr(true) << "\n";
            ss << alignedRead.id << " corrected " << corrected.subPath(left, corrected.size() - right).covStr(true)
               << "\n";
//        ss << alignedRead.id << " rc  initial  " << initial.RC().str(true) << "\n";
//        ss << alignedRead.id << " rc corrected " << corrected.RC().str(true) << "\n";
            if (ss.size() > 100000) {
                dump(ss);
            }
        }

        template<class Traits>
        void logInvalidate(AlignedRead<Traits> &alignedRead, const std::string &message) {
            CountingSS &ss = logs[omp_get_thread_num()];
            ss << alignedRead.id << " invalidated " << message << ")\n";
            ss << alignedRead.id << "    final    " << alignedRead.path.unpack().covStr(true) << "\n";
            if (ss.size() > 100000) {
                dump(ss);
            }
        }
    };

    template<class Traits>
    class RecordStorage {
    protected:
        std::vector<AlignedRead<Traits>> reads;
        std::unordered_map<typename Traits::Vertex::ConstVertexId, VertexRecord<Traits>> data;
        ReadLogger *readLogger;
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge::EdgeId EdgeId;
        typedef typename Traits::Vertex::VertexId VertexId;

        size_t min_len;
        size_t max_len;
        bool track_suffixes;
        bool track_cov;
        bool log_changes;

    private:
        void
        processPath(const ag::CompactPath<Traits> &cpath, const std::function<void(Vertex &, const Sequence &)> &task,
                    const std::function<void(Segment<Edge>)> &edge_task = [](Segment<Edge>) {}) const;

    public:
        RecordStorage(AssemblyGraph<Traits> &graph, size_t _min_len, size_t _max_len,
                      bool _track_cov = false, bool log_changes = false,
                      bool track_suffixes = true) :
                min_len(_min_len), max_len(_max_len), track_cov(_track_cov), readLogger(nullptr),
                log_changes(log_changes), track_suffixes(track_suffixes) {
            for(Vertex &v : graph.vertices()) {
                data.emplace(std::piecewise_construct,
                             std::forward_as_tuple(v.getId()),
                             std::forward_as_tuple(v));
            }
        }

        void setReadLogger(ReadLogger &_readLogger) {
            readLogger = &_readLogger;
        }

        typedef typename std::vector<AlignedRead<Traits>>::iterator iterator;
        typedef typename std::vector<AlignedRead<Traits>>::const_iterator const_iterator;

        RecordStorage &operator=(const RecordStorage &other) = delete;

        RecordStorage(const RecordStorage &other) = delete;

        RecordStorage &operator=(RecordStorage &&other)  noexcept = default;

        RecordStorage(RecordStorage &&other) = default;

        const VertexRecord<Traits> &getRecord(const Vertex &v) const;

        iterator begin() { return reads.begin(); }

        iterator end() { return reads.end(); }

        const_iterator begin() const { return reads.begin(); }

        const_iterator end() const { return reads.end(); }

        AlignedRead<Traits> &operator[](size_t ind) { return reads[ind]; }

        const AlignedRead<Traits> &operator[](size_t ind) const { return reads[ind]; }

        size_t getMinLen() const { return min_len; }

        size_t getMaxLen() const { return max_len; }

        bool isTrackingCov() const { return track_cov; }

        bool isTrackingSuffixes() const { return track_suffixes; }

        size_t size() const { return reads.size(); }

        std::function<std::string(const Edge &)> labeler() const;

        void addSubpath(const CompactPath<Traits> &cpath);

        void removeSubpath(const CompactPath<Traits> &cpath);

        AlignedRead<Traits>& addRead(std::string readId, ag::GraphPath<Traits> alignment = {});

        AlignedRead<Traits>& addRead(AlignedRead<Traits>&& ar);

        void delayedInvalidateRead(AlignedRead<Traits> &read, const std::string &message);

        void reroute(AlignedRead<Traits> &alignedRead, const GraphPath<Traits> &initial, const GraphPath<Traits> &corrected,
                     const std::string &message);

        void reroute(AlignedRead<Traits> &alignedRead, const GraphPath<Traits> &corrected, const std::string &message);

        bool apply(AlignedRead<Traits> &alignedRead);

        void
        delayedInvalidateBad(logging::Logger &logger, size_t threads, double threshold, const std::string &message);

        void
        delayedInvalidateBad(logging::Logger &logger, size_t threads, const std::function<bool(const Edge &)> &is_bad,
                             const std::string &message);

        void invalidateSubreads(logging::Logger &logger, size_t threads);

        void trackSuffixes(logging::Logger &logger, size_t threads);

        void untrackSuffixes();

        //    void updateExtensionSize(logging::Logger &logger, size_t threads, size_t new_max_extension);
        void applyCorrections(logging::Logger &logger, size_t threads);

        void printReadAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const;

        void printReadFasta(logging::Logger &logger, const std::experimental::filesystem::path &path) const;

        void printReadPaths(logging::Logger &logger, const std::experimental::filesystem::path &aln_path,
                                                     const std::experimental::filesystem::path &gfa_path,
                                                     const std::experimental::filesystem::path &rp_path,
                                                     size_t k) const;

        void printFullAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const;

        ReadLogger &getLogger() { return *readLogger; }

        void flush() { if(readLogger != nullptr) readLogger->flush(); }

        void Save(std::ostream &os) const;

        void Load(std::istream &is, const IdIndex<Vertex> &index);
    };

    template<class Traits>
    void
    SaveAllReads(const std::experimental::filesystem::path &fname, const std::vector<RecordStorage<Traits> *> &recs);

    template<class Traits>
    void LoadAllReads(const std::experimental::filesystem::path &fname, const std::vector<RecordStorage<Traits> *> &recs,
                      const IdIndex<typename Traits::Vertex> &index);


    template<class Traits>
    void AlignedRead<Traits>::correct(CompactPath<Traits> &&cpath) {
        VERIFY_MSG(!corrected, "Attempt to correct path while previous correction was not yet applied");
        corrected_path = std::move(cpath);
        corrected = true;
    }

    template<class Traits>
    void AlignedRead<Traits>::applyCorrection() {
        if (corrected)
            path = std::move(corrected_path);
        corrected_path = {};
        corrected = false;
    }

    template<class Traits>
    void AlignedRead<Traits>::delayedInvalidate() {
        VERIFY(!corrected);
        corrected_path = {};
        corrected = true;
    }

    template<class Traits>
    void VertexRecord<Traits>::addPath(const Sequence &seq) {
        lock();
        cov += 1;
        for (std::pair<Sequence, size_t> &path: paths) {
            if (path.first == seq) {
                path.second += 1;
                if (path.second == 1)
                    zero_cnt -= 1;
                unlock();
                return;
            }
        }
        paths.emplace_back(seq, 1);
        unlock();
    }

    template<class Traits>
    void VertexRecord<Traits>::removePath(const Sequence &seq) {
        lock();
        bool found = false;
        for (std::pair<Sequence, size_t> &path: paths) {
            if (path.first == seq) {
                found = true;
                VERIFY(path.second > 0);
                path.second -= 1;
                cov -= 1;
                if (path.second == 0)
                    zero_cnt += 1;
                break;
            }
        }
        if (!found) {
            std::cout << "Error" << std::endl;
            unlock();
            std::cout << seq << std::endl;
            std::cout << this->str() << std::endl;
        }
        VERIFY(found);
        if (zero_cnt > paths.size() / 3) {
            std::vector<std::pair<Sequence, size_t>> new_paths;
            for (std::pair<Sequence, size_t> &rec: paths) {
                if (rec.second != 0) {
                    new_paths.emplace_back(std::move(rec.first), rec.second);
                }
            }
            std::swap(paths, new_paths);
            zero_cnt = 0;
        }
        unlock();
    }

    template<class Traits>
    bool VertexRecord<Traits>::isDisconnected(const typename Traits::Edge &edge) const {
        if (edge.getFinish().outDeg() == 0)
            return false;
        for (const std::pair<Sequence, size_t> &rec: paths) { // NOLINT(readability-use-anyofallof)
            if (rec.second > 0 && rec.first.startsWith(edge.nuclLabel()) && rec.first.size() > 1)
                return false;
        }
        return true;
    }

    template<class Traits>
    size_t VertexRecord<Traits>::countStartsWith(const Sequence &seq) const {
//        lock();
        size_t cnt = 0;
        for (const std::pair<Sequence, size_t> &rec: paths) {
            if (rec.first.startsWith(seq)) {
                cnt += rec.second;
            }
        }
//        unlock();
        return cnt;
    }

    template<class Traits>
    std::vector<GraphPath<Traits>> VertexRecord<Traits>::getBulgeAlternatives(const Vertex &end, double threshold) const {
//        lock();
        std::vector<std::pair<Sequence, size_t>> candidates;
        for (const auto &extension: paths) {
            GraphPath<Traits> unpacked = CompactPath<Traits>(v, extension.first).unpack();
            for (size_t i = 1; i <= unpacked.size(); i++) {
                if (end == unpacked.getVertex(i)) {
                    candidates.emplace_back(extension.first.Subseq(0, i), extension.second);
                }
            }
        }
//        unlock();
        if (candidates.empty())
            return {};
        std::sort(candidates.begin(), candidates.end());
        std::vector<GraphPath<Traits>> res;
        size_t cnt = 0;
        for (size_t i = 0; i < candidates.size(); i++) {
            if (i > 0 && candidates[i - 1].first != candidates[i].first) {
                if (cnt >= threshold)
                    res.emplace_back(CompactPath<Traits>(v, candidates[i - 1].first).unpack());
                cnt = 0;
            }
            cnt += candidates[i].second;
        }
        if (cnt >= threshold)
            res.emplace_back(CompactPath<Traits>(v, candidates.back().first).unpack());
        return std::move(res);
    }

    template<class Traits>
    CompactPath <Traits>
    VertexRecord<Traits>::getFullUniqueExtension(const Sequence &start, size_t min_good_cov, size_t max_bad_cov,
                                                 size_t max_size) const {
        Sequence res = start;
        while (res.size() < max_size) {
            unsigned char next = getUniqueExtension(res, min_good_cov, max_bad_cov);
            if (next == (unsigned char) -1)
                break;
            std::vector<unsigned char> ext = {next};
            res = res + Sequence(ext);
        }
        return {v, res};
    }

    template<class Traits>
    unsigned char
    VertexRecord<Traits>::getUniqueExtension(const Sequence &start, size_t min_good, size_t max_bad) const {
        std::vector<size_t> counts(4);
        for (const auto &extension: paths) {
            if (extension.first.size() > start.size() && extension.first.startsWith(start)) {
                counts[extension.first[start.size()]] += extension.second;
            }
        }
        size_t bad = 0;
        size_t good = 0;
        size_t res = 0;
        for (size_t c = 0; c < 4; c++) {
            if (counts[c] <= max_bad)
                bad++;
            if (counts[c] >= min_good) {
                good++;
                res = c;
            }
        }
        if (bad != 3 || good != 1)
            return (unsigned char) (-1);
        return (unsigned char) (res);
    }

    template<class Traits>
    std::vector<GraphPath<Traits>> VertexRecord<Traits>::getTipAlternatives(size_t len, double threshold) const {
        len += std::max<size_t>(30, len / 20);
//        lock();
        std::vector<std::pair<Sequence, size_t>> candidates;
        for (const auto &extension: paths) {
            GraphPath<Traits> unpacked = CompactPath<Traits>(v, extension.first).unpack();
            if (unpacked.truncLen() >= len) {
                unpacked.cutBack(unpacked.truncLen() - len);
                candidates.emplace_back(CompactPath<Traits>(unpacked).cpath(), extension.second);
            }
        }
//        unlock();
        if (candidates.empty())
            return {};
        std::sort(candidates.begin(), candidates.end());
        std::vector<GraphPath<Traits>> res;
        size_t cnt = 0;
        for (size_t i = 0; i < candidates.size(); i++) {
            if (i > 0 && candidates[i - 1].first != candidates[i].first) {
                if (cnt > threshold) {
                    GraphPath<Traits> cp = CompactPath<Traits>(v, candidates[i - 1].first).unpack();
                    cp.cutBack(cp.truncLen() - len);
                    res.emplace_back(cp);
                }
                cnt = 0;
            }
            cnt += candidates[i].second;
        }
        if (cnt > threshold) {
            GraphPath<Traits> cp = CompactPath<Traits>(v, candidates.back().first).unpack();
            cp.cutBack(cp.truncLen() - len);
            res.emplace_back(cp);
        }
        return std::move(res);
    }

    template<class Traits>
    std::string VertexRecord<Traits>::str() const {
        std::stringstream ss;
        lock();
        for (const auto &path: paths) {
            ss << path.first << " " << path.second << std::endl;
        }
        unlock();
        return ss.str();
    }


    template<class Traits>
    std::function<std::string(const typename Traits::Edge &)> RecordStorage<Traits>::labeler() const {
        if (track_suffixes)
            return [this](const Edge &edge) {
                const VertexRecord<Traits> &rec = getRecord(edge.getStart());
                std::stringstream ss;
                size_t cnt = 0;
                for (const auto &ext: rec) {
                    if (ext.first.startsWith(edge.nuclLabel())) {
                        if (cnt < 30)
                            ss << ext.first << "(" << ext.second << ")\\n";
                        cnt++;
                    }
                }
                if (cnt > 30) {
                    ss << "and another " << (cnt - 30) << " records\\n";
                }
                return ss.str();
            };
        else
            return [](const Edge &edge) {
                return std::string();
            };
    }

    template<class Traits>
    AlignedRead<Traits>& RecordStorage<Traits>::addRead(std::string readId, ag::GraphPath<Traits> alignment) {
        reads.emplace_back(std::move(readId), std::move(alignment));
        addSubpath(reads.back().path);
        addSubpath(reads.back().path.RC());
        return reads.back();
    }

    template<class Traits>
    AlignedRead<Traits>& RecordStorage<Traits>::addRead(AlignedRead<Traits> &&ar) {
        reads.emplace_back(std::move(ar));
        addSubpath(reads.back().path);
        addSubpath(reads.back().path.RC());
        return reads.back();
    }

    template<class Traits>
    void RecordStorage<Traits>::delayedInvalidateRead(AlignedRead<Traits> &read,
                                                      const std::string &message) { // NOLINT(readability-convert-member-functions-to-static)
        if (log_changes && readLogger != nullptr)
            readLogger->logInvalidate(read, message);
        read.delayedInvalidate();
    }

    template<class Traits>
    void RecordStorage<Traits>::delayedInvalidateBad(logging::Logger &logger, size_t threads, double threshold,
                                                     const std::string &message) {
        const std::function<bool(const Edge &)> &is_bad = [threshold](const Edge &edge) {
            return edge.getCoverage() < threshold;
        };
        delayedInvalidateBad(logger, threads, is_bad, message);
    }

    template<class Traits>
    void
    RecordStorage<Traits>::delayedInvalidateBad(logging::Logger &logger, size_t threads, const std::function<bool(
            const Edge &)> &is_bad, const std::string &message) {
        ParallelCounter cnt(threads);
        omp_set_num_threads(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(cnt, message, is_bad)
        for (size_t i = 0; i < reads.size(); i++) {
            AlignedRead<Traits> &alignedRead = reads[i];
            if (!alignedRead.valid())
                continue;
            GraphPath<Traits>al = alignedRead.path.unpack();
            size_t l = 0;
            size_t r = al.size();
            while (l < al.size() && is_bad(al[l].contig())) {
                l++;
            }
            while (r > 0 && is_bad(al[r - 1].contig())) {
                r--;
            }
            bool middle_bad = false;
            size_t len = 0;
            for (size_t j = l; j < r; j++) {
                if (is_bad(al[j].contig())) {
                    middle_bad = true;
                    break;
                }
                len += al[j].size();
            }
            if (middle_bad || l >= r) {
                delayedInvalidateRead(alignedRead, message);
                cnt += 1;
            } else if (r - l != al.size()) {
                GraphPath<Traits>sub = al.subPath(l, r);
                if (sub.truncLen() < 1000)
                    delayedInvalidateRead(alignedRead, message);
                else {
                    std::string extra_message =
                            message + "_EndsClipped_" + itos(al.subPath(0, l).truncLen()) + "_" + itos(
                                    al.subPath(r, al.size()).truncLen());
                    reroute(alignedRead, al, al.subPath(l, r), extra_message);
                }
                cnt += 1;
            }
        }
        if (cnt.get() > 0)
            logger.info() << "Could not correct " << cnt.get() << " reads. They were removed or truncated."
                          << std::endl;
    }

    template<class Traits>
    void RecordStorage<Traits>::invalidateSubreads(logging::Logger &logger, size_t threads) {
        for (AlignedRead<Traits> &alignedRead: reads) {
            const VertexRecord<Traits> &rec = this->getRecord(alignedRead.path.start());
            size_t cnt = rec.countStartsWith(alignedRead.path.cpath());
            VERIFY_MSG(cnt >= 1, "This function assumes that suffixes are stored for complete reads")
            if (rec.countStartsWith(alignedRead.path.cpath()) >= 2) {
                delayedInvalidateRead(alignedRead, "Subread");
            }
        }
        logger.info() << "Uncorrected reads were removed." << std::endl;
    }

    template<class Traits>
    void RecordStorage<Traits>::addSubpath(const CompactPath <Traits> &cpath) {
        if (!cpath.valid())
            return;
        std::function< void(Vertex
                            &,
                            const Sequence &)> vertex_task = [](Vertex &v, const Sequence &s) {};
        if (track_suffixes)
            vertex_task = [this](Vertex &v, const Sequence &s) {
                data.find(v.getId())->second.addPath(s);
            };
        std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg) {};
        if (track_cov)
            edge_task = [](Segment<Edge> seg) {
                seg.contig().incCov(seg.size());
            };
        processPath(cpath, vertex_task, edge_task);
    }

    template<class Traits>
    void RecordStorage<Traits>::removeSubpath(const CompactPath <Traits> &cpath) {
        if (!cpath.valid())
            return;
        std::function< void(Vertex
                            &,
                            const Sequence &)> vertex_task = [](Vertex &v, const Sequence &s) {};
        if (track_suffixes)
            vertex_task = [this](Vertex &v, const Sequence &s) {
                data.find(v.getId())->second.removePath(s);
            };
        std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg) {};
        if (track_cov)
            edge_task = [](Segment<Edge> seg) {
                seg.contig().incCov(size_t(-seg.size()));
            };
        processPath(cpath, vertex_task, edge_task);
    }

    template<class Traits>
    bool RecordStorage<Traits>::apply(AlignedRead<Traits> &alignedRead) {
        if (!alignedRead.checkCorrected())
            return false;
        GraphPath<Traits>path = alignedRead.path.unpack();
        this->removeSubpath(alignedRead.path);
        this->removeSubpath(alignedRead.path.RC());
        alignedRead.applyCorrection();
        this->addSubpath(alignedRead.path);
        this->addSubpath(alignedRead.path.RC());
        return true;
    }

    template<class Traits>
    void RecordStorage<Traits>::reroute(AlignedRead<Traits> &alignedRead, const GraphPath<Traits> &initial,
                                        const GraphPath<Traits> &corrected,
                                        const string &message) {
        if (log_changes && readLogger != nullptr) {
            readLogger->logRerouting(alignedRead, initial, corrected, message);
        } if (corrected.truncLen() < 500) {
            delayedInvalidateRead(alignedRead, "Deleted_after_" + message);
        } else {
            alignedRead.correct(CompactPath<Traits>(corrected));
        }
    }

    template<class Traits>
    void RecordStorage<Traits>::reroute(AlignedRead<Traits> &alignedRead, const GraphPath<Traits> &corrected,
                                        const string &message) {
        reroute(alignedRead, alignedRead.path.unpack(), corrected, message);
    }

    template<class Traits>
    void RecordStorage<Traits>::applyCorrections(logging::Logger &logger, size_t threads) {
        if (size() > 10000)
            logger.info() << "Applying corrections to reads" << std::endl;
        omp_set_num_threads(threads);
        ParallelCounter cnt(threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(cnt)
        for (size_t i = 0; i < reads.size(); i++) { // NOLINT(modernize-loop-convert)
            if (apply(reads[i]))
                cnt += 1;
        }
        flush();
        if (size() > 10000)
            logger.info() << "Applied correction to " << cnt.get() << " reads" << std::endl;
    }

    template<class Traits>
    void RecordStorage<Traits>::printReadAlignments(logging::Logger &logger,
                                                    const std::experimental::filesystem::path &path) const {
        logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
        std::string acgt = "ACGT";
        std::ofstream os;
        os.open(path);
        for (const AlignedRead<Traits> &read: reads) {
            const CompactPath<Traits> &al = read.path;
            if (!al.valid())
                continue;
            os << read.id << " " << al.start().getInnerId()
               << " " << al.cpath().str() << "\n";
            CompactPath<Traits> rc_al = al.RC();
            os << "-" << read.id << " " << rc_al.start().getInnerId() << " " << rc_al.cpath().str() << "\n";
        }
        os.close();
    }

    template<class Traits>
    void RecordStorage<Traits>::printReadFasta(logging::Logger &logger,
                                               const std::experimental::filesystem::path &path) const {
        logger.info() << "Printing reads to fasta file " << path << std::endl;
        std::string acgt = "ACGT";
        std::ofstream os;
        os.open(path);
        for (const AlignedRead<Traits> &read: reads) {
            const CompactPath<Traits> &al = read.path;
            if (!al.valid())
                continue;
            os << ">" << read.id << "\n" << read.path.unpack().Seq() << "\n";
        }
        os.close();
    }

    template<class Traits>
    void RecordStorage<Traits>::printReadPaths(logging::Logger &logger,
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

    template<class Traits>
    void RecordStorage<Traits>::printFullAlignments(logging::Logger &logger,
                                                    const std::experimental::filesystem::path &path) const {
        logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
        std::string acgt = "ACGT";
        std::ofstream os;
        os.open(path);
        for (const AlignedRead<Traits> &read: reads) {
            const CompactPath<Traits> &al = read.path;
            if (!al.valid())
                continue;
            os << read.id << " " << read.path.unpack().covStr(true) << "\n";
            os << "-" << read.id << " " << read.path.unpack().RC().covStr(true) << "\n";
        }
        os.close();
    }

//void RecordStorage<Traits>::updateExtensionSize(logging::Logger &logger, size_t threads, size_t new_max_extension) {
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

    template<class Traits>
    const VertexRecord<Traits> &RecordStorage<Traits>::getRecord(const Vertex &v) const {
        VERIFY(track_suffixes);
        auto it = data.find(v.getId());
        VERIFY(it != data.end());
        return it->second;
    }

    template<class Traits>
    void RecordStorage<Traits>::trackSuffixes(logging::Logger &logger, size_t threads) {
        logger.info() << "Collecting and storing read suffixes" << std::endl;
        VERIFY(!track_suffixes);
        track_suffixes = true;
        omp_set_num_threads(threads);
        std::function< void(Vertex&, const Sequence &)> vertex_task = [this](Vertex &v, const Sequence &s) {
            data.find(v.getId())->second.addPath(s);
        };
        std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg) {};
#pragma omp parallel for default(none) shared(vertex_task, edge_task)
        for (size_t i = 0; i < reads.size(); i++) {
            if (reads[i].valid()) {
                processPath(reads[i].path, vertex_task, edge_task);
                processPath(reads[i].path.RC(), vertex_task, edge_task);
            }
        }
    }

    template<class Traits>
    void RecordStorage<Traits>::untrackSuffixes() {
        if (track_suffixes) {
            track_suffixes = false;
            for (auto &it: this->data) {
                it.second.clear();
            }
        }
    }

    template<class Traits>
    void RecordStorage<Traits>::Save(std::ostream &os) const {
        os << size() << "\n";
        for (const AlignedRead<Traits> &alignedRead: *this) {
            os << alignedRead << "\n";
        }
    }

    template<class Traits>
    void RecordStorage<Traits>::Load(std::istream &is, const IdIndex<Vertex> &index) {
        size_t sz;
        is >> sz;
        for (size_t i = 0; i < sz; i++) {
            addRead(AlignedRead<Traits>::Load(is, index));
        }
    }

    template<class Traits>
    void SaveAllReads(const std::experimental::filesystem::path &fname,
                      const std::vector<RecordStorage<Traits> *> &recs) {
        std::ofstream os;
        os.open(fname);
        os << recs.size() << "\n";
        for (RecordStorage<Traits> *rs: recs) {
            rs->Save(os);
        }
        os.close();
    }

    template<class Traits>
    void
    LoadAllReads(const std::experimental::filesystem::path &fname, const std::vector<RecordStorage<Traits> *> &recs,
                 const IdIndex<typename Traits::Vertex> &index) {
        std::ifstream is;
        is.open(fname);
        size_t sz;
        is >> sz;
        VERIFY(sz == recs.size())
        for (RecordStorage<Traits> *recordStorage: recs) {
            recordStorage->Load(is, index);
        }
        is.close();
    }

    template<class Traits>
    void RecordStorage<Traits>::processPath(const CompactPath <Traits> &cpath, const std::function<void(Vertex &,const Sequence &)> &task,
                                                        const std::function<void(Segment<Edge>)> &edge_task) const {
        GraphPath<Traits>al = cpath.unpack();
        for (Segment<Edge> seg: al) {
            edge_task(seg);
        }
        size_t j = 1;
        size_t clen = al[0].contig().truncSize();
        for (size_t i = 1; i <= al.size(); i++) {
            clen -= al[i - 1].contig().truncSize();
            VERIFY(j >= i);
            while (j < al.size() && (clen < max_len || clen == 0)) {
                clen += al[j].contig().truncSize();
                j++;
            }
            if (clen >= min_len) {
                task(al.getVertex(i - 1), cpath.cpath().Subseq(i - 1, j));
            }
        }
    }

}
