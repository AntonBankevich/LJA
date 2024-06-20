#pragma once
#include "read_alignment_storage.hpp"
namespace ag {
    template<class Traits>
    class ReadLogger : public AlignedReadStorageListener<Traits> {
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
        ReadLogger(AlignedReadStorageFire<Traits> &fire, size_t threads, const std::experimental::filesystem::path &out_file) :
                            AlignedReadStorageListener<Traits>(fire, "ReadLogger"), logs(threads), os() {
            os.open(out_file);
        }

        ~ReadLogger();
        ReadLogger(ReadLogger &&other) = default;

        ReadLogger &operator=(ReadLogger &&other) = default;

        ReadLogger(const ReadLogger &other) = delete;

        ReadLogger &operator=(const ReadLogger &other) = delete;

        void flush();

        void logRerouting(AlignedRead<Traits> &alignedRead, const std::string &message) {
            const GraphPath<Traits> &initial = alignedRead.getPath();
            const GraphPath<Traits> &corrected = alignedRead.getCorrected();
            CountingSS &ss = logs[omp_get_thread_num()];
            size_t left = 0;
            PathPosition<Traits> left_pos_initial = initial.firstPosition();//left position
            PathPosition<Traits> left_pos_corrected = corrected.firstPosition();//left position
            size_t right = 0;
            PathPosition<Traits> right_pos_initial = initial.lastPosition();//size - right position
            PathPosition<Traits> right_pos_corrected = corrected.lastPosition();//size - right position
            size_t left_len = 0;
            size_t right_len = 0;
            while (left_pos_initial != initial.lastPosition() && left_pos_corrected != corrected.lastPosition()) {
                if (left_pos_initial.nextEdge() != left_pos_corrected.nextEdge())
                    break;
                left_len += left_pos_initial.nextEdge().rc().truncSize();
                left++;
                ++left_pos_initial;
                ++left_pos_corrected;
            }
            if(left_len > 0)
                left_len -= initial.leftCut();
            while (left_pos_initial != right_pos_initial && left_pos_corrected != right_pos_corrected) {
                if (right_pos_initial.prevEdge() != right_pos_corrected.prevEdge())
                    break;
                right_len += right_pos_initial.prevEdge().truncSize();
                right++;
                --right_pos_initial;
                --right_pos_corrected;
            }
            if(right_len > 0)
                right_len -= initial.rightCut();
            ss << alignedRead.getId() << " " << message << " " << left << "(" << left_len << ") " << right << "("
               << right_len << ")\n";
            ss << alignedRead.getId() << "  initial  " << initial.subPath(left_pos_initial, right_pos_initial).str() << "\n";
            ss << alignedRead.getId() << " corrected " << corrected.subPath(left_pos_corrected, right_pos_corrected).str()
               << "\n";
//        ss << alignedRead.getId() << " rc  initial  " << initial.RC().str(true) << "\n";
//        ss << alignedRead.getId() << " rc corrected " << corrected.RC().str(true) << "\n";
            if (ss.size() > 100000) {
                dump(ss);
            }
        }

        virtual void fireAddRead(const AlignedRead<Traits> &read) {
            CountingSS &ss = logs[omp_get_thread_num()];
            ss << read.getId() << "  new  " << read.getPath().str() << "\n";
            if (ss.size() > 100000) {
                dump(ss);
            }
        }

        void fireDelayedRerouteRead(AlignedRead<Traits> &read, const std::string &message) {
            logRerouting(read, message);
        };

        void logInvalidate(AlignedRead<Traits> &alignedRead, const std::string &message) {
            CountingSS &ss = logs[omp_get_thread_num()];
            ss << alignedRead.getId() << " invalidated " << message << ")\n";
            ss << alignedRead.getId() << "    final    " << alignedRead.path.unpack().str() << "\n";
            if (ss.size() > 100000) {
                dump(ss);
            }
        }
        void fireDelayedInvalidateRead(AlignedRead<Traits> &read, const std::string &message) {
            logRerouting(read, message);
        };

        void fireAppliedCorrections(size_t cnt) override {
            flush();
        }
    };

    template<class Traits>
    typename ReadLogger<Traits>::CountingSS &ReadLogger<Traits>::CountingSS::operator<<(const string &s) {
        log << s;
        len += s.size();
        return *this;
    }

    template<class Traits>
    typename ReadLogger<Traits>::CountingSS &ReadLogger<Traits>::CountingSS::operator<<(const size_t &s) {
        log << s;
        len += 5;
        return *this;
    }

    template<class Traits>
    void ReadLogger<Traits>::CountingSS::clear() {
        log = std::stringstream();
        len = 0;
    }

    template<class Traits>
    void ReadLogger<Traits>::dump(ReadLogger::CountingSS &sublog) {
#pragma omp critical
        {
            os << sublog.str();
        };
        sublog.clear();
    }

    template<class Traits>
    void ReadLogger<Traits>::flush() {
        for (CountingSS &sublog: logs) {
            dump(sublog);
        }
    }

    template<class Traits>
    ReadLogger<Traits>::~ReadLogger() {
        flush();
        os.close();
    }
}