#pragma once
#include "read_alignment_storage.hpp"
namespace ag {
    class ReadLogger : public AlignedReadStorageListener {
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
        ReadLogger(AlignedReadStorageFire &fire, size_t threads, const std::experimental::filesystem::path &out_file);
        virtual ~ReadLogger();
        ReadLogger(ReadLogger &&other) = default;
        ReadLogger &operator=(ReadLogger &&other) = default;
        ReadLogger(const ReadLogger &other) = delete;
        ReadLogger &operator=(const ReadLogger &other) = delete;

        void flush();
        void logRerouting(AlignedRead &alignedRead, const std::string &message);

        virtual void fireAddRead(const AlignedRead &read);
        void fireDelayedRerouteRead(AlignedRead &read, const std::string &message);;
        void logInvalidate(AlignedRead &alignedRead, const std::string &message);
        void fireDelayedInvalidateRead(AlignedRead &read, const std::string &message);;
        void fireAppliedCorrections(size_t cnt) override {flush();}
    };
}