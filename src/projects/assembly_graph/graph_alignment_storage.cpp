#include "graph_alignment_storage.hpp"

namespace ag {
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
        for (CountingSS &sublog: logs) {
            dump(sublog);
        }
    }

    ReadLogger::~ReadLogger() {
        flush();
        os.close();
    }
}