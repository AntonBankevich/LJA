#include "read_logger.hpp"
using namespace ag;

typename ReadLogger::CountingSS &ReadLogger::CountingSS::operator<<(const string &s) {
    log << s;
    len += s.size();
    return *this;
}

typename ReadLogger::CountingSS &ReadLogger::CountingSS::operator<<(const size_t &s) {
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

void ReadLogger::logRerouting(AlignedRead &alignedRead, const string &message) {
    const GraphPath &initial = alignedRead.getPath();
    const GraphPath &corrected = alignedRead.getCorrected();
    CountingSS &ss = logs[omp_get_thread_num()];
    size_t left = 0;
    PathPosition left_pos_initial = initial.firstPosition();//left position
    PathPosition left_pos_corrected = corrected.firstPosition();//left position
    size_t right = 0;
    PathPosition right_pos_initial = initial.lastPosition();//size - right position
    PathPosition right_pos_corrected = corrected.lastPosition();//size - right position
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

void ReadLogger::fireAddRead(const AlignedRead &read) {
    CountingSS &ss = logs[omp_get_thread_num()];
    ss << read.getId() << "  new  " << read.getPath().str() << "\n";
    if (ss.size() > 100000) {
        dump(ss);
    }
}

void ReadLogger::fireDelayedRerouteRead(AlignedRead &read, const string &message) {
    logRerouting(read, message);
}

void ReadLogger::logInvalidate(AlignedRead &alignedRead, const string &message) {
    CountingSS &ss = logs[omp_get_thread_num()];
    ss << alignedRead.getId() << " invalidated " << message << ")\n";
    ss << alignedRead.getId() << "    final    " << alignedRead.getPath().str() << "\n";
    if (ss.size() > 100000) {
        dump(ss);
    }
}

void ReadLogger::fireDelayedInvalidateRead(AlignedRead &read, const string &message) {
    logRerouting(read, message);
}

ReadLogger::ReadLogger(AlignedReadStorageFire &fire, size_t threads,
                       const std::experimental::filesystem::path &out_file) :
        AlignedReadStorageListener(fire, "ReadLogger"), logs(threads), os() {
    os.open(out_file);
}
