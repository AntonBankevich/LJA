#include "aligned_read.hpp"
using namespace ag;

AlignedReadDirection AlignedRead::forward() { return {*this, false}; }


AlignedReadDirection AlignedRead::backward() { return {*this, true}; }


void AlignedRead::setPath(GraphPath  new_path) {
    path = std::move(new_path);
}


void AlignedRead::correct(GraphPath  cpath) {
    VERIFY_MSG(!corrected, "Attempt to correct path while previous correction was not yet applied");
    corrected_path = std::move(cpath);
    corrected = true;
}


void AlignedRead::resetEdgeCodes() {
    path.resetEdgeCodes();
}


std::ostream &operator<<(std::ostream &os, const AlignedReadDirection &dir) {
    os << dir.getRead().getId() << "_" << (dir.isForward() ? "F" : "B") << ":";
    bool first = true;
    for (auto &edge: dir) {
        if (!first)
            os << ",";
        else
            first = false;
        os << edge.getId();
    }
    return os;
}

void AlignedRead::applyCorrection() {
    if (corrected)
        setPath(std::move(corrected_path));
    corrected_path = {};
    corrected = false;
}


void AlignedRead::delayedInvalidate() {
    VERIFY(!corrected);
    corrected_path = {};
    corrected = true;
}
