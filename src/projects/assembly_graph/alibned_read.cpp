#include "aligned_read.hpp"
using namespace ag;

AlignedReadDirection AlignedRead::forward() { return {*this, false}; }


AlignedReadDirection AlignedRead::backward() { return {*this, true}; }


//    A path set anew has fresh decks, numbered anew, so its suffix ends start at its current ends. Graph listeners change
//    paths in place by pushes and pops instead, which keep the deck positions; the exceptions are applying a correction
//    (a SuffixTracker adds the corrected path with its ends) and DBG edge splitting (single-edge paths, which hold no
//    suffixes).
void AlignedRead::setPath(GraphPath  new_path) {
    path = std::move(new_path);
    resetSuffixEnds(path);
}


void AlignedRead::correct(GraphPath  cpath) {
    VERIFY_MSG(!corrected, "Attempt to correct path while previous correction was not yet applied");
    corrected_path = std::move(cpath);
    corrected = true;
}


void AlignedRead::resetEdgeCodes() {
    path.resetEdgeCodes();
    resetSuffixEnds(path);
}

void AlignedRead::resetSuffixEnds(const GraphPath &new_path) {
    suffix_fend = new_path.getFSplits().end();
    suffix_rend = new_path.getRSplits().end();
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
