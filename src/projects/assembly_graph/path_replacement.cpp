#include "path_replacement.hpp"
#include "assembly_graph.hpp"

using namespace ag;

namespace {
//    Offsets of the vertices' starts: each edge moves the start by its rc's trunc size.
    std::vector<size_t> Positions(const RAGraphPath &path) {
        std::vector<size_t> res = {0};
        for (Edge &e : path.edges())
            res.push_back(res.back() + e.rc().truncSize());
        return res;
    }
}

PathReplacement::PathReplacement(const RAGraphPath &old_path, const RAGraphPath &new_path, size_t dominant_index,
                                 const std::unordered_map<EdgeId, EdgeId> &outer) :
        old_path(old_path), new_path(new_path), outer(outer), dominant_index(dominant_index),
        old_positions(Positions(old_path)), new_positions(Positions(new_path)) {
    VERIFY(old_path.size() > 0 && new_path.size() <= 2 && dominant_index <= new_path.size());
    for (size_t i = 0; i <= old_path.size() && dominant_old_index == size_t(-1); i++)
        if (old_path.getVertex(i) == getDominant())
            dominant_old_index = i;
    VERIFY(!startReplaced() || old_path.getStart().outDeg() == 1);
    VERIFY(!endReplaced() || old_path.getFinish().inDeg() == 1);
}

bool PathReplacement::survives(size_t i) const {
    Vertex &v = old_path.getVertex(i);
    for (Vertex &w : new_path.vertices())
        if (v == w)
            return true;
    return false;
}

Sequence PathReplacement::oldCodes() const {
    SequenceBuilder sb;
    for (Edge &e : old_path.edges())
        sb.append(e.getCode());
    return sb.BuildSequence();
}
