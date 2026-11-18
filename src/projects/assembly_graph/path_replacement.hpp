#pragma once
#include "random_access_paths.hpp"
#include <unordered_map>

namespace ag {
//    The geometry of a path replacement (ResolutionListener::fireReplacePath). The old path P = p0 -> ... -> pn is replaced
//    by the new path Q, which spells the same sequence: Q = [q0] -VP-> D -VS-> [qe], each edge optional. D, the dominant
//    vertex, is given by its index in Q and may be a vertex of P. A vertex of P survives if it is a vertex of Q; the others
//    are removed and their information goes to D. Positions are offsets of a vertex's start in the spelled sequence.
//    Contract: a replaced start of P has a single outgoing edge, a replaced finish a single incoming edge. If both are
//    replaced, an edge from P's finish to P's start (a loop of Q) is both an edge into the start and out of the finish
//    and has a single replacement.
    class PathReplacement {
    private:
        const RAGraphPath &old_path;
        const RAGraphPath &new_path;
        const std::unordered_map<EdgeId, EdgeId> &outer;
        size_t dominant_index;
        size_t dominant_old_index = size_t(-1);
        std::vector<size_t> old_positions;
        std::vector<size_t> new_positions;
    public:
        PathReplacement(const RAGraphPath &old_path, const RAGraphPath &new_path, size_t dominant_index,
                        const std::unordered_map<EdgeId, EdgeId> &outer);

        const RAGraphPath &oldPath() const {return old_path;}
        const RAGraphPath &newPath() const {return new_path;}
        Vertex &getDominant() const {return new_path.getVertex(dominant_index);}
        size_t dominantPosition() const {return new_positions[dominant_index];}
//        The index of D in P, or size_t(-1) if D is a new vertex.
        size_t dominantOldIndex() const {return dominant_old_index;}
//        Q's edge into D (virtually prefix) and out of D (virtually suffix).
        bool hasVP() const {return dominant_index > 0;}
        bool hasVS() const {return dominant_index < new_path.size();}
        Edge &vpEdge() const {return new_path.getEdge(dominant_index - 1);}
        Edge &vsEdge() const {return new_path.getEdge(dominant_index);}
        bool startReplaced() const {return old_path.getStart() != new_path.getStart();}
        bool endReplaced() const {return old_path.getFinish() != new_path.getFinish();}
//        Whether P's i-th vertex is a vertex of Q.
        bool survives(size_t i) const;
//        Offset of the start of P's i-th vertex.
        size_t oldPosition(size_t i) const {return old_positions[i];}
//        The replacement of an edge into P's start or out of P's finish (only if that end is replaced).
        Edge &replacement(Edge &old_edge) const {return *outer.at(old_edge.getId());}
//        Whether an edge out of a replaced finish was already replaced as an edge into the replaced start (the loop).
        bool closesLoop(const Edge &out) const {return startReplaced() && out.getFinish() == old_path.getStart();}
//        The concatenated codes of P's edges. They go to Q's edge into D if it has one, else to the replacements of the
//        edges into a replaced start (after their own codes), else to Q's edge out of D, else to the replacements of the
//        edges out of a replaced finish (before their own codes).
        Sequence oldCodes() const;
        bool codesToReplacedStart() const {return !hasVP() && startReplaced();}
    };
}
