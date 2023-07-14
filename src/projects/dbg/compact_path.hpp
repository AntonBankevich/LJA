#pragma once
#include "sparse_dbg.hpp"
#include "graph_algorithms.hpp"
#include "graph_stats.hpp"
#include <parallel/algorithm>
#include <utility>
namespace dbg {

    class CompactPath {
    private:
        Vertex *_start;
        Sequence _edges;
        size_t _first_skip;
        size_t _last_skip;
    public:
        CompactPath() : _start(nullptr), _first_skip(0), _last_skip(0) {
        }

        CompactPath(Vertex &start, Sequence edges, size_t first_skip = 0, size_t last_skip = 0) :
                _start(&start), _edges(std::move(edges)), _first_skip(first_skip), _last_skip(last_skip) {
        }

        explicit CompactPath(const GraphPath &path) :
                _start(&path.getVertex(0)), _first_skip(path.leftSkip()), _last_skip(path.rightSkip()) {
            std::vector<char> edges;
            for (Edge &edge: path.edges()) {
                edges.push_back(edge.truncSeq()[0]);
            }
            _edges = Sequence(edges);
        }

        static CompactPath Subpath(const GraphPath &path, size_t left, size_t right) {
            std::vector<char> edges;
            for (size_t i = left; i < right; i++) {
                edges.push_back(path[i].contig().truncSeq()[0]);
            }
            return {path.getVertex(left), Sequence(edges), path[left].left, path[right - 1].RC().left};
        }

        bool valid() const {
            VERIFY(_start != nullptr || size() == 0);
            return _start != nullptr;
        }

        GraphPath getAlignment() const {
            if (!valid())
                return {};
            GraphPath res;
            Vertex *cur = _start;
            for (size_t i = 0; i < _edges.size(); i++) {
                VERIFY(cur->hasOutgoing(_edges[i]));
                Edge &edge = cur->getOutgoing(_edges[i]);
                res += edge;
                cur = &edge.getFinish();
            }
            res.cutFront(_first_skip);
            res.cutBack(_last_skip);
            return std::move(res);
        }

        CompactPath RC() const {
            if (!valid())
                return {};
            return CompactPath(getAlignment().RC());
        }

//    Vertex &start() {
//        return *_start;
//    }

        Vertex &start() const {
            return *_start;
        }

        const Sequence &cpath() const {
            return _edges;
        }

        size_t leftSkip() const {
            return _first_skip;
        }

        size_t rightSkip() const {
            return _last_skip;
        }

        size_t size() const {
            return _edges.size();
        }

        unsigned char operator[](size_t ind) const {
            return _edges[ind];
        }

        static CompactPath Load(std::istream &os, SparseDBG &dbg) {
            hashing::htype hash;
            bool canonical;
            size_t left = 0;
            size_t right = 0;
            std::string path;
            os >> hash >> canonical >> path >> left >> right;
            path = path.substr(2);
            if(hash == 0 && path.size() == 0) {
                return {};
            }
            return {dbg.getVertex(hash, canonical), Sequence(path), left, right};
        }
    };
}
inline std::ostream& operator<<(std::ostream  &os, const dbg::CompactPath &cpath) {
    if(cpath.valid())
        return os << cpath.start().hash() << " " << cpath.start().isCanonical() << " P:" << cpath.cpath() << " " << cpath.leftSkip() << " " << cpath.rightSkip();
    else
        return os << "0 0 P: 0 0";
}