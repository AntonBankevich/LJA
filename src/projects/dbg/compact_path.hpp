#pragma once
#include "paths.hpp"
#include "id_index.hpp"
#include "sequences/sequence.hpp"
#include <parallel/algorithm>
#include <utility>
namespace ag {

    template<class Traits>
    class CompactPath {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
    private:
        Vertex *_start;
        Sequence _edges = {};
        size_t _first_skip;
        size_t _last_skip;
    public:
        CompactPath() : _start(nullptr), _first_skip(0), _last_skip(0) {
        }

        CompactPath(Vertex &start, Sequence edges, size_t first_skip = 0, size_t last_skip = 0) :
                _start(&start), _edges(std::move(edges)), _first_skip(first_skip), _last_skip(last_skip) {
        }

        explicit CompactPath(const GraphPath<Traits> &path) :
                _start(&path.getVertex(0)), _first_skip(path.leftSkip()), _last_skip(path.rightSkip()) {
            SequenceBuilder sb;
            for (Edge &edge: path.edges()) {
                sb.append(edge.nuclLabel());
            }
            _edges = sb.BuildSequence();
        }

        static CompactPath Subpath(const GraphPath<Traits> &path, size_t left, size_t right) {
            SequenceBuilder sb;
            for (size_t i = left; i < right; i++) {
                sb.append(path[i].contig().nuclLabel());
            }
            return {path.getVertex(left), sb.BuildSequence(), path[left].left, path[right - 1].RC().left};
        }

        bool valid() const {
            VERIFY(_start != nullptr || size() == 0);
            return _start != nullptr;
        }

        GraphPath<Traits> unpack() const {
            if (!valid())
                return {};
            GraphPath<Traits> res(*_start);
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
            return CompactPath(unpack().RC());
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
        static CompactPath Load(std::istream &os, IdIndex<Vertex> &index) {
            typename Vertex::id_type id;
            size_t left = 0;
            size_t right = 0;
            std::string path;
            os >> id >> path >> left >> right;
            path = path.substr(2);
            if(id == 0 && path.empty()) {
                return {};
            }
            return {index.getById(id), Sequence(path), left, right};
        }
    };
}

template<class Traits>
inline std::ostream& operator<<(std::ostream  &os, const ag::CompactPath<Traits> &cpath) {
    if(cpath.valid())
        return os << cpath.start().getInnerId() << " P:" << cpath.cpath() << " " << cpath.leftSkip() << " " << cpath.rightSkip();
    else
        return os << "0 P: 0 0";
}