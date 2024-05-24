#pragma once
#include "paths.hpp"
#include "common/id_index.hpp"
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
        size_t cut_left;
        size_t cut_right;
    public:
        CompactPath() : _start(nullptr), cut_left(0), cut_right(0) {
        }

        CompactPath(Vertex &start, Sequence edges, size_t cut_left = 0, size_t cut_right = 0) :
                _start(&start), _edges(std::move(edges)), cut_left(cut_left), cut_right(cut_right) {
        }

        explicit CompactPath(const GraphPath <Traits> &path) :
                _start(&path.start()), cut_left(path.cutLeft()), cut_right(path.cutRight()) {
            SequenceBuilder sb;
            for (Edge &edge: path.edges()) {
                sb.append(edge.nuclLabel());
            }
            _edges = sb.BuildSequence();
        }

        static CompactPath Subpath(const GraphPath <Traits> &path, size_t left, size_t right) {
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

        GraphPath <Traits> unpack() const {
            if (!valid())
                return {};
            GraphPath<Traits> res(*_start);
            while(res.finish().outDeg() == 1 && res.finish().front().truncSize() == 0) {
                res += res.finish().front();
            }
            for (size_t i = 0; i < _edges.size(); i++) {
                VERIFY(res.finish().hasOutgoing(_edges[i]));
                Edge &edge = res.finish().getOutgoing(_edges[i]);
                res += edge;
                while(res.finish().outDeg() == 1 && res.finish().front().truncSize() == 0) {
                    res += res.finish().front();
                }
            }
            while(res.size() > 0 && res.backEdge().truncSize() == 0 && res.finish().size() <= cut_right)
                res.pop_back();
            res.cutFront(cut_left);
            res.cutBack(cut_right);
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

        size_t cutLeft() const {
            return cut_left;
        }

        size_t cutRight() const {
            return cut_right;
        }

        size_t size() const {
            return _edges.size();
        }

        unsigned char operator[](size_t ind) const {
            return _edges[ind];
        }

        static CompactPath Load(std::istream &os, const IdIndex<Vertex> &index) {
            typename Vertex::id_type id;
            size_t left = 0;
            size_t right = 0;
            std::string path;
            os >> id >> path >> left >> right;
            path = path.substr(2);
            if (id == 0 && path.empty()) {
                return {};
            }
            return {index.getById(id), Sequence(path), left, right};
        }
    };


    template<class Traits>
    inline std::ostream &operator<<(std::ostream &os, const ag::CompactPath<Traits> &cpath) {
        if (cpath.valid())
            return os << cpath.start().getInnerId() << " P:" << cpath.cpath() << " " << cpath.cutLeft() << " "
                      << cpath.cutRight();
        else
            return os << "0 P: 0 0";
    }
}