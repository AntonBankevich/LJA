#pragma once

#include "paths.hpp"
#include "sparse_dbg.hpp"

namespace dbg {
    template<class U, class V>
    class PerfectAlignment {
    public:
        Segment<U> seg_from;
        Segment<V> seg_to;
        PerfectAlignment(const Segment<U> &seg_from_, const Segment<V> &seg_to_) : seg_from(seg_from_), seg_to(seg_to_) {
            VERIFY(seg_from_.size() == seg_to_.size());
        }
        size_t size() {return seg_from.size();}
        PerfectAlignment RC() const {
            return {seg_from.RC(), seg_to.RC()};
        }
        bool operator<(const PerfectAlignment<U, V> &other) const {
            if(seg_to != other.seg_to)
                return seg_to < other.seg_to;
            else
                return seg_from < other.seg_from;
        }
    };

    class KmerIndex {
    public:
        typedef std::unordered_map<hashing::htype , Vertex*, hashing::alt_hasher<hashing::htype>> vertex_map_type;
        typedef std::unordered_map<hashing::htype, Vertex*, hashing::alt_hasher<hashing::htype>>::iterator vertex_iterator_type;
        typedef std::unordered_map<hashing::htype, EdgePosition, hashing::alt_hasher<hashing::htype>> anchor_map_type;

    private:
        hashing::RollingHash hasher_;
        mutable vertex_map_type v;
        anchor_map_type anchors;
        size_t w = 0;
        bool anchors_filled = false;

        PerfectAlignment<Contig, dbg::Edge> extendLeft(const hashing::MovingKWH &kwh, Contig &contig) const;
        PerfectAlignment<Contig, dbg::Edge> extendRight(const hashing::MovingKWH &kwh, Contig &contig) const;

    public:
        explicit KmerIndex(hashing::RollingHash _hasher) : hasher_(_hasher) {
        }

        explicit KmerIndex(SparseDBG &dbg);

        const hashing::RollingHash &hasher() const {return hasher_;}

        void addVertex(Vertex &vert) {
            if(vert.isCanonical())
                v[vert.getHash()] = &vert;
            else
                v[vert.getHash()] = &vert.rc();
        }

        bool containsVertex(const hashing::htype &hash) const {return v.find(hash) != v.end();}
        Vertex &getVertex(hashing::htype hash, bool canonical = true) const {
            return canonical ? *v.find(hash)->second : v.find(hash)->second->rc();
        }
        Vertex &getVertex(const hashing::KWH &kwh) const;
        Vertex &getVertex(const Sequence &seq) const;
        Vertex &getVertex(const Vertex &other_graph_vertex) const;
        bool isAnchor(hashing::htype hash) const {return anchors.find(hash) != anchors.end();}
        EdgePosition getAnchor(const hashing::KWH &kwh) const;
        bool alignmentReady() const {return anchors_filled;}
        size_t minReadLen() const {
            if(anchors_filled)
                return w + hasher().getK() - 1;
            else
                return size_t(-1);
        }

        std::vector<hashing::MovingKWH> extractVertexPositions(const Sequence &seq, size_t max = size_t(-1)) const {
            std::vector<hashing::MovingKWH> res;
            hashing::MovingKWH kwh(hasher(), seq, 0);
            while (true) {
                if (containsVertex(kwh.hash())) {
                    res.emplace_back(kwh);
                }
                if (!kwh.hasNext() || res.size() == max)
                    break;
                kwh = kwh.next();
            }
            return std::move(res);
        }

        void noAnchors() {
            anchors_filled = true;
        }

        void fillAnchors(logging::Logger &logger, size_t threads, SparseDBG &dbg, size_t w);
        void fillAnchors(logging::Logger &logger, size_t threads, SparseDBG &dbg, size_t w,
                         const std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> &to_add);


    public:

        dbg::GraphPath align(const dbg::EdgePosition &pos, const Sequence &seq) const;

        dbg::GraphPath align(const Sequence &seq, dbg::Edge *edge_to, size_t pos_to);

        dbg::GraphPath align(const Sequence &seq, const std::string &name = "") const;

        std::vector<PerfectAlignment<Contig, dbg::Edge>> carefulAlign(Contig &contig) const;

        std::vector<PerfectAlignment<dbg::Edge, dbg::Edge>> oldEdgeAlign(
                dbg::Edge &contig) const;

        std::vector<PerfectAlignment<Contig, dbg::Edge>> sparseAlign(Contig &contig) const;
    };
}

namespace std {
    template<class U, class V>
    std::ostream &operator<<(std::ostream &os, const dbg::PerfectAlignment<U, V> &al) {
        return os << al.seg_from << "->" << al.seg_to;
    }
}