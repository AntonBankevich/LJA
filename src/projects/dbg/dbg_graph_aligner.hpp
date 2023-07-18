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

    class GraphAligner {
    private:
        dbg::SparseDBG &dbg;

        PerfectAlignment<Contig, dbg::Edge> extendLeft(const hashing::KWH &kwh, Contig &contig) const;

        PerfectAlignment<Contig, dbg::Edge> extendRight(const hashing::KWH &kwh, Contig &contig) const;

    public:
        explicit GraphAligner(dbg::SparseDBG &dbg) : dbg(dbg) {
            VERIFY(dbg.alignmentReady());
        }

        DBGGraphPath align(const dbg::EdgePosition &pos, const Sequence &seq) const;

        DBGGraphPath align(const Sequence &seq, dbg::Edge *edge_to, size_t pos_to);

        DBGGraphPath align(const Sequence &seq, const std::string &name = "") const;

        std::vector<PerfectAlignment<Contig, dbg::Edge>> carefulAlign(Contig &contig) const;

        std::vector<PerfectAlignment<dbg::Edge, dbg::Edge>> oldEdgeAlign(
                dbg::Edge &contig) const;

        std::vector<PerfectAlignment<Contig, dbg::Edge>> sparseAlign(Contig &contig) const;
    };
}

namespace std {
    template<class U, class V>
    std::ostream &operator<<(std::ostream &os, const dbg::PerfectAlignment<U, V> &al) {
        return os << al.seg_from << "->" << al.seg_to << std::endl;
    }
}