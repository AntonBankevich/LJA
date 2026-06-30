#pragma once
#include "alignment_chain.hpp"
#include "assembly_graph_base.hpp"
#include "graph_listeners.hpp"
#include "graph_paths.hpp"
#include "visualization.hpp"
#include "common/double_linked_list.hpp"
#include "sequences/contigs.hpp"
#include "common/double_linked_list.hpp"

namespace ag {
    // TODO: make reversable
    struct AlignmentFragment {
        Segment<Contig> seg;
        VertexId vertex = {};
        EdgeId edge = {};
        size_t cut_left = 0;
        size_t cut_right = 0;

        AlignmentFragment(Segment<Contig> seg, VertexId vertex = {}, EdgeId edge = {}, size_t cut_left = 0, size_t cut_right = 0) : seg(seg), vertex(vertex), edge(edge), cut_left(cut_left), cut_right(cut_right) {
            VERIFY(!vertex.valid() || !edge.valid());
            if (vertex.valid()) {
                VERIFY(seg.size() == vertex->size() - cut_left - cut_right);
            } else if (edge.valid()) {
                VERIFY(seg.size() == edge->fullSize() - cut_left - cut_right);
            }
        }
    public:
        AlignmentFragment() = default;
        static AlignmentFragment InnerFragment(Segment<Contig> seg, Vertex &vertex, size_t from);
        static AlignmentFragment EdgeSegment(Segment<Contig> seg, Edge &edge, size_t from);

        bool operator==(const AlignmentFragment &other) const {
            return seg == other.seg && vertex == other.vertex && edge == other.edge && cut_left == other.cut_left && cut_right == other.cut_right;
        };
        bool operator!=(const AlignmentFragment &other) const {
            return !(*this == other);
        };

        bool valid() const {return seg.valid();}
        bool isEdgeSegment() const {return !vertex.valid() && edge.valid();}
        bool isInnerFragment() const {return vertex.valid() && !edge.valid();}
        bool isLinked(const AlignmentFragment &other) const;
        bool checkOverlap(const AlignmentFragment &other) const;
        bool checkSeqMatch() const;

        // Only run merge if overlap is true
        AlignmentFragment merge(const AlignmentFragment &other) const;
    };
}

namespace std {
    inline std::ostream &operator<<(std::ostream &os, const ag::AlignmentFragment &fragment) {
        if (fragment.isEdgeSegment()) {
            return os<<fragment.seg << "->" << fragment.edge->getId() << "[" << fragment.cut_left << "," << (fragment.edge->fullSize() - fragment.cut_right) << "]";
        } else if (fragment.isInnerFragment()) {
            return os<<fragment.seg << "->" << fragment.vertex->getId() << "[" << fragment.cut_left << "," << (fragment.vertex->size() - fragment.cut_right) << "]";
        } else {
            VERIFY(!fragment.valid());
            return os << "InvalidFragment";
        }
    }
}