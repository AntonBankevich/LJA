#pragma once

#include <assembly_graph/paths.hpp>
#include <alignment/alignment_form.hpp>
#include <alignment/ksw_wrapper.hpp>


template<class Traits>
class GraphAlignment {
private:
    Segment<Contig> segment;
    ag::GraphPath<Traits> path;
    AlignmentForm alignment;
public:
    GraphAlignment(Segment<Contig> segment, ag::GraphPath<Traits> path, AlignmentForm alignment) :
                segment(segment), path(std::move(path)), alignment(std::move(alignment)) {
    }

    void operator+=(const GraphAlignment<Traits> &other) {
        alignment += other.alignment;
        path += other.path;
        segment = segment + other.segment;
    }

    GraphAlignment operator+(const GraphAlignment<Traits> &other) const {
        GraphAlignment<Traits> res = *this;
        res += other;
        return std::move(res);
    }
};