#pragma once

#include <dbg/paths.hpp>
#include <alignment/alignment_form.hpp>
#include <ksw2/ksw_wrapper.hpp>


template<class Graph>
class GraphAlignment {
private:
    Segment<Contig> segment;
    GraphPath<Graph> path;
    AlignmentForm alignment;
public:
    GraphAlignment(Segment<Contig> segment, GraphPath<Graph> path, AlignmentForm alignment) :
                segment(segment), path(std::move(path)), alignment(std::move(alignment)) {
    }

    void operator+=(const GraphAlignment<Graph> &other) {
        alignment += other.alignment;
        path += other.path;
        segment = segment + other.segment;
    }

    GraphAlignment operator+(const GraphAlignment<Graph> &other) const {
        GraphAlignment<Graph> res = *this;
        res += other;
        return std::move(res);
    }
};