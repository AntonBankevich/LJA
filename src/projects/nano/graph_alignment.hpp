#pragma once

#include <assembly_graph/random_access_paths.hpp>
#include <alignment/alignment_form.hpp>
#include <alignment/ksw_wrapper.hpp>



class GraphAlignment {
private:
    Segment<Contig> segment;
    ag::GraphPath path;
    AlignmentForm alignment;
public:
    GraphAlignment(Segment<Contig> segment, ag::GraphPath path, AlignmentForm alignment) :
                segment(segment), path(std::move(path)), alignment(std::move(alignment)) {
    }

    void operator+=(const GraphAlignment &other) {
        alignment += other.alignment;
        path += other.path;
        segment = segment + other.segment;
    }

    GraphAlignment operator+(const GraphAlignment &other) const {
        GraphAlignment res = *this;
        res += other;
        return std::move(res);
    }
};