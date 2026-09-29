#pragma once

#include "error_correction.hpp"

namespace ag {
//    Splits a read's path into tips and internal segments delimited by vertices marked
//    VertexReliability::reliable (see Vertex::isReliable()), then reassembles the read from
//    per-part corrections returned by the pure virtual methods below. Some other algorithm is
//    expected to have already marked reliable vertices on the graph before this one runs -- this
//    class only consumes the marker, it never inspects DBG- or Supregraph-specific state, so it
//    works unchanged on either graph type.
    class AbstractReliableVertexSplittingCorrectionAlgorithm : public ag::AbstractCorrectionAlgorithm {
    public:
        using ag::AbstractCorrectionAlgorithm::AbstractCorrectionAlgorithm;

        std::string correctRead(const std::string &name, ag::GraphPath &read_path) final;

    protected:
//        segment.getStart() and segment.getFinish() are both reliable vertices (possibly the same
//        vertex, if two reliable vertices are adjacent). The returned path must start and finish at
//        exactly the same two vertices as segment.
        virtual std::pair<ag::GraphPath, std::string> correctSegment(const ag::GraphPath &segment) = 0;

//        tip.getStart() is a reliable vertex; tip extends away from it towards a read end that has
//        no further reliable vertex before it. Incoming tips (at the very start of the read) are
//        passed in RC'd form, so tip.getStart() is always the reliable anchor and the tip always
//        extends "forward" from it, regardless of which end of the read it originally came from.
//        The returned path must start at tip.getStart().
        virtual std::pair<ag::GraphPath, std::string> correctTip(const ag::GraphPath &tip) = 0;
    };
}
