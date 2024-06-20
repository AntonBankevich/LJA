#pragma once

#include "assembly_graph/assembly_graph.hpp"
#include <assembly_graph/random_access_paths.hpp>
#include "sequences/contigs.hpp"

namespace spg {
    class SPGVertex;
    class SPGEdge;

    typedef Position<SPGEdge> EdgePosition;
    typedef Segment<SPGEdge> EdgeSegment;

    class SPGVertexData {
    private:
        bool cyclic;
        bool inf_left;
        bool inf_right;
    public:
        SPGVertexData(bool cyclic, bool inf_left, bool inf_right) : cyclic(cyclic), inf_left(inf_left),
                                                                    inf_right(inf_right) {
        }

        bool isCyclic() const { return cyclic; }

        bool isInfLeft() const { return inf_left; }

        bool isInfRight() const { return inf_right; }

        SPGVertexData RC() const {
            return {cyclic, inf_right, inf_left};
        }
    };

    class SPGEdgeData {
    public:
        SPGEdgeData RC() const { return {}; }
    };


    struct SPGTraits {
        typedef SPGVertex Vertex;
        typedef SPGEdge Edge;
        typedef SPGEdgeData EdgeData;
        typedef SPGVertexData VertexData;
    };

//    TODO: switch to new path without random access
    typedef typename ag::GraphPath<SPGTraits> GraphPath;
//    typedef typename old::ag::PathPosition<SPGTraits> PathPosition;
    typedef typename ag::PathDirection<SPGTraits> PathDirection;
    typedef ag::VertexResolutionResult<SPGTraits> VertexResolutionResult;


    class SPGVertex : public ag::BaseVertex<SPGTraits>, public SPGVertexData {
    public:
        explicit SPGVertex(id_type id, Sequence seq, SPGVertexData data) : ag::BaseVertex<SPGTraits>(id,
                                                                                                     std::move(seq)),
                                                                           SPGVertexData(data) {}

        explicit SPGVertex(id_type id, bool canonical, SPGVertexData data) : ag::BaseVertex<SPGTraits>(id, canonical),
                                                                             SPGVertexData(data) {}

//        SPGVertex(): seq(""), id(0), label("") {VERIFY(false);}
        SPGVertex(const SPGVertex &) = delete;
        ~SPGVertex() override {
            VERIFY_MSG(fire_destroy, getId());
        }

        bool isCore();


//        SPGVertex(SPGVertex && v) noexcept : seq(std::move(v.seq)), id(v.id), label(std::move(v.label)) {VERIFY (v.outgoing.empty() && v._rc == nullptr);}
    };


    class SPGEdge : public ag::BaseEdge<SPGTraits>, public SPGEdgeData {
    public:
        explicit SPGEdge(id_type id, SPGVertex &start, SPGVertex &end, Sequence _seq, SPGEdgeData = {}) :
                ag::BaseEdge<SPGTraits>(id, start, end, std::move(_seq)), SPGEdgeData() {}
        size_t intCov() const {return 0;}
        double getCoverage() const {return 0;}
        SPGEdge(const SPGEdge &) = delete;
        ~SPGEdge() override {
            VERIFY(fire_destroy);
        }
    };

    typedef SPGEdge Edge;
    typedef SPGVertex Vertex;
    typedef SPGEdge::EdgeId EdgeId;
    typedef SPGVertex::VertexId VertexId;
    typedef SPGEdge::ConstEdgeId ConstEdgeId;
    typedef SPGVertex::ConstVertexId ConstVertexId;
}