#pragma once

#include "assembly_graph/vertex_resolution.hpp"
#include <vector>
#include <assembly_graph/visualization.hpp>

namespace spg {
    using ag::VertexResolutionPlan;

//    TODO: enforce the multiplexing contract here, for every rule: reads stored in AlignedReadStorages attached to the
//     graph are never broken by multiplexing, i.e. a plan must contain every (incoming, outgoing) pair some stored read
//     passes through. AlignedReadStorageMaintenance::fireResolveVertex relies on it and does not handle missing pairs.
    class DecisionRule {
        bool checkForkForward(ag::Vertex &v) {
            return v.outDeg() == 1 && v.inDeg() > 1 && (v.frontVertex().outDeg() != 1 ||
                (v.frontVertex().front().isSuffix() && v.frontVertex().frontVertex().outDeg() != 1));
        }
    public:
        DecisionRule() {}
        virtual VertexResolutionPlan judgeNontrivial(ag::Vertex & vertex) = 0;
        virtual VertexResolutionPlan judge(ag::Vertex &v);
        virtual bool judgeFlip(ag::GraphPath path) = 0;

        virtual void check() {};// I do not remember what this method is for and there are no implementations. Depricated.

        virtual ~DecisionRule() = default;
    };

}