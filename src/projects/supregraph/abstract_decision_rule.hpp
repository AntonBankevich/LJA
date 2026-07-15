#pragma once

#include "assembly_graph/vertex_resolution.hpp"
#include <vector>
#include <assembly_graph/visualization.hpp>

namespace spg {
    using ag::VertexResolutionPlan;

    class DecisionRule {
        bool checkForkForward(ag::Vertex &v) {
            return v.outDeg() == 1 && v.inDeg() > 1 && (v.front().getFinish().outDeg() != 1 ||
                (v.front().getFinish().front().isSuffix() && v.front().getFinish().front().getFinish().outDeg() != 1));
        }
    public:
        DecisionRule() {}
        virtual VertexResolutionPlan judgeNontrivial(ag::Vertex &v) = 0;
        virtual VertexResolutionPlan judge(ag::Vertex &v);

        virtual void check() {};// I do not remember what this method is for and there are no implementations. Depricated.

        virtual ~DecisionRule() = default;
    };

}