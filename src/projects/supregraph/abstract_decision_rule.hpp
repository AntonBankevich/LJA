#pragma once

#include "assembly_graph/vertex_resolution.hpp"
#include <vector>
#include <assembly_graph/visualization.hpp>

namespace spg {
    using ag::VertexResolutionPlan;

    class DecisionRule {
    public:
        virtual VertexResolutionPlan judge(ag::Vertex &v) = 0;
        virtual void check() {};// I do not remember what this method is for and there are no implementations. Depricated.

        virtual ~DecisionRule() = default;
    };

    class RandomDecisionRule : public DecisionRule {
    public:
        VertexResolutionPlan judge(ag::Vertex &v) override;
    };

}