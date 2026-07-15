#include "abstract_decision_rule.hpp"

ag::VertexResolutionPlan spg::DecisionRule::judge(ag::Vertex &v) {
    if (v.outDeg() == 1 && !checkForkForward(v))
        return {v};
    if (v.rc().outDeg() == 1 && !checkForkForward(v.rc()))
        return {v};
    return judgeNontrivial(v);
}