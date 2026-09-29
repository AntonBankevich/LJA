#include "abstract_decision_rule.hpp"

//TODO: move filling from read set here.
//TODO: add edge detachment to resolution plan and move final check here.
ag::VertexResolutionPlan spg::DecisionRule::judge(ag::Vertex &v) {
    // if (v.outDeg() == 1 && !checkForkForward(v))
    //     return {v};
    // if (v.rc().outDeg() == 1 && !checkForkForward(v.rc()))
    //     return {v};
    VERIFY(v.inDeg() != 1);
    VERIFY(v.outDeg() != 1);
    // if (v.inDeg() == 1 && v.outDeg() > 1) {
    //     if (judgeFlip(v.incFrontVertex()))
    //         return VertexResolutionPlan::SimplePlan(v);
    //     else
    //         return {v};
    // }
    // if (v.inDeg() > 1 && v.outDeg() == 1) {
    //     if (judgeFlip(v.frontVertex()))
    //         return VertexResolutionPlan::SimplePlan(v);
    //     else
    //         return {v};
    // }
    if (v.inDeg() <= 1 && v.outDeg() <= 1)
        return VertexResolutionPlan::SimplePlan(v);
    if (v.inDeg() == 0 || v.outDeg() == 0 && v.inDeg() + v.outDeg() > 1)
        return {v};
    VERIFY(v.inDeg() > 1);
    VERIFY(v.outDeg() > 1);
    return judgeNontrivial(v);
}