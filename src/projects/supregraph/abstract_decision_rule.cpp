#include "abstract_decision_rule.hpp"

//TODO: move filling from read set here.
//TODO: add edge detachment to resolution plan and move final check here.
ag::VertexResolutionPlan spg::DecisionRule::judge(ag::Vertex &v) {
    // if (v.outDeg() == 1 && !checkForkForward(v))
    //     return {v};
    // if (v.rc().outDeg() == 1 && !checkForkForward(v.rc()))
    //     return {v};
//    These two VERIFYs used to assume judge() is only ever called with both sides branching -- not true
//    at small scale (e.g. n-in-1-out/1-in-n-out core vertices), so removed for now; the dispatch below
//    already handles every degree combination without them.
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
    return judgeNontrivial(v);
}