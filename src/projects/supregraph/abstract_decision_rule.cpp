#include "abstract_decision_rule.hpp"

//TODO: move filling from read set here.
//TODO: add edge detachment to resolution plan and move final check here.
ag::VertexResolutionPlan spg::DecisionRule::judge(ag::Vertex &v) {
    // if (v.outDeg() == 1 && !checkForkForward(v))
    //     return {v};
    // if (v.rc().outDeg() == 1 && !checkForkForward(v.rc()))
    //     return {v};
    if (v.outDeg() == 1 && v.frontVertex().outDeg() == 1)
        return {v};
    if (v.rc().outDeg() == 1 && v.rc().frontVertex().outDeg() == 1)
        return {v};
    return judgeNontrivial(v);
}