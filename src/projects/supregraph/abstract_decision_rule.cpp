#include "abstract_decision_rule.hpp"

ag::VertexResolutionPlan spg::DecisionRule::judge(ag::Vertex &v) {
    if (v.outDeg() == 1 && !checkForkForward(v))
        return {v};
    if (v.rc().outDeg() == 1 && !checkForkForward(v.rc()))
        return {v};
    return judgeNontrivial(v);
}

ag::VertexResolutionPlan spg::RandomDecisionRule::judgeNontrivial(ag::Vertex &v) {
    VertexResolutionPlan res(v);
    auto out_it = v.begin();
    auto inc = v.incoming();
    auto in_it = inc.begin();
    while (out_it != v.end() || in_it != inc.end()) {
        if (out_it == v.end()) --out_it;
        if (in_it == inc.end()) --in_it;
        res.add(*in_it, *out_it);
        ++in_it;
        ++out_it;
    }
    return std::move(res);
}
