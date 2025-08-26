#include "abstract_decision_rule.hpp"

ag::VertexResolutionPlan spg::RandomDecisionRule::judge(ag::Vertex &v) {
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
