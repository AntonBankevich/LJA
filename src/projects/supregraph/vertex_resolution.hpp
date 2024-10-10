#pragma once

#include "supregraph_base.hpp"
#include "assembly_graph/vertex_resolution.hpp"
#include <vector>

namespace spg {
    typedef ag::InOutEdgePair<SPGTraits> InOutEdgePair;
    typedef ag::VertexResolutionResult<SPGTraits> VertexResolutionResult;

    class VertexResolutionPlan {
    private:
        VertexId v;
        mutable bool sorted = true;
        mutable std::vector<InOutEdgePair> edge_pairs;

        void sort() const;
    public:
        VertexResolutionPlan(Vertex &v) : v(v.getId()) {} // NOLINT(google-explicit-constructor)
        VertexResolutionPlan RC() const {
            VertexResolutionPlan res(v->rc());
            for(const InOutEdgePair &ep : edge_pairs) {
                res.add(ep.RC());
            }
            return std::move(res);
        }

        Vertex &getCore() const {return *v;}
        void add(const InOutEdgePair &edgePair);
        void add(Edge &edge1, Edge &edge2) {add({edge1, edge2});}

        bool empty() const {return edge_pairs.empty();}
        bool incConnected(Edge &edge) const;
        bool outConnected(Edge &edge) const;
        bool allConnected() const;
        IterableStorage<std::vector<InOutEdgePair>::const_iterator> connections() const;
        IterableStorage<SkippingIterator<std::vector<InOutEdgePair>::const_iterator>> connectionsUnique() const;
    };

    inline std::ostream &operator<<(std::ostream &stream, const VertexResolutionPlan &vr) {
        stream << "VRResult." << vr.getCore().getId() << ":";
        for(const InOutEdgePair & it : vr.connections()) {
            stream << "(" << it.incoming().getId() << "|" << it.outgoing().getId() << ")";
        }
        return stream;
    }

    class DecisionRule {
    public:
        virtual VertexResolutionPlan judge(Vertex &v) = 0;
        virtual void check() {};

        virtual ~DecisionRule() = default;
    };

    class RandomDecisionRule : public DecisionRule {
    public:
        VertexResolutionPlan judge(Vertex &v) override {
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
    };
}