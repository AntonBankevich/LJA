#pragma once

#include <assembly_graph/data_structures/suffix_tracker.hpp>
#include "abstract_decision_rule.hpp"
#include "read_storage.hpp"

namespace spg {
//    class ChainRule : public DecisionRule {
//    private:
//        PathIndex *storage;
//        size_t k;
//
//        size_t getDiveSize(Edge &edge);
//    public:
//        ChainRule(const ChainRule &other) = delete;
//        ChainRule(PathIndex &storage, size_t k) : storage(&storage), k(k) {}
//
//        VertexResolutionPlan judge(Vertex &v) override;
//    };

    class AndreyRule: public DecisionRule {
    private:
        ag::SuffixTracker * suffixes;
        UniqueVertexStorage const * unique_storage;
        ag::EdgeId getUniqueDisconnectedInc(const VertexResolutionPlan &plan);
    public:
        AndreyRule(const AndreyRule &other) = delete;
        explicit AndreyRule(ag::SuffixTracker &suffixes, const spg::UniqueVertexStorage &unique_storage) :
                suffixes(&suffixes), unique_storage(&unique_storage) {}

        void loopHeuristic(VertexResolutionPlan &res) const;
        void uniqueHeuristic(VertexResolutionPlan &res);
        void noChoiceHeuristic(VertexResolutionPlan &res);

        VertexResolutionPlan judge(Vertex &v) override;
        void check() override {
        }
    };

    class ObviousRule: public DecisionRule {
    private:
        ag::SuffixTracker * suffixes;
        size_t complex_support = 4;
        size_t simple_support = 2;
        size_t max_complex_length = 4000;
        size_t max_simple_length = 7000;
    public:
        ObviousRule(const AndreyRule &other) = delete;
        explicit ObviousRule(ag::SuffixTracker &suffixes) :
                suffixes(&suffixes) {}

        VertexResolutionPlan judge(Vertex &v) override;

        void check() override {
        }
    };
}