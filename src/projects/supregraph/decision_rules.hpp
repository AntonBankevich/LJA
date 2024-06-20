#pragma once
#include "vertex_resolution.hpp"
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
        PathIndex * storage;
        UniqueVertexStorage const * unique_storage;
        EdgeId getUniqueDisconnectedInc(const VertexResolutionPlan &plan);
    public:
        AndreyRule(const AndreyRule &other) = delete;
        explicit AndreyRule(PathIndex &read_storage, const spg::UniqueVertexStorage &unique_storage) :
                storage(&read_storage), unique_storage(&unique_storage) {}

        void loopHeuristic(VertexResolutionPlan &res) const;
        void uniqueHeuristic(VertexResolutionPlan &res);
        void noChoiceHeuristic(VertexResolutionPlan &res);

        VertexResolutionPlan judge(Vertex &v) override;
        void check() override {
            storage->checkReadIndexConsistency();
        }
    };
}