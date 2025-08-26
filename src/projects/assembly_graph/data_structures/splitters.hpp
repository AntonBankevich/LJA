#pragma once
#include "component.hpp"
namespace ag {
    class AbstractSplitter {
    public:
        virtual std::vector<Component> split(const Component &component) const = 0;
        std::vector<Component> splitGraph(AssemblyGraph &dbg) const {return split(Component(dbg));}
    };

    class ConditionSplitter : public AbstractSplitter {
    private:
        std::function<bool(const Edge &)> splitEdge;
    public:
        explicit ConditionSplitter(std::function<bool(const Edge &)> splitEdge) : splitEdge(std::move(splitEdge)) {}
        std::vector<Component> split(const Component &comp) const override;
    };

    class CCSplitter : public ConditionSplitter {
    public:
        explicit CCSplitter() : ConditionSplitter([](const Edge &){return false;}) {}
    };

    class LengthSplitter : public ConditionSplitter {
    public:
        explicit LengthSplitter(size_t min_len) :
                ConditionSplitter([min_len](const Edge& edge){return edge.fullSize() > min_len;}){
        }
    };
}