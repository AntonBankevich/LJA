#pragma once
#include "component.hpp"
namespace ag {
    template<class Traits>
    class AbstractSplitter {
    public:
        virtual std::vector<Component<Traits>> split(const Component<Traits> &component) const = 0;
        std::vector<Component<Traits>> splitGraph(AssemblyGraph<Traits> &dbg) const {return split(Component<Traits>(dbg));}
    };

    template<class Traits>
    class ConditionSplitter : public AbstractSplitter<Traits> {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex::VertexId VertexId;
        typedef typename Traits::Edge::EdgeId EdgeId;
    private:
        std::function<bool(const typename Traits::Edge &)> splitEdge;
    public:
        explicit ConditionSplitter(std::function<bool(const typename Traits::Edge &)> splitEdge) : splitEdge(std::move(splitEdge)) {}
        std::vector<Component<Traits>> split(const Component<Traits> &comp) const override;
    };

    template<class Traits>
    class CCSplitter : public ConditionSplitter<Traits> {
    private:
        std::function<bool(const typename Traits::Edge &)> splitEdge;
    public:
        explicit CCSplitter() : ConditionSplitter<Traits>([](const typename Traits::Edge &){return false;}) {}
    };

    template<class Traits>
    class LengthSplitter : public ConditionSplitter<Traits> {
    public:
        explicit LengthSplitter(size_t min_len) :
                ConditionSplitter<Traits>([min_len](const typename Traits::Edge& edge){return edge.fullSize() > min_len;}){
        }
    };

    template<class Traits>
    std::vector<Component<Traits>> ConditionSplitter<Traits>::split(const Component<Traits> &comp) const {
        AssemblyGraph<Traits> &dbg = comp.graph();
        std::vector<Component<Traits>> res;
        std::unordered_set<VertexId> visited;
        size_t size = 0;
        for (Vertex &v : comp.verticesUnique()) {
            std::vector<VertexId> queue;
            if (visited.find(v.getId()) != visited.end())
                continue;
            queue.push_back(v.getId());
            queue.push_back(v.rc().getId());
            std::vector<VertexId> component;
            while (!queue.empty()) {
                VertexId vert = queue.back();
                queue.pop_back();
                if (visited.find(vert) != visited.end())
                    continue;
                visited.insert(vert);
                component.emplace_back(vert);
                for (Edge &edge : *vert) {
                    if (!splitEdge(edge) && comp.contains(edge.getFinish())) {
                        queue.emplace_back(edge.getFinish().getId());
                        queue.emplace_back(edge.getFinish().rc().getId());
                    }
                }
            }
            res.emplace_back(dbg, component.begin(), component.end());
            size += res.back().uniqueSize();
        }
        VERIFY(size == comp.uniqueSize());
        return std::move(res);
    }
}