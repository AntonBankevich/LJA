#include "splitters.hpp"
using namespace ag;

std::vector<Component> ConditionSplitter::split(const Component &comp) const {
    AssemblyGraph &dbg = comp.getGraph();
    std::vector<Component> res;
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
