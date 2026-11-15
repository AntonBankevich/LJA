#include "random_graph_utils.hpp"

namespace ag_test {
    std::string CollapsedGenome(size_t len, uint32_t seed) {
        static const char bases[] = {'A', 'C', 'G', 'T'};
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> dist(0, 2);
        std::string g(1, bases[seed % 4]);
        while (g.size() < len) {
            int c = dist(rng);
            if (bases[c] >= g.back()) c++;
            g += bases[c];
        }
        return g;
    }

    bool ResolveMatching(ag::AssemblyGraph &graph, size_t k) {
        ag::Vertex *core = nullptr;
        for (ag::Vertex &v : graph.vertices())
            if (core == nullptr && v.size() == k && v.inDeg() == 2 && v.outDeg() == 2)
                core = &v;
        if (core == nullptr)
            return false;
        std::vector<ag::Edge *> in, out;
        for (ag::Edge &e : core->incoming())
            in.push_back(&e);
        for (ag::Edge &e : *core)
            out.push_back(&e);
        ag::VertexResolutionPlan plan(*core);
        plan.add(*in[0], *out[0]);
        plan.add(*in[1], *out[1]);
        graph.resolveVertex(*core, plan);
        graph.removeMarked();
        return true;
    }

    bool ResolveByReads(ag::AssemblyGraph &graph, ag::AlignedReadStorage &storage, size_t k) {
        ag::Vertex *core = nullptr;
        for (ag::Vertex &v : graph.vertices())
            if (core == nullptr && v.size() == k && v.inDeg() >= 2 && v.outDeg() >= 2)
                core = &v;
        if (core == nullptr)
            return false;
        ag::VertexResolutionPlan plan(*core);
        for (ag::AlignedRead &read : storage)
            if (read.valid() && !read.getPath().empty()) {
                std::vector<ag::Edge *> edges;
                for (ag::Edge &e : read.getPath().edges())
                    edges.push_back(&e);
                for (size_t i = 0; i + 1 < edges.size(); i++)
                    if (edges[i]->getFinish() == *core)
                        plan.add(*edges[i], *edges[i + 1], 1);
                    else if (edges[i]->getFinish() == core->rc())
                        plan.add(edges[i + 1]->rc(), edges[i]->rc(), 1);
            }
        graph.resolveVertex(*core, plan);
        graph.removeMarked();
        return true;
    }

    size_t FlipAll(ag::AssemblyGraph &graph) {
        std::vector<ag::VertexId> cores;
        for (ag::Vertex &v : graph.vertices())
            if (v.inDeg() == 1 && v.incFront().isSuffix() && v.outDeg() != 1)
                cores.push_back(v.getId());
        size_t flipped = 0;
        for (ag::VertexId vid : cores)
            if (!vid->marked() && vid->inDeg() == 1 && vid->incFront().isSuffix() && vid->outDeg() != 1) {
                graph.flipVertexLockFree(vid->incFront());
                flipped++;
            }
        graph.removeMarked();
        return flipped;
    }
}
