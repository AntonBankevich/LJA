#pragma once


#include "assembly_graph/graph_listeners.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include <unordered_set>

namespace spg {
    using ag::Vertex;
    using ag::Edge;
    using ag::VertexId;

//    TODO: make this concurrent by moving uniqueness indicator into vertex itself and locking it every time we need access.
    class UniqueVertexStorage : public ag::ResolutionListener {
    private:
        std::unordered_set<ag::ConstVertexId> unique;
        size_t unique_threshold;

//        TODO: redo this using corporeal indicator
        VertexId nextOutAfterDelete(Vertex &cur, Vertex &deleted_core);
        VertexId nextInAfterDelete(Vertex &cur, Vertex &deleted_core);
        void propagateUniquenessForward(Vertex &uv, Vertex &deleted_core);
        void propagateUniqueness(Vertex &uv, Vertex &deleted_core);

    public:
        template<class I>
        UniqueVertexStorage(ag::AssemblyGraph &spg, I begin, I end, size_t unique_threshold = 1000000);
        UniqueVertexStorage(ag::AssemblyGraph &spg, const std::function<bool(Vertex &)> &is_unique, size_t unique_threshold = 1000000);
        explicit UniqueVertexStorage(ag::AssemblyGraph &spg, size_t unique_threshold = 1000000);

        UniqueVertexStorage(UniqueVertexStorage &&) = default;
        UniqueVertexStorage(const UniqueVertexStorage &) = delete;

        void add(const Vertex &v);
        void remove(const Vertex &v);

        bool isUnique(const Vertex &v) const;
        std::function<std::string(const Vertex &)> getColorer(const std::string &default_color, const std::string &unique_color) {
            return [this, default_color, unique_color](const Vertex &v) {
                return isUnique(v) ? unique_color : default_color;
            };
        }

        void fireAddVertex(Vertex &vertex) override {if (vertex.innerSize() >= unique_threshold) add(vertex);}
        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;
        void fireMergePath(const ag::RAGraphPath &path, Vertex &vertex) override;
        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &vertex) override {};
        void fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) override {VERIFY(false);};
        void fireMergeLoop(const ag::GraphPath &path, Vertex &vertex) override;
        void fireDeleteVertex(spg::Vertex &v) override {unique.erase(v.getId());}
    };
}

template<class I>
spg::UniqueVertexStorage::UniqueVertexStorage(ag::AssemblyGraph &spg, I begin, I end, size_t unique_threshold) :
            ag::ResolutionListener(spg, "UniqueVertexStorage"), unique_threshold(unique_threshold) {
    for (Vertex &vertex : spg.vertices()) {fireAddVertex(vertex);}
    for(;begin != end; ++begin) {
        add(*begin);
    }
}