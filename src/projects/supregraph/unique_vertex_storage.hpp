#pragma once

#include "supregraph.hpp"
#include "assembly_graph/graph_listeners.hpp"
#include <unordered_set>

namespace spg {
//    TODO: make this concurrent by moving uniqueness indicator into vertex itself and locking it every time we need access.
    class UniqueVertexStorage : public ag::ResolutionListener<SPGTraits> {
    private:
        std::unordered_set<ConstVertexId> unique;

//        TODO: redo this using corporeal indicator
        VertexId nextOutAfterDelete(Vertex &cur, Vertex &deleted_core);
        VertexId nextInAfterDelete(Vertex &cur, Vertex &deleted_core);
        void propagateUniquenessForward(Vertex &uv, Vertex &deleted_core);
        void propagateUniqueness(Vertex &uv, Vertex &deleted_core);

    public:
        template<class I>
        UniqueVertexStorage(SupreGraph &spg, I begin, I end);
        UniqueVertexStorage(SupreGraph &spg, const std::function<bool(Vertex &)> &is_unique);
        explicit UniqueVertexStorage(SupreGraph &spg) : ag::ResolutionListener<SPGTraits>(spg, "UniqueVertexStorage") {}
        UniqueVertexStorage(UniqueVertexStorage &&) = default;
        UniqueVertexStorage(const UniqueVertexStorage &) = delete;

        void add(const Vertex &v);
        void remove(const Vertex &v);

        bool isUnique(Vertex &v) const;
        std::function<std::string(Vertex &)> getColorer(const std::string &default_color, const std::string &unique_color) {
            return [this, default_color, unique_color](Vertex &v) {
                return isUnique(v) ? unique_color : default_color;
            };
        }

        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
        void fireMergePath(const std::vector<EdgeId> &path, Vertex &vertex) override;
        void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) override {VERIFY(false);};
        void fireMergeLoop(const ag::GraphPath<SPGTraits> &path, Vertex &vertex) override;
        void fireDeleteVertex(spg::Vertex &v) override {unique.erase(v.getId());}
    };
}

template<class I>
spg::UniqueVertexStorage::UniqueVertexStorage(spg::SupreGraph &spg, I begin, I end) : ag::ResolutionListener<SPGTraits>(spg, "UniqueVertexStorage") {
    for(;begin != end; ++begin) {
        add(*begin);
    }
}