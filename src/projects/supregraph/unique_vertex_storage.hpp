#pragma once

#include "supregraph.hpp"
#include "assembly_graph/listeners.hpp"
#include <unordered_set>

namespace spg {
//    TODO: make this concurrent by moving uniqueness indicator into vertex itself and locking it every time we need access.
    class UniqueVertexStorage : ag::ResolutionListener<SPGTraits> {
    private:
        std::unordered_set<ConstVertexId> unique;

        VertexId nextOutAfterDelete(Vertex &cur, Vertex &deleted_core);
        VertexId nextInAfterDelete(Vertex &cur, Vertex &deleted_core);
        void propagateUniquenessForward(Vertex &uv, Vertex &deleted_core);
        void propagateUniqueness(Vertex &uv, Vertex &deleted_core);

    public:
        template<class I>
        UniqueVertexStorage(SupreGraph &spg, I begin, I end);
        UniqueVertexStorage(SupreGraph &spg) : ResolutionListener(spg) {}
        UniqueVertexStorage(UniqueVertexStorage &&) = default;
        UniqueVertexStorage(const UniqueVertexStorage &) = delete;

        void add(const Vertex &v);
        void remove(const Vertex &v);

        bool isUnique(Vertex &v) const;

        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
        void fireMergePath(const GraphPath &path, Vertex &vertex) override;
        void fireMergeLoop(const GraphPath &path, Vertex &vertex) override;
        void fireDeleteVertex(spg::Vertex &v) override {unique.erase(v.getId());}
    };
}

template<class I>
spg::UniqueVertexStorage::UniqueVertexStorage(spg::SupreGraph &spg, I begin, I end) : ResolutionListener(spg) {
    for(;begin != end; ++begin) {
        add(*begin);
    }
}
