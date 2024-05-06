#pragma once

#include "supregraph.hpp"
#include <unordered_set>

namespace spg {
//    TODO: make this concurrent by moving uniqueness indicator into vertex itself and locking it every time we need access.
    class UniqueVertexStorage : ResolutionListener {
    private:
        std::unordered_set<ConstVertexId> unique;

        VertexId nextOutAfterDelete(Vertex &cur, Vertex &deleted_core);
        VertexId nextInAfterDelete(Vertex &cur, Vertex &deleted_core);
        void propagateUniquenessForward(Vertex &uv, Vertex &deleted_core);
        void propagateUniqueness(Vertex &uv, Vertex &deleted_core);

    public:
        template<class I>
        UniqueVertexStorage(I begin, I end) {
            for(;begin != end; ++begin) {
                add(*begin);
            }
        }
        UniqueVertexStorage() = default;
        UniqueVertexStorage(UniqueVertexStorage &&) = default;
        UniqueVertexStorage(const UniqueVertexStorage &) = delete;

        void add(Vertex &v);
        void remove(Vertex &v);

        bool isUnique(Vertex &v) const;


        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
    };
}