#pragma once
#include "supregraph.hpp"

namespace spg {
    class UniqueVertexStorage : ResolutionListener {
    private:
        std::unordered_set<ConstVertexId> unique;
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

        void add(VertexId vid) {
            unique.emplace(vid);
            unique.emplace(vid->rc().getId());
        }

        void remove(VertexId vid) {
            unique.erase(vid);
            unique.erase(vid->rc().getId());
        }

        bool isUnique(Vertex &v) const {
            return unique.find(v.getId()) != unique.end();
        }

        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) = 0;
    };
}