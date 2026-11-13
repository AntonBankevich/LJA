#pragma once


#include "spg_listener.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "libcuckoo/cuckoohash_map.hh"
#include <unordered_set>

namespace spg {
    using ag::Vertex;
    using ag::Edge;
    using ag::VertexId;

//    Merges may run in parallel, so the set is concurrent (a map whose values are unused).
    class UniqueVertexStorage : public SPGListener {
    private:
        libcuckoo::cuckoohash_map<ag::ConstVertexId, bool> unique;
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
//        Once a vertex is unique, every future superstring of it must be unique too.
        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireMergePath(const ag::RAGraphPath &path, Vertex &vertex) override;
        void fireMergeLoop(const ag::GraphPath &path, Vertex &vertex) override;
        void fireDeleteVertex(spg::Vertex &v) override {unique.erase(v.getId());}
    };
}

template<class I>
spg::UniqueVertexStorage::UniqueVertexStorage(ag::AssemblyGraph &spg, I begin, I end, size_t unique_threshold) :
            SPGListener(spg, "UniqueVertexStorage"), unique_threshold(unique_threshold) {
    for (Vertex &vertex : spg.vertices()) {fireAddVertex(vertex);}
    for(;begin != end; ++begin) {
        add(*begin);
    }
}