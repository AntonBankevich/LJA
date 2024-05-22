#pragma once

#include "vertex_resolution.hpp"
#include "supregraph_base.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/paths.hpp"
#include "assembly_graph/compact_path.hpp"
#include "sequences/contigs.hpp"

namespace spg {
    class ResolutionFire;

    class ResolutionListener {
    private:
        ResolutionFire *fire;
    public:
        ResolutionListener(ResolutionFire &fire);

        ResolutionListener(ResolutionListener &&other) noexcept;

        ResolutionListener(const ResolutionListener &other) = delete;

        virtual ~ResolutionListener();

        virtual void fireResolveVertex(spg::Vertex &core, const spg::VertexResolutionResult &resolution) = 0;

        virtual void fireMergePath(const GraphPath &path, Vertex &new_vertex) {};

        virtual void fireAddVertex(spg::Vertex &v) {}

        virtual void fireAddEdge(spg::Edge &e) {}

        virtual void fireDeleteVertex(spg::Vertex &v) {}

        virtual void fireDeleteEdge(spg::Edge &e) {}
    };

    class ResolutionFire {
        std::vector<ResolutionListener *> listeners;
    public:
        void addListener(ResolutionListener &listener) {
            listeners.emplace_back(&listener);
        }

        void replaceListener(ResolutionListener &old_listener, ResolutionListener &new_listener) {
            auto it = std::find(listeners.begin(), listeners.end(), &old_listener);
            VERIFY(it != listeners.end());
            *it = &new_listener;
        }

        void removeListener(ResolutionListener &listener) {
            listeners.erase(std::find(listeners.begin(), listeners.end(), &listener));
        }

        void fireResolveVertex(spg::Vertex &core, const spg::VertexResolutionResult &resolution) {
            for (auto *listener: listeners) {
                listener->fireResolveVertex(core, resolution);
            }
        }

        void fireMergePath(const GraphPath &path, Vertex &new_vertex) {
            for (auto *listener: listeners) {
                listener->fireMergePath(path, new_vertex);
            }
        }

        void fireAddVertex(spg::Vertex &vertex) {
            VERIFY(!vertex.fire_create)
            vertex.fire_create = true;
            for (auto *listener: listeners) {
                listener->fireAddVertex(vertex);
            }
        }

        void fireAddEdge(spg::Edge &edge) {
            VERIFY(!edge.fire_create)
            edge.fire_create = true;
            for (auto *listener: listeners) {
                listener->fireAddEdge(edge);
            }
        }

        void fireDeleteVertex(spg::Vertex &vertex) {
            VERIFY(vertex.fire_create);
            VERIFY(!vertex.fire_destroy);
            vertex.fire_destroy = true;
            for (auto *listener: listeners) {
                listener->fireDeleteVertex(vertex);
            }
        }

        void fireDeleteEdge(spg::Edge &edge) {
            VERIFY(edge.fire_create);
            VERIFY(!edge.fire_destroy);
            edge.fire_destroy = true;
            for (auto *listener: listeners) {
                listener->fireDeleteEdge(edge);
            }
        }
    };
}