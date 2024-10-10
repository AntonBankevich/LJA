#pragma once

#include "assembly_graph_base.hpp"
#include "paths.hpp"
#include "compact_path.hpp"
#include "sequences/contigs.hpp"
#include "vertex_resolution.hpp"

namespace ag {
    template<class Traits>
    class ResolutionFire;

//    TODO: this could be split into layers and specified for specific graphs' operations.
//     Looks like too much work that is not necessary and may complicate the code.
    template<class Traits>
    class ResolutionListener {
        friend class ResolutionFire<Traits>;
    private:
        ResolutionFire<Traits> *fire;
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        ResolutionListener(ResolutionFire<Traits> &fire);
        ResolutionListener(ResolutionListener &&other) noexcept;
        ResolutionListener(const ResolutionListener &other) = delete;
        virtual ~ResolutionListener();
        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) {};
        virtual void fireMergePath(const GraphPath <Traits> &path, Vertex &new_vertex) {};
        virtual void fireMergeLoop(const GraphPath <Traits> &path, Vertex &new_vertex) {};
        virtual void fireMergePathToEdge(const GraphPath <Traits> &path, Edge &new_edge) {};
        virtual void fireAddVertex(Vertex &v) {}
        virtual void fireAddEdge(Edge &e) {}
        virtual void fireDeleteVertex(Vertex &v) {}
        virtual void fireDeleteEdge(Edge &e) {}
    };

    template<class Traits>
    class ResolutionFire {
        std::vector<ResolutionListener<Traits> *> listeners;
        typedef typename Traits::Edge Edge;
    protected:
        typedef typename Traits::Vertex Vertex;
    public:
        ResolutionFire() = default;

        ResolutionFire(ResolutionFire &&other) {
            *this = std::move(other);
        }

        ResolutionFire &operator=(ResolutionFire &&other) {
            listeners = std::move(other.listeners);
            for(ResolutionListener<Traits> *listener : listeners) {
                listener->fire = this;
            }
            return *this;
        }

        void addListener(ResolutionListener<Traits> &listener) {
            listeners.emplace_back(&listener);
        }

        void replaceListener(ResolutionListener<Traits> &old_listener, ResolutionListener<Traits> &new_listener) {
            auto it = std::find(listeners.begin(), listeners.end(), &old_listener);
            VERIFY(it != listeners.end());
            *it = &new_listener;
        }

        void removeListener(ResolutionListener<Traits> &listener) {
            listeners.erase(std::find(listeners.begin(), listeners.end(), &listener));
        }

        void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) {
            for (auto *listener: listeners) {
                listener->fireResolveVertex(core, resolution);
            }
        }

        void fireMergePath(const GraphPath <Traits> &path, Vertex &new_vertex) {
            for (auto *listener: listeners) {
                listener->fireMergePath(path, new_vertex);
            }
        }

        void fireMergeLoop(const GraphPath <Traits> &path, Vertex &new_vertex) {
            for (auto *listener: listeners) {
                listener->fireMergeLoop(path, new_vertex);
            }
        }

        void fireMergePathToEdge(const GraphPath <Traits> &path, Edge &new_edge) {
            for (auto *listener: listeners) {
                listener->fireMergePathToEdge(path, new_edge);
            }
        }

        void fireAddVertex(Vertex &vertex) {
            VERIFY_MSG(!vertex.fire_create, "Vertex already fired: " << vertex);
            vertex.fire_create = true;
            for (auto *listener: listeners) {
                listener->fireAddVertex(vertex);
            }
        }

        void fireAddEdge(Edge &edge) {
            VERIFY_MSG(!edge.fire_create, "Edge already fired: " << edge)
            edge.fire_create = true;
            for (auto *listener: listeners) {
                listener->fireAddEdge(edge);
            }
        }

        void fireDeleteVertex(Vertex &vertex) {
            VERIFY_MSG(vertex.fire_create, "Vertex not fired: " << vertex);
            VERIFY_MSG(!vertex.fire_destroy, "Vertex already destroyed: " << vertex);
            vertex.fire_destroy = true;
            for (auto *listener: listeners) {
                listener->fireDeleteVertex(vertex);
            }
        }

        void fireDeleteEdge(Edge &edge) {
            VERIFY_MSG(edge.fire_create, "Edge not fired: " << edge);
            VERIFY_MSG(!edge.fire_destroy, "Edge already destroyed: " << edge);
            edge.fire_destroy = true;
            for (auto *listener: listeners) {
                listener->fireDeleteEdge(edge);
            }
        }
    };

    template<class Traits>
    ResolutionListener<Traits>::ResolutionListener(ResolutionFire<Traits> &fire) : fire(&fire) {
        this->fire->addListener(*this);
    }

    template<class Traits>
    ResolutionListener<Traits>::~ResolutionListener() {
        if (fire != nullptr)
            fire->removeListener(*this);
        fire = nullptr;
    }

    template<class Traits>
    ResolutionListener<Traits>::ResolutionListener(ResolutionListener<Traits> &&other) noexcept {
        fire = other.fire;
        other.fire = nullptr;
        fire->replaceListener(other, *this);
    }
}

namespace ag {
    template<class Traits>
    class LoggingListener : public ResolutionListener<Traits> {
    private:
        std::ostream *outp = nullptr;
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
    public:
        LoggingListener(ResolutionFire<Traits> &fire, std::ostream &out) : ResolutionListener<Traits>(fire), outp(&out) {}

        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) {
            *outp << "FireResolveVertex " << resolution << std::endl;
        };

        virtual void fireMergePath(const GraphPath <Traits> &path, Vertex &new_vertex) {
            *outp << "FireMergePath " << new_vertex << " " << path.lenStr() << std::endl;
        };

        virtual void fireMergeLoop(const GraphPath <Traits> &path, Vertex &new_vertex) {
            *outp << "FireMergeLoop " << new_vertex << " " << path.lenStr() << std::endl;
        };

        virtual void fireMergePathToEdge(const GraphPath <Traits> &path, Edge &new_edge) {
            *outp << "FireMergePathToEdge " << new_edge << " " << path.lenStr() << std::endl;
        };

        virtual void fireAddVertex(Vertex &vertex) {
            *outp << "Fire Add Vertex " << vertex << std::endl;
        }

        virtual void fireAddEdge(Edge &edge) {
            *outp << "Fire Add Edge " << edge << std::endl;
        }

        virtual void fireDeleteVertex(Vertex &vertex) {
            *outp << "Fire delete vertex " << vertex << std::endl;
        }

        virtual void fireDeleteEdge(Edge &edge) {
            *outp << "Fire delete edge " << edge << std::endl;
        }
    };
}
