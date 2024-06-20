#pragma once

#include "assembly_graph_base.hpp"
#include "graph_paths.hpp"
#include "sequences/contigs.hpp"
#include "vertex_resolution.hpp"
#include "common/fire_listeners.hpp"
#include "common/logging.hpp"
#include "alignment/alignment_form.hpp"

namespace ag {

    template<class Traits>
    class ResolutionFire;

//    TODO: this could be split into layers and specified for specific graphs' operations.
//     Looks like too much work that is not necessary and may complicate the code.
    template<class Traits>
    class ResolutionListener : public AbstractListener {
    public:
        typedef typename Traits::Edge Edge;
        typedef typename Traits::Vertex Vertex;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Vertex::VertexId VertexId;

        explicit ResolutionListener(ResolutionFire<Traits> &fire, const std::string &name);
        ResolutionListener(ResolutionListener &&other) noexcept = default;
        ResolutionListener &operator=(ResolutionListener &&other) noexcept = default;
        ResolutionListener(const ResolutionListener &other) = delete;
        ResolutionListener &operator=(const ResolutionListener &other) = delete;
        virtual void fireAddVertex(Vertex &v) {}
        virtual void fireAddEdge(Edge &e) {}
        virtual void fireDeleteVertex(Vertex &v) {}
        virtual void fireDeleteEdge(Edge &e) {}
        virtual void fireAddSupreVertex(Vertex &v, Edge &e) {}

        virtual void fireMergePath(const std::vector<EdgeId> &path, Vertex &new_vertex) {}
        virtual void fireMergeLoop(const ag::GraphPath <Traits> &path, Vertex &new_vertex) {}
        virtual void fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) {}
        virtual void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &left_al, const AlignmentForm &right_al) {}
        virtual void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) {}

        virtual void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {}

        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) {};
    };

    template<class Traits>
    class ResolutionFire : public AbstractFire {
        std::vector<ResolutionListener<Traits> *> listeners;
        typedef typename Traits::Edge Edge;
        typedef typename Edge::EdgeId EdgeId;
    protected:
        typedef typename Traits::Vertex Vertex;
    public:
        ResolutionFire() = default;

        ResolutionFire(ResolutionFire &&other)   noexcept = default;
        ResolutionFire &operator=(ResolutionFire &&other)  noexcept = default;

        void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireResolveVertex(core, resolution);
            }
        }

        void fireMergePath(const std::vector<EdgeId> &path, Vertex &new_vertex) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireMergePath(path, new_vertex);
            }
        }

        void fireMergeLoop(const GraphPath <Traits> &path, Vertex &new_vertex) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireMergeLoop(path, new_vertex);
            }
        }

        void fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireMergePathToEdge(path, new_edge);
            }
        }

        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &left_al, const AlignmentForm &right_al) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireMergeTipsToEdge(new_edge, left, right, left_al, right_al);
            }
        }

        void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireSplitEdge(edge, split);
            }
        }

        void fireAddVertex(Vertex &vertex) {
            VERIFY_MSG(!vertex.fire_create, "Vertex already fired: " << vertex);
            vertex.fire_create = true;
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireAddVertex(vertex);
            }
        }

        void fireAddEdge(Edge &edge) {
            VERIFY_MSG(!edge.fire_create, "Edge already fired: " << edge)
            edge.fire_create = true;
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireAddEdge(edge);
            }
        }

        void fireDeleteVertex(Vertex &vertex) {
            VERIFY_MSG(vertex.fire_create, "Vertex not fired: " << vertex);
            VERIFY_MSG(!vertex.fire_destroy, "Vertex already destroyed: " << vertex);
            vertex.fire_destroy = true;
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireDeleteVertex(vertex);
            }
        }

        void fireDeleteEdge(Edge &edge) {
            VERIFY_MSG(edge.fire_create, "Edge not fired: " << edge);
            VERIFY_MSG(!edge.fire_destroy, "Edge already destroyed: " << edge);
            edge.fire_destroy = true;
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireDeleteEdge(edge);
            }
        }

        void fireAddSupreVertex(Vertex &v, Edge &e) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireAddSupreVertex(v, e);
            }
        }

        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph<Traits> &graph) {
            for (ResolutionListener<Traits> *listener: getListeners<ResolutionListener<Traits>>()) {
                listener->fireResetEdgeCodes(logger, threads, graph);
            }
        }
    };

    template<class Traits>
    ResolutionListener<Traits>::ResolutionListener(ResolutionFire<Traits> &fire, const std::string &name) : AbstractListener(fire, name) {}
//TODO: Make parallel logger with buffers like read logger
    template<class Traits>
    class LoggingListener : public ResolutionListener<Traits> {
    private:
        std::ostream *outp = nullptr;
        typedef typename Traits::Edge Edge;
        typedef typename Edge::EdgeId EdgeId;
        typedef typename Traits::Vertex Vertex;
        omp_lock_t writelock = {};
        void lock() { omp_set_lock(&writelock); }
        void unlock() { omp_unset_lock(&writelock); }
    public:
        LoggingListener(ResolutionFire<Traits> &fire, std::ostream &out) : ResolutionListener<Traits>(fire, "LoggingListener"), outp(&out) {}

        void fireResolveVertex(Vertex &core, const VertexResolutionResult<Traits> &resolution) override {
            lock();
            *outp << "FireResolveVertex " << resolution << std::endl;
            unlock();
        };

        void fireMergePath(const std::vector<EdgeId> &path, Vertex &new_vertex) override {
            lock();
            *outp << "FireMergePath " << new_vertex << " " << GraphPath<Traits>(path).str() << std::endl;
            unlock();
        };

        void fireMergeLoop(const GraphPath <Traits> &path, Vertex &new_vertex) override {
            lock();
            *outp << "FireMergeLoop " << new_vertex << " " << path.str() << std::endl;
            unlock();
        };

        void fireMergePathToEdge(const std::vector<EdgeId> &path, Edge &new_edge) override {
            lock();
            *outp << "FireMergePathToEdge " << new_edge << " " << GraphPath<Traits>(path).str() << std::endl;
            unlock();
        };

        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override {
            lock();
            *outp << "FireMergeTipsToEdge " << new_edge << " " << left_al.targetLength() << std::endl;
            unlock();
        };

        void fireSplitEdge(Edge &edge, const std::vector<EdgeId> &split) override {
            lock();
            *outp << "FireSplitEdge " << edge << " " << split << std::endl;
            unlock();
        }

        void fireAddSupreVertex(Vertex &v, Edge &e) override {
            lock();
            *outp << "fireAddSupreVertex " << v << " " << e << std::endl;
            unlock();
        }

        void fireAddVertex(Vertex &vertex) override {
            lock();
            *outp << "Fire Add Vertex " << vertex << std::endl;
            unlock();
        }

        void fireAddEdge(Edge &edge) override {
            lock();
            *outp << "Fire Add Edge " << edge << std::endl;
            unlock();
        }

        void fireDeleteVertex(Vertex &vertex) override {
            lock();
            *outp << "Fire delete vertex " << vertex << std::endl;
            unlock();
        }

        void fireDeleteEdge(Edge &edge) override {
            lock();
            *outp << "Fire delete edge " << edge << std::endl;
            unlock();
        }
    };
}
