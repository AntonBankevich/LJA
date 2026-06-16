#pragma once

#include "assembly_graph_base.hpp"
#include "graph_paths.hpp"
#include "sequences/contigs.hpp"
#include "vertex_resolution.hpp"
#include "common/fire_listeners.hpp"
#include "common/logging.hpp"
#include "alignment/alignment_form.hpp"

namespace ag {

    class ResolutionFire;

//    TODO: this could be split into layers and specified for specific graphs' operations.
//     Looks like too much work that is not necessary and may complicate the code.
    class ResolutionListener : public AbstractListener {
    public:
        ResolutionListener(ResolutionFire &fire, const std::string &name);
        ResolutionListener(ResolutionListener &&other) noexcept = default;
        ResolutionListener &operator=(ResolutionListener &&other) noexcept = default;
        ResolutionListener(const ResolutionListener &other) = delete;
        ResolutionListener &operator=(const ResolutionListener &other) = delete;
        virtual void fireAddVertex(Vertex &v) {}
        virtual void fireAddEdge(Edge &e) {}
        virtual void fireDeleteVertex(Vertex &v) {}
        virtual void fireDeleteEdge(Edge &e) {}
        virtual void fireEdgeToSupreVertex(Vertex &v, Edge &e) {}

        virtual void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {}
        virtual void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {}
        virtual void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {}
        virtual void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &left_al, const AlignmentForm &right_al) {}
        virtual void fireSplitEdge(Edge &edge, const RAGraphPath &split) {}

        virtual void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {}

        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {};
    };

    class ResolutionFire : public AbstractFire {
        std::vector<ResolutionListener *> listeners;
    public:
        ResolutionFire() = default;

        ResolutionFire(ResolutionFire &&other)   noexcept = default;
        ResolutionFire &operator=(ResolutionFire &&other)  noexcept = default;

        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireResolveVertex(core, resolution);
                if(core != core.rc())
                    listener->fireResolveVertex(core.rc(), resolution.RC());
            }
        }

        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireMergePath(path, new_vertex);
                if(new_vertex != new_vertex.rc())
                    listener->fireMergePath(path.RC(), new_vertex.rc());
            }
        }

        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireMergeLoop(path, new_vertex);
                if(new_vertex != new_vertex.rc())
                    listener->fireMergeLoop(path.RC(), new_vertex.rc());
            }
        }

        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireMergePathToEdge(path, new_edge);
                if(new_edge != new_edge.rc())
                    listener->fireMergePathToEdge(path.RC(), new_edge.rc());
            }
        }

        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &left_al, const AlignmentForm &right_al) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireMergeTipsToEdge(new_edge, left, right, left_al, right_al);
                if(new_edge != new_edge.rc())
                    listener->fireMergeTipsToEdge(new_edge.rc(), right.rc(), left.rc(), right_al.RC(), left_al.RC());
            }
        }

        void fireSplitEdge(Edge &edge, const RAGraphPath &split) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireSplitEdge(edge, split);
                if(edge != edge.rc())
                    listener->fireSplitEdge(edge.rc(), split.RC());
            }
        }

        void fireAddVertex(Vertex &vertex) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireAddVertex(vertex);
                if(vertex != vertex.rc())
                    listener->fireAddVertex(vertex.rc());
            }
        }

        void fireAddEdge(Edge &edge) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireAddEdge(edge);
                if(edge != edge.rc())
                    listener->fireAddEdge(edge.rc());
            }
        }

        void fireDeleteVertex(Vertex &vertex) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireDeleteVertex(vertex);
                if(vertex != vertex.rc())
                    listener->fireDeleteVertex(vertex.rc());
            }
        }

        void fireDeleteEdge(Edge &edge) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireDeleteEdge(edge);
                if(edge != edge.rc())
                    listener->fireDeleteEdge(edge.rc());
            }
        }

        void fireEdgeToSupreVertex(Vertex &v, Edge &e) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireEdgeToSupreVertex(v, e);
                if(v != v.rc())
                    listener->fireEdgeToSupreVertex(v.rc(), e.rc());
            }
        }

        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireResetEdgeCodes(logger, threads, graph);
            }
        }
    };

//TODO: Make parallel logger with buffers like read logger

    class LoggingListener : public ResolutionListener {
    private:
        std::ostream *outp = nullptr;
        omp_lock_t writelock = {};
        void lock() { omp_set_lock(&writelock); }
        void unlock() { omp_unset_lock(&writelock); }
    public:
        LoggingListener(ResolutionFire &fire, std::ostream &out) : ResolutionListener(fire, "LoggingListener"), outp(&out) {}

        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override {
            lock();
            *outp << "FireResolveVertex " << resolution << std::endl;
            unlock();
        };

        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override {
            lock();
            *outp << "FireMergePath " << new_vertex << " " << GraphPath(path).str() << std::endl;
            unlock();
        };

        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) override {
            lock();
            *outp << "FireMergeLoop " << new_vertex << " " << path.str() << std::endl;
            unlock();
        };

        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override {
            lock();
            *outp << "FireMergePathToEdge " << new_edge << " " << GraphPath(path).str() << std::endl;
            unlock();
        };

        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override {
            lock();
            *outp << "FireMergeTipsToEdge " << new_edge << " " << left_al.targetLength() << std::endl;
            unlock();
        };

        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override {
            lock();
            *outp << "FireSplitEdge " << edge << " " << split.RC() << std::endl;
            unlock();
        }

        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override {
            lock();
            *outp << "fireEdgeToSupreVertex " << v << " " << e << std::endl;
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
