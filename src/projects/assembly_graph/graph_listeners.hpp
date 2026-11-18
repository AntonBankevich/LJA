#pragma once

#include "assembly_graph_base.hpp"
#include "graph_paths.hpp"
#include "sequences/contigs.hpp"
#include "vertex_resolution.hpp"
#include "common/fire_listeners.hpp"
#include "common/logging.hpp"
#include "alignment/alignment_form.hpp"
#include "path_replacement.hpp"
#include <memory>
#include <unordered_map>

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

//        A path P is replaced by a path Q spelling the same sequence (mergePath, flipVertex, extendVertex); see
//        PathReplacement for the geometry. Information of P's vertices that are not in Q goes to Q's dominant vertex,
//        information of P's edges to Q's edges, and that of the edges into P's start (out of P's finish), if Q starts
//        (ends) elsewhere, to their replacements. Codes concatenate, so paths through P keep their codes. Fired forward
//        and with both paths reverse-complemented (once if P is its own rc or contains a self-rc vertex that is removed).
        virtual void fireReplacePath(const PathReplacement &replacement) {}
        virtual void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {}
//        DBG only.
        virtual void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {}
        virtual void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &left_al, const AlignmentForm &right_al) {}
        virtual void fireSplitEdge(Edge &edge, const RAGraphPath &split) {}

        virtual void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {}
//        Checks the invariants of whatever this listener maintains, reporting every violation to os.
//        Run by AssemblyGraph::checkConsistency, so it may assume the graph structure itself is consistent.
        virtual bool checkConsistency(std::ostream &os) {return true;}
//        The graph this listener is attached to (AssemblyGraph is the only ResolutionFire).
        AssemblyGraph &getGraph();

        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {};
    };

    class ResolutionFire : public AbstractFire {
        std::vector<ResolutionListener *> listeners;
    public:
        ResolutionFire() = default;

        ResolutionFire(ResolutionFire &&other)   noexcept = default;
        ResolutionFire &operator=(ResolutionFire &&other)  noexcept = default;

//        Mirrored arguments are computed once per event, not once per listener.
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
            std::unique_ptr<VertexResolutionResult> rc;
            if(core != core.rc())
                rc = std::make_unique<VertexResolutionResult>(resolution.RC());
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireResolveVertex(core, resolution);
                if(rc)
                    listener->fireResolveVertex(core.rc(), *rc);
            }
        }

        void fireReplacePath(const RAGraphPath &old_path, const RAGraphPath &new_path, size_t dominant_index,
                             const std::unordered_map<EdgeId, EdgeId> &outer);

        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) {
            std::unique_ptr<GraphPath> rc;
            if(new_vertex != new_vertex.rc())
                rc = std::make_unique<GraphPath>(path.RC());
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireMergeLoop(path, new_vertex);
                if(rc)
                    listener->fireMergeLoop(*rc, new_vertex.rc());
            }
        }

        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {
            std::unique_ptr<RAGraphPath> rc;
            if(new_edge != new_edge.rc())
                rc = std::make_unique<RAGraphPath>(path.RC());
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireMergePathToEdge(path, new_edge);
                if(rc)
                    listener->fireMergePathToEdge(*rc, new_edge.rc());
            }
        }

        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right, const AlignmentForm &left_al, const AlignmentForm &right_al) {
            std::unique_ptr<AlignmentForm> rc_left, rc_right;
            if(new_edge != new_edge.rc()) {
                rc_left = std::make_unique<AlignmentForm>(right_al.RC());
                rc_right = std::make_unique<AlignmentForm>(left_al.RC());
            }
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireMergeTipsToEdge(new_edge, left, right, left_al, right_al);
                if(rc_left)
                    listener->fireMergeTipsToEdge(new_edge.rc(), right.rc(), left.rc(), *rc_left, *rc_right);
            }
        }

        void fireSplitEdge(Edge &edge, const RAGraphPath &split) {
            std::unique_ptr<RAGraphPath> rc;
            if(edge != edge.rc())
                rc = std::make_unique<RAGraphPath>(split.RC());
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                listener->fireSplitEdge(edge, split);
                if(rc)
                    listener->fireSplitEdge(edge.rc(), *rc);
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

        bool fireCheckConsistency(std::ostream &os) {
            bool ok = true;
            for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
                if (!listener->checkConsistency(os))
                    ok = false;
            }
            return ok;
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

        void fireReplacePath(const PathReplacement &replacement) override {
            lock();
            *outp << "FireReplacePath " << replacement.oldPath() << " -> " << replacement.newPath() << std::endl;
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
