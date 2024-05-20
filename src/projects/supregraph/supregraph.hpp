#pragma once
#include "vertex_resolution.hpp"
#include "supregraph_base.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "sequences/contigs.hpp"
#include "assembly_graph/paths.hpp"
#include "assembly_graph/compact_path.hpp"

namespace spg {

    class ResolutionFire;
    class ResolutionListener {
    private:
        ResolutionFire *fire;
    public:
        ResolutionListener(ResolutionFire &fire);
        ResolutionListener(ResolutionListener &&other)  noexcept;
        ResolutionListener(const ResolutionListener &other) = delete;

        virtual ~ResolutionListener();
        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) = 0;
        virtual void fireDeleteVertex(Vertex &v) {}
        virtual void fireDeleteEdge(Edge &e) {}
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
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
            for(auto *listener:listeners) {
                listener->fireResolveVertex(core, resolution);
            }
        }
        void fireDeleteVertex(Vertex &vertex) {
            for(auto *listener:listeners) {
                listener->fireDeleteVertex(vertex);
            }
        }

        void fireDeleteEdge(Edge &edge) {
            for(auto *listener:listeners) {
                listener->fireDeleteEdge(edge);
            }
        }
    };

    class SupreGraph : public ag::AssemblyGraph<SPGTraits>, public ResolutionFire {
    public:
        SupreGraph() = default;
        SupreGraph(SupreGraph &&other) = default;
        SupreGraph &operator=(SupreGraph &&other) = default;
        SupreGraph(const SupreGraph &) = delete;
        Vertex &addSPGVertex(Sequence seq, bool cyclic, bool inf_left, bool inf_right, Vertex::id_type id = Vertex::id_type()) {
            return addVertex(std::move(seq), SPGVertexData(cyclic, inf_left, inf_right), id);
        }

        Vertex &outerEdgeToVertex(Edge &edge);

        void IsolateAndMark(Vertex &v);

//        Graph should be in normal form. resolution should contain all new edges to be created including reverse-complement
//Returns all new vertices including rc
        VertexResolutionResult resolveVertex(Vertex &core, const VertexResolutionPlan &resolution) {
            VERIFY(core.isCore() && core.inDeg() > 0 && core.outDeg() > 0);
            VertexResolutionResult result(core);
            for(const EdgePair &p : resolution.connectionsUnique()) {
                VERIFY(p.first->getFinish() == core);
                if(p.middle() != core)
                    continue;
                VERIFY(p.second->getStart() == core);
                VERIFY(p.first->isSuffix());
                VERIFY(p.second->isPrefix());
                Sequence seq = p.getSeq();
                Vertex &newv = addSPGVertex(seq, false, false, false);
                result.add(newv, p);
                p.first->getStart().addSPEdgeLockFree(newv);
                newv.addSPEdgeLockFree(p.second->getFinish());
            }
            fireResolveVertex(core, result);
            if(core != core.rc())
                fireResolveVertex(core.rc(), result.RC());
            for(Edge &edge : core) {
                fireDeleteEdge(edge);
            }
            if(core != core.rc())
                for(Edge &edge : core.incoming()) {
                    fireDeleteEdge(edge);
                }
            for(Edge &edge : core.rc()) {
                fireDeleteEdge(edge);
            }
            if(core != core.rc())
                for(Edge &edge : core.rc().incoming()) {
                    fireDeleteEdge(edge);
                }
            fireDeleteVertex(core);
            if(core != core.rc())
                fireDeleteVertex(core.rc());
            IsolateAndMark(core);
            return std::move(result);
        }
    };

}