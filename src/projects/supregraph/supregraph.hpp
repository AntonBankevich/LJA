#pragma once
#include "vertex_resolution.hpp"
#include "supregraph_base.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "sequences/contigs.hpp"
#include "assembly_graph/paths.hpp"
#include "assembly_graph/compact_path.hpp"

namespace spg {

    class ResolutionListener {
    public:
        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) = 0;
        virtual ~ResolutionListener() = default;
    };

    class ResolutionFire {
        std::vector<ResolutionListener *> listeners;
    public:
        void addListener(ResolutionListener &listener) {
            listeners.emplace_back(&listener);
        }
        void removeListener(ResolutionListener &listener) {
            listeners.erase(std::find(listeners.begin(), listeners.end(), &listener));
        }
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {
            for(auto *listener:listeners) {
                listener->fireResolveVertex(core, resolution);
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
        std::vector<VertexId> resolveVertex(Vertex &core, const VertexResolutionPlan &resolution) {
            VERIFY(core.isCore() && core.inDeg() > 0 && core.outDeg() > 0);
            VertexResolutionResult result(core);
            std::vector<VertexId> res;
            for(const EdgePair &p : resolution.connectionsUnique()) {
                VERIFY(p.first->getFinish() == core);
                if(p.first->getFinish() != core)
                    continue;
                VERIFY(p.second->getStart() == core);
                VERIFY(p.first->isSuffix());
                VERIFY(p.second->isPrefix());
                Sequence seq = p.first->rc().truncSeq().rc() + core.getSeq() + p.second->truncSeq();
                if(core == core.rc() && seq > !seq)
                    continue;
                Vertex &newv = addSPGVertex(seq, false, false, false);
                result.add(newv, p);
                p.first->getStart().addSPEdgeLockFree(newv);
                newv.addSPEdgeLockFree(p.second->getFinish());
                res.emplace_back(newv.getId());
                if(newv != newv.rc())
                    res.emplace_back(newv.rc().getId());
            }
            fireResolveVertex(core, result);
            IsolateAndMark(core);
            return std::move(res);
        }
    };

    template<class Traits>
    class SPGConverter {
    private:
        typedef ag::AssemblyGraph<Traits> OGraph;
        typedef typename Traits::Edge OEdge;
        typedef typename Traits::Vertex OVertex;
        typedef typename OVertex::VertexId OVertexId;
        typedef typename OEdge::EdgeId OEdgeId;
        typedef ag::GraphPath<Traits> OPath;
        std::unordered_map<OVertexId, Segment<Vertex>> vmap;

        OPath maxExtension(OVertex &vertex) const;
        VertexId processUnbranching(SupreGraph &g, const OVertex &v) const;
        void constructVertices(OGraph &other, SupreGraph &g);
        void addSPEdges(OGraph &other) const;
        void addOuterEdges(OGraph &other);

    public:
        SPGConverter() = default;

        //        Before using this method merge unbranching paths.
//        This method only works for graphs with edges distinguished by their first letters.
//        TODO: rewrite this using a sequence of operations: outerEdgeToVertex, normalizeLoops, extendSequences
        SupreGraph Convert(OGraph &other);

//        This method does not allow to convert loops paths properly
        GraphPath convertPath(OPath &path) const;
    };

    template<class Traits>
    typename SPGConverter<Traits>::OPath SPGConverter<Traits>::maxExtension(OVertex &vertex) const {
        OPath path(vertex);
        while(path.finish().outDeg() == 1) {
            path += path.finish().front();
            if(path.finish() == path.start() || path.finish() == path.getVertex(path.size() / 2))
                break;
        }
        if(path.empty())
            return std::move(path);
        size_t pos = path.size() - 1;
        while(pos > 0 && path.getVertex(pos) != path.finish())
            pos--;
        while(pos > 0 && path.backEdge() == path.getEdge(pos - 1)) {
            path.pop_back();
            pos--;
        }
        return std::move(path);
    }

    template<class Traits>
    VertexId SPGConverter<Traits>::processUnbranching(SupreGraph &g, const OVertex &v) const {
        OPath fpath = OPath::WalkForward(v.front());
        VERIFY(fpath.finish() == v || fpath.finish() == v.rc());
        Sequence loop;
        if(fpath.finish() == v && v != v.rc()) {
            if(v < v.rc())
                return {};
            loop = v.front().truncSeq();
        } else {
            if(v == v.rc()) {
                loop = v.front().getSeq();
                VERIFY(loop == !loop);
                loop = loop.Subseq(v.size() / 2, loop.size() / 2);
            } else {
                loop = v.rc().front().getSeq() + v.front().truncSeq() + v.rc().front().truncSeq();
                VERIFY(loop == !loop);
                loop = loop.Subseq(v.rc().front().fullSize() / 2, loop.size() / 2);
            }
            loop = !loop + loop;
            VERIFY(v.rc().front() == v.rc().front().rc())
        }
        Vertex &newv = g.addSPGVertex(loop, true, false, false, v.getInnerId());
        newv.addSPEdgeLockFree(newv);
        return newv.getId();
    }

    template<class Traits>
    void SPGConverter<Traits>::constructVertices(SPGConverter::OGraph &other, SupreGraph &g) {
        for(OVertex &v : other.vertices()) {
            if(v.inDeg() == 1) {
                if(v.outDeg() == 1) processUnbranching(g, v);
            } else {
                OPath right = maxExtension(v);
                std::cout << right.str() << std::endl;
                OVertex * f = &right.finish().rc();
                if((right.empty() || (right.size() == 1 && right.finish().inDeg() == 1)) && vmap.find(right.finish().rc().getId()) != vmap.end()) {
//                        Only possible for core vertices
                    OVertexId rv = right.finish().rc().getId();
                    vmap[v.getId()] = {vmap[rv].contig().rc(), 0, v.size()};
                    VERIFY_MSG(v.getSeq() == vmap[v.getId()].fullSeq(), "1");
                    vmap[v.rc().getId()] = vmap[v.getId()].RC();
                    VERIFY_MSG(v.rc().getSeq() == vmap[v.rc().getId()].fullSeq(), "2");
                } else {
                    bool inf_right = false;
                    for(size_t i = 0; i < right.size(); i++)
                        if(right.getVertex(i) == right.finish()) {
                            inf_right = true;
                            break;
                        }
                    Sequence vseq = right.Seq();
                    if(inf_right)
                        vseq = vseq.Subseq(0, vseq.size() - right.finish().size());
                    Vertex &newv = g.addSPGVertex(vseq, false, false, inf_right, v.getInnerId());
                    vmap[v.getId()] = {newv, 0, v.size()};
                    VERIFY_MSG(v.getSeq() == vmap[v.getId()].fullSeq(), "3");
                    vmap[v.rc().getId()] = vmap[v.getId()].RC();
                    VERIFY_MSG(v.rc().getSeq() == vmap[v.rc().getId()].fullSeq(), "4");
                }
            }
        }
    }

    template<class Traits>
    void SPGConverter<Traits>::addSPEdges(SPGConverter::OGraph &other) const {
        for(OVertex &v : other.vertices()) {
            if(v.inDeg() == 1 || v.outDeg() != 1)
                continue;
            VERIFY(vmap.find(v.getId()) != vmap.end());
            OVertex &u = v.front().getFinish();
            if(u.inDeg() == 1)
                continue;
            VERIFY(vmap.find(u.getId()) != vmap.end());
            Vertex& from = vmap.at(v.getId()).contig();
            Vertex& to = vmap.at(u.getId()).contig();
            if(from.isInfRight()) {
                VERIFY(to.isInfRight());
                size_t shift = v.front().rc().truncSize();
                VERIFY(to.getSeq().startsWith(from.getSeq().Subseq(shift)));
                Sequence seq = from.getSeq().Subseq(0, shift) + to.getSeq();
                from.addEdgeLockFree(to, seq);
            } else
                from.addSPEdgeLockFree(to);
        }
    }

    template<class Traits>
    void SPGConverter<Traits>::addOuterEdges(SPGConverter::OGraph &other) {
        for(OEdge &e : other.edgesUnique()) {
            if(e.getStart().outDeg() != 1 && e.getFinish().inDeg() != 1) {
                VERIFY(vmap.find(e.getStart().rc().getId()) != vmap.end());
                VERIFY(vmap.find(e.getFinish().getId()) != vmap.end());
                Vertex& rcfrom = vmap[e.getStart().rc().getId()].contig();
                Vertex& from = rcfrom.rc();
                Vertex& to = vmap[e.getFinish().getId()].contig();
                Sequence seq = from.getSeq() + e.truncSeq() + to.getSeq().Subseq(e.getFinish().size());
                Edge &new_edge = from.addEdgeLockFree(to, seq);
            }
        }
    }

    template<class Traits>
    SupreGraph SPGConverter<Traits>::Convert(SPGConverter::OGraph &other) {
        SupreGraph g;
        constructVertices(other, g);
        for(Vertex &v : g.vertices()) {
            VERIFY(!v.getSeq().empty());
        }
        addSPEdges(other);
        addOuterEdges(other);
        return std::move(g);
    }

    template<class Traits>
    GraphPath SPGConverter<Traits>::convertPath(SPGConverter::OPath &path) const {
        if(!path.valid())
            return {};
        if(!path.start().isJunction())
            return {};
        Segment<Vertex> seg = vmap.at(path.start().getId());
        GraphPath res(seg.contig(), seg.left, seg.contig().size() - seg.right);
        if(path.size() == 0)
            return std::move(res);
        VERIFY(path.start().getSeq() == res.Seq());
        OPath opath;
        for(OEdge &edge : path.edges()) {
            opath += edge;
            size_t pos = 0;
            while(pos < edge.truncSize()) {
                if(res.rightSkip() != 0) {
                    size_t l = std::min(edge.truncSize() - pos, res.rightSkip());
                    pos += l;
                    res.uniqueExtendBack(l);
                } else {
                    while(res.finish().outDeg() == 1 && res.finish().front().isSuffix()) {
                        res += res.finish().front();
                    }
                    res += res.finish().getOutgoing(edge.truncSeq()[pos]);
                    pos += 1;
                    res.cutBack(res.backEdge().truncSize() - 1);
                }
            }
            VERIFY(opath.Seq() == res.Seq());
        }
        size_t l = path.leftSkip() + res.leftSkip();
        size_t r = path.rightSkip() + res.rightSkip();
        res.uniqueExtendBack(res.rightSkip());
        res.uniqueExtendFront(res.leftSkip());
        while(res.finish().outDeg() == 1 && res.finish().front().truncSize() == 0)
            res += res.finish().front();
        res = res.RC();
        while(res.finish().outDeg() == 1 && res.finish().front().truncSize() == 0)
            res += res.finish().front();
        res = res.RC();
        res.cutFront(l);
        res.cutBack(r);
        return std::move(res);
    }
};