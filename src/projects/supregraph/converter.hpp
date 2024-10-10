#pragma once

#include <dbg/graph_alignment_storage.hpp>
#include "vertex_resolution.hpp"
#include "supregraph_base.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "sequences/contigs.hpp"
#include "assembly_graph/paths.hpp"
#include "assembly_graph/compact_path.hpp"
#include "read_storage.hpp"

namespace spg {
    template<class Traits>
    class SPGConverter {
    private:
        typedef ag::AssemblyGraph<Traits> OGraph;
        typedef typename Traits::Edge OEdge;
        typedef typename Traits::Vertex OVertex;
        typedef typename OVertex::VertexId OVertexId;
        typedef typename OEdge::EdgeId OEdgeId;
        typedef ag::GraphPath<Traits> OPath;
        std::unordered_map<OVertexId, Segment<spg::Vertex>> vmap;
        std::unordered_map<OEdgeId, Segment<spg::Vertex>> outer_map;

        OPath maxExtension(OVertex &vertex) const;
        spg::VertexId processUnbranching(spg::SupreGraph &g, const OVertex &v) const;
        void constructInnerVertices(OGraph &other, spg::SupreGraph &g);
        void constructOuterVertices(SPGConverter::OGraph &other, spg::SupreGraph &g);
        void addSPEdges(spg::SupreGraph &g, OGraph &other) const;

    public:
        SPGConverter() = default;

        //        Before using this method merge unbranching paths.
//        This method only works for graphs with edges distinguished by their first letters.
//        TODO: rewrite this using a sequence of operations: outerEdgeToVertex, normalizeLoops, extendSequences
        spg::SupreGraph convert(OGraph &other);

//        This method does not allow to convert loops paths properly
        spg::GraphPath convertPath(const OPath &path) const;

        spg::PathStorage convertLibs(const std::vector<dbg::ReadAlignmentStorage *> &storages, SupreGraph &spg);

        Vertex &map(OVertex &v) const {
            return vmap.at(v.getId()).contig();
        }
        Vertex &map(OEdge &e) const {
            VERIFY(e.getStart().outDeg() > 1 && e.getFinish().inDeg() > 1);
            return outer_map.at(e.getId()).contig();
        }
    };

    template<class Traits>
    typename SPGConverter<Traits>::OPath SPGConverter<Traits>::maxExtension(OVertex &vertex) const {
        OPath path(vertex);
        while (path.finish().outDeg() == 1) {
            path += path.finish().front();
            if (path.finish() == path.start() || path.finish() == path.getVertex(path.size() / 2))
                break;
        }
        if (path.empty())
            return std::move(path);
        size_t pos = path.size() - 1;
        while (pos > 0 && path.getVertex(pos) != path.finish())
            pos--;
        while (pos > 0 && path.backEdge() == path.getEdge(pos - 1)) {
            path.pop_back();
            pos--;
        }
        return std::move(path);
    }

    template<class Traits>
    spg::VertexId SPGConverter<Traits>::processUnbranching(spg::SupreGraph &g, const OVertex &v) const {
        OPath fpath = OPath::WalkForward(v.front());
        VERIFY(fpath.finish() == v || fpath.finish() == v.rc());
        Sequence loop;
        if (fpath.finish() == v && v != v.rc()) {
            if (v < v.rc())
                return {};
            loop = v.front().truncSeq();
        } else {
            if (v == v.rc()) {
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
        spg::Vertex &newv = g.addSPGVertex(loop, true, false, false);
        g.addSPEdgeLockFree(newv, newv);
        return newv.getId();
    }

    template<class Traits>
    void SPGConverter<Traits>::constructInnerVertices(SPGConverter::OGraph &other, spg::SupreGraph &g) {
        for (OVertex &v: other.vertices()) {
            if (v.inDeg() == 1) {
                if (v.outDeg() == 1) processUnbranching(g, v);
            } else {
                OPath right = maxExtension(v);
                OVertex *f = &right.finish().rc();
                if ((right.empty() || (right.size() == 1 && right.finish().inDeg() == 1)) &&
                    vmap.find(right.finish().rc().getId()) != vmap.end()) {
//                        Only possible for core vertices
                    OVertexId rv = right.finish().rc().getId();
                    vmap[v.getId()] = {vmap[rv].contig().rc(), 0, v.size()};
                    VERIFY(v.getSeq() == vmap[v.getId()].fullSeq());
                    vmap[v.rc().getId()] = vmap[v.getId()].RC();
                    VERIFY(v.rc().getSeq() == vmap[v.rc().getId()].fullSeq());
                } else {
                    bool inf_right = false;
                    for (size_t i = 0; i < right.size(); i++)
                        if (right.getVertex(i) == right.finish()) {
                            inf_right = true;
                            break;
                        }
                    Sequence vseq = right.Seq();
                    if (inf_right)
                        vseq = vseq.Subseq(0, vseq.size() - right.finish().size());
                    spg::Vertex &newv = g.addSPGVertex(vseq, false, false, inf_right);
                    vmap[v.getId()] = {newv, 0, v.size()};
                    VERIFY(v.getSeq() == vmap[v.getId()].fullSeq());
                    vmap[v.rc().getId()] = vmap[v.getId()].RC();
                    VERIFY(v.rc().getSeq() == vmap[v.rc().getId()].fullSeq());
                }
            }
        }
    }

    template<class Traits>
    void SPGConverter<Traits>::constructOuterVertices(SPGConverter::OGraph &other, spg::SupreGraph &g) {
        for (OEdge &e: other.edgesUnique()) {
            if(!e.isOuter())
                continue;
            OPath left = maxExtension(e.getStart().rc());
            OPath right = maxExtension(e.getFinish());
            OPath vPath = left.RC() + e + right;
            Sequence vseq = vPath.Seq();
            spg::Vertex &newv = g.addSPGVertex(vseq, false, false, false);
            outer_map[e.getId()] = {newv, left.truncLen(), newv.size() - right.truncLen()};
        }
    }

    template<class Traits>
    void SPGConverter<Traits>::addSPEdges(spg::SupreGraph &g, SPGConverter::OGraph &other) const {
        for (OVertex &v: other.vertices()) {
            if (v.inDeg() == 1 || v.outDeg() != 1)
                continue;
            VERIFY(vmap.find(v.getId()) != vmap.end());
            OVertex &u = v.front().getFinish();
            if (u.inDeg() == 1)
                continue;
            VERIFY(vmap.find(u.getId()) != vmap.end());
            spg::Vertex &from = vmap.at(v.getId()).contig();
            spg::Vertex &to = vmap.at(u.getId()).contig();
            if (from.isInfRight()) {
                VERIFY(to.isInfRight());
                size_t shift = v.front().rc().truncSize();
                VERIFY(to.getSeq().startsWith(from.getSeq().Subseq(shift)));
                Sequence seq = from.getSeq().Subseq(0, shift) + to.getSeq();
                g.addEdgeLockFree(from, to, seq);
            } else
                g.addSPEdgeLockFree(from, to);
        }
        for (OEdge &e: other.edgesUnique()) {
            if(!e.isOuter())
                continue;
            Vertex &outer = outer_map.at(e.getId()).contig();
            Vertex &start = vmap.at(e.getStart().getId()).contig();
            Vertex &finish = vmap.at(e.getFinish().getId()).contig();
            g.addSPEdgeLockFree(start, outer);
            if(start != start.rc())
                g.addSPEdgeLockFree(outer, finish);
        }
    }

    template<class Traits>
    spg::SupreGraph SPGConverter<Traits>::convert(SPGConverter::OGraph &other) {
        spg::SupreGraph g;
        constructInnerVertices(other, g);
        constructOuterVertices(other, g);
        for (spg::Vertex &v: g.vertices()) {
            VERIFY(!v.getSeq().empty());
        }
        addSPEdges(g, other);
        return std::move(g);
    }

    template<class Traits>
    spg::GraphPath SPGConverter<Traits>::convertPath(const SPGConverter::OPath &path) const {
        if (!path.valid())
            return {};
        if (!path.start().isJunction())
            return {};
        Segment<spg::Vertex> seg = vmap.at(path.start().getId());
        spg::GraphPath res(seg.contig(), seg.left, seg.contig().size() - seg.right);
        if (path.size() == 0)
            return std::move(res);
        VERIFY(path.start().getSeq() == res.Seq());
//        OPath opath(path.start());
//        VERIFY(opath.Seq() == res.Seq());
        for (OEdge &edge: path.edges()) {
//            opath += edge;
            res.fastExtend(edge.truncSeq());
//            VERIFY(opath.Seq() == res.Seq());
        }
        res.cutFront(path.cutLeft());
        res.cutBack(path.cutRight());
        res.normalize();
//
//        size_t l = path.cutLeft() + res.cutLeft();
//        size_t r = path.cutRight()() + res.cutRight()();
//        res.uniqueExtendBack(res.cutRight()());
//        res.uniqueExtendFront(res.cutLeft());
//        while(res.finish().outDeg() == 1 && res.finish().front().truncSize() == 0)
//            res += res.finish().front();
//        res = res.RC();
//        while(res.finish().outDeg() == 1 && res.finish().front().truncSize() == 0)
//            res += res.finish().front();
//        res = res.RC();
//        res.cutFront(l);
//        res.cutBack(r);
        return std::move(res);
    }

    template<class Traits>
    spg::PathStorage
    SPGConverter<Traits>::convertLibs(const std::vector<dbg::ReadAlignmentStorage *> &storages, SupreGraph &spg) {
        spg::PathStorage res(spg);
        for (const dbg::ReadAlignmentStorage *storageIt: storages) {
            const dbg::ReadAlignmentStorage &storage = *storageIt;
            for (const ag::AlignedRead<dbg::DBGTraits> &alignedRead : storage) {
                convertPath(alignedRead.path.unpack());
                res.addRead(alignedRead.id, convertPath(alignedRead.path.unpack()));
            }
        }
        return std::move(res);
    }
}