//
// Created by anton on 7/22/20.
//

#pragma once
#include "assembly_graph/data_structures/component.hpp"
#include "assembly_graph/random_access_paths.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include <common/oneline_utils.hpp>
#include <common/iterator_utils.hpp>
#include <common/object_id.hpp>
#include <utility>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <forward_list>
#include <assembly_graph/ag_algorithms.hpp>

namespace dbg {

//    TODO: this class should be constructed as a mixture of multiple classes each representing possible piece of informations
//that is to be stored in the edge. Corresponding information should be able to support itself during various graph
//operations.
    class DBGEdgeData {
    protected:
        size_t cov = 0;
    public:
        DBGEdgeData RC() const {
            return *this;
        }

        template<class I>
        static DBGEdgeData Merge(I begin, I end) {
            size_t cov = 0;
            for(;begin != end; ++begin) {
                cov += begin->cov;
            }
            DBGEdgeData res;
            res.cov = cov;
            return std::move(res);
        }

    };

    class HashListener;
    class DBGVertexData {
        friend class HashListener;
    protected:
        std::list<Sequence> hanging{};
        hashing::htype hash;
    public:
        static const hashing::htype default_hash;
        DBGVertexData(hashing::htype hash = default_hash) : hash(hash) {}
        DBGVertexData RC() const {
            return *this;
        }
        hashing::htype getHash() const {return hash;}
    };

    class DBGVertex;
    class DBGEdge;

    struct DBGTraits {
        typedef DBGVertexData VertexData;
        typedef DBGEdgeData EdgeData;
        typedef DBGVertex Vertex;
        typedef DBGEdge Edge;
    };


    class DBGEdge : public ag::BaseEdge<DBGTraits>, public DBGEdgeData {
    public:
        DBGEdge(id_type id, Vertex &_start, Vertex &_end, Sequence _seq, DBGEdgeData data) :
                BaseEdge<DBGTraits>(id, _start, _end, std::move(_seq)), DBGEdgeData(std::move(data)) {}
        DBGEdge() {
//            TODO: Remove this!!! It exists only for Andreys code compilation but that code should be purged
        }
        DBGEdge(DBGEdge &&) = delete;
        DBGEdge(const DBGEdge &) = delete;
        mutable bool is_reliable = false;
        void incCov(int64_t delta) {
#pragma omp atomic
            cov += delta;
            VERIFY(cov < size_t(-1) >> 2)
        }
        size_t intCov() const {return cov;}
        void setCov(size_t val) {cov = val;}
        double getCoverage() const {return double(cov) / truncSize();}
    };

    class DBGVertex : public ag::BaseVertex<DBGTraits>, public DBGVertexData {
    public:
        DBGVertex(id_type id, bool canonical, DBGVertexData data) : BaseVertex<DBGTraits>(id, canonical), DBGVertexData(std::move(data)) {}
        DBGVertex(id_type id, Sequence seq, DBGVertexData data) : BaseVertex<DBGTraits>(id, std::move(seq)), DBGVertexData(std::move(data)) {}
    };

//    class DBGMaintainence : public ag::ResolutionListener<DBGTraits> {
//    private:
//        hashing::RollingHash hasher;
//    public:
//        DBGMaintainence(ag::ResolutionFire<DBGTraits> &dbg, hashing::RollingHash &hasher) :
//                                ag::ResolutionListener<DBGTraits>(dbg), hasher(hasher) {}
//        void fireAddVertex(Vertex &v) override {
//            VERIFY(!v.getSeq().empty() || v.getHash() != Vertex::default_hash);
//            if(v.getHash() == Vertex::default_hash) {
//                v.hash = hasher.hash(v.getSeq(), 0);
//            }
//        }
//        void fireMergePath(const std::vector<EdgeId> &path, Vertex &new_vertex) override {VERIFY(false);}
//        void fireMergeLoop(const ag::GraphPath <DBGTraits> &path, Vertex &new_vertex) override {VERIFY(false);}
//        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult<DBGTraits> &resolution) override {VERIFY(false);};
//        void fireAddSupreVertex(Vertex &v, Edge &e) override {VERIFY(false);}
//    };

    class HashListener : ag::ResolutionListener<DBGTraits> {
    private:
        hashing::RollingHash hasher;
    public:
        HashListener(ag::ResolutionFire<DBGTraits> &fire, const hashing::RollingHash &hasher) :
                            ag::ResolutionListener<DBGTraits>(fire, "HashListener"), hasher(hasher) {}

        void fireAddVertex(Vertex &v) override {
            if(v.getHash() == DBGVertexData::default_hash && !v.getSeq().empty()) {
                v.hash = hashing::MovingKWH(hasher, v.getSeq(), 0).hash();
            }
        }
    };

    typedef DBGEdge Edge;
    typedef DBGVertex Vertex;
    typedef ag::EdgePosition<DBGTraits> EdgePosition;
    typedef DBGEdge::EdgeId EdgeId;
    typedef DBGVertex::VertexId VertexId;
    typedef DBGEdge::ConstEdgeId ConstEdgeId;
    typedef DBGVertex::ConstVertexId ConstVertexId;
    typedef ag::GraphPath<DBGTraits> GraphPath;
    typedef ag::Component<DBGTraits> Component;
    typedef ag::PathHelper<DBGTraits> PathHelper;
    typedef ag::PathPosition<DBGTraits> PathPosition;

    class SparseDBG : public ag::AssemblyGraph<DBGTraits> {
    private:
        hashing::RollingHash hasher_;
        HashListener hashListener;
    public:
        explicit SparseDBG(const hashing::RollingHash &_hasher) : hasher_(_hasher), hashListener(*this, _hasher) {}
        SparseDBG(SparseDBG &&other) = default;
        SparseDBG &operator=(SparseDBG &&other) = default;
        SparseDBG(const SparseDBG &other) noexcept = delete;
        template <class I>
        SparseDBG(I begin, I end, const hashing::RollingHash &_hasher) : hasher_(_hasher), hashListener(*this, _hasher) {
            std::vector<hashing::htype> all(begin, end);
            std::sort(all.begin(), all.end());
            all.erase(std::unique(all.begin(), all.end()), all.end());
            for(hashing::htype hash : all) {
                addKmerVertex(hash);
            }
        }

        const hashing::RollingHash &hasher() const {return hasher_;}
        size_t getK() const {return hasher().getK();}

        Vertex &addKmerVertex(const hashing::KWH &kwh, Vertex::id_type id = 0) {
            return AssemblyGraph<DBGTraits>::addVertex(kwh.getSeq(getK()), VertexData(kwh.hash()), id);
        }

        Vertex &addKmerVertex(const Sequence &kmer, Vertex::id_type id = 0) {
            return addKmerVertex(hashing::MovingKWH(hasher_, kmer, 0), id);
        }

        Vertex &addKmerVertex(hashing::htype hash, Vertex::id_type id = 0) {
            return AssemblyGraph<DBGTraits>::addVertexPair(VertexData(hash), id);
        }
    };

    inline std::string SaveEdgeName(const Edge &edge) {
        VERIFY((edge.getFinish().rc().getInnerId() > 0) == edge.getFinish().rc().isCanonical());
        return edge.getInnerId().str() + "_" + edge.rc().getInnerId().str();
    }

}
