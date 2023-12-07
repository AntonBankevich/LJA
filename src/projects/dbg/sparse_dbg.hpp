//
// Created by anton on 7/22/20.
//

#pragma once
#include "paths.hpp"
#include "assembly_graph.hpp"
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include "compact_path.hpp"
#include <common/oneline_utils.hpp>
#include <common/iterator_utils.hpp>
#include <common/object_id.hpp>
#include <utility>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <forward_list>

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

    class DBGVertexData {
    protected:
        std::list<Sequence> hanging{};
        hashing::htype hash;
    public:
        DBGVertexData(hashing::htype hash) : hash(hash) {}
        DBGVertexData RC() const {
            return *this;
        }
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
        mutable bool is_reliable = false;
        void incCov(int delta) {
#pragma omp atomic
            cov += delta;
        }
        size_t intCov() const {return cov;}
        double getCoverage() const {return double(cov) / truncSize();}
    };

    class DBGVertex : public ag::BaseVertex<DBGTraits>, public DBGVertexData {
    public:
        DBGVertex(id_type id, bool canonical, DBGVertexData data) : BaseVertex<DBGTraits>(id, canonical), DBGVertexData(std::move(data)) {}
        DBGVertex(id_type id, Sequence seq, DBGVertexData data) : BaseVertex<DBGTraits>(id, std::move(seq)), DBGVertexData(std::move(data)) {}
        const std::list<Sequence> &getHanging() const {return hanging;}
        void addOutgoingSequence(const Sequence &new_seq) {
            lock();
            for (Edge &edge : *this) {
                if (new_seq.size() < edge.truncSize() && edge.truncSeq().Subseq(new_seq.size()) == new_seq) {
                    unlock();
                    return;
                }
            }
            for (Sequence &out : hanging) {
                if (new_seq.size() > out.size() && new_seq.Subseq(0, out.size())== out) {
                    out = new_seq;
                    unlock();
                    return;
                } else if(new_seq.size() <= out.size() && out.Subseq(0, new_seq.size())== new_seq) {
                    unlock();
                    return;
                }
            }
            hanging.emplace_back(new_seq);
            unlock();
        }

        virtual void fireAddEdge(Edge &edge) override {
            for(auto it = hanging.begin(); it != hanging.end(); ++it) {
                if(it->size() <= edge.truncSize() && edge.truncSeq().Subseq(0, it->size()) == *it) {
                    hanging.erase(it);
                    break;
                }
            }
        }

        void clearHanging() {hanging.clear();}
        hashing::htype getHash() const {return hash;}
    };

    typedef DBGEdge Edge;
    typedef DBGVertex Vertex;
    typedef ag::EdgePosition<DBGTraits> EdgePosition;
    typedef DBGEdge::EdgeId EdgeId;
    typedef DBGVertex::VertexId VertexId;
    typedef DBGEdge::ConstEdgeId ConstEdgeId;
    typedef DBGVertex::ConstVertexId ConstVertexId;
    typedef ag::GraphPath<dbg::DBGTraits> GraphPath;
    typedef ag::CompactPath<dbg::DBGTraits> CompactPath;

    class SparseDBG : public ag::AssemblyGraph<DBGTraits> {
    private:
        hashing::RollingHash hasher_;
    public:
        explicit SparseDBG(const hashing::RollingHash &_hasher) : hasher_(_hasher) {}
        SparseDBG(SparseDBG &&other) = default;
        SparseDBG &operator=(SparseDBG &&other) = default;
        SparseDBG(const SparseDBG &other) noexcept = delete;
        template <class I>
        SparseDBG(I begin, I end, const hashing::RollingHash &_hasher) : hasher_(_hasher) {
            std::vector<hashing::htype> all(begin, end);
            std::sort(all.begin(), all.end());
            all.erase(std::unique(all.begin(), all.end()), all.end());
            for(hashing::htype hash : all) {
                addKmerVertex(hash);
            }
        }

        const hashing::RollingHash &hasher() const {return hasher_;}

        Vertex &addKmerVertex(const Sequence &kmer, Vertex::id_type id = 0) {
            return AssemblyGraph<DBGTraits>::addVertex(kmer, DBGVertexData(hashing::KWH(hasher_, kmer, 0).hash()), id);
        }

        Vertex &addKmerVertex(const hashing::KWH &kwh, Vertex::id_type id = 0) {
            return AssemblyGraph<DBGTraits>::addVertex(kwh.getSeq(), VertexData(kwh.hash()), id);
        }

        Vertex &addKmerVertex(hashing::htype hash, Vertex::id_type id = 0) {
            return AssemblyGraph<DBGTraits>::addVertexPair(VertexData(hash), id);
        }
    };


}
