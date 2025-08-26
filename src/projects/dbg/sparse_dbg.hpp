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
    using ag::Vertex;
    using ag::Edge;
    using ag::VertexId;
    using ag::EdgeId;
    using ag::VertexData;
    using ag::EdgeData;

    class SparseDBG;

    class HashListener : public ag::ResolutionListener {
    private:
        hashing::RollingHash hasher;
    public:
        HashListener(SparseDBG &dbg, const hashing::RollingHash &hasher);

        void fireAddVertex(ag::Vertex &v) override;
    };

    class SparseDBG : public ag::AssemblyGraph {
        friend class HashListener;
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

        void disableHashing() {hashListener.detach();}

        const hashing::RollingHash &hasher() const {return hasher_;}
        size_t getK() const {return hasher().getK();}

        Vertex &addKmerVertex(const hashing::KWH &kwh, Vertex::id_type id = 0) {
            Vertex &res = AssemblyGraph::addVertex(kwh.getSeq(getK()), VertexData::DBGData(kwh.hash()), id);
            return res;
        }

        Vertex &addKmerVertex(const Sequence &kmer, Vertex::id_type id = 0) {
            return addKmerVertex(hashing::MovingKWH(hasher_, kmer, 0), id);
        }

        Vertex &addKmerVertex(hashing::htype hash, Vertex::id_type id = 0) {
            return this->addVertexPair(VertexData::DBGData(hash), id);
        }
    };

}
