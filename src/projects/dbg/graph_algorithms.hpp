#pragma once
#include "dbg_graph_aligner.hpp"
#include "assembly_graph/compact_path.hpp"
#include "assembly_graph/paths.hpp"
#include "sparse_dbg.hpp"
namespace dbg {

    class DbgConstructionHelper {
    private:
        hashing::RollingHash hash;
        const hashing::RollingHash &hasher() const {
            return hash;
        }
    public:
        DbgConstructionHelper(const hashing::RollingHash &hash) : hash(hash){
        }
        void checkConsistency(size_t threads, logging::Logger &logger, SparseDBG &dbg) const;
        void checkDBGConsistency(size_t threads, logging::Logger &logger, SparseDBG &dbg) const;
        void checkSeqFilled(size_t threads, logging::Logger &logger, SparseDBG &dbg) const;

        void processRead(SparseDBG &dbg, KmerIndex &index, const Sequence &seq) const;
        void processFullEdgeSequence(SparseDBG &dbg, KmerIndex &index, const Sequence &old_seq) const;

        SparseDBG Subgraph(std::vector<Segment<Edge>> &pieces) const;

        void addAllKmers(SparseDBG &dbg, const std::vector<Sequence> &new_seqs, KmerIndex &index) const;

        void checkIndexConsistency(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex& index) const;
    };

    template<class Iterator>
    void FillSparseDBGEdges(SparseDBG &sdbg, Iterator begin, Iterator end, logging::Logger &logger, size_t threads,
                            const size_t min_read_size) {
        typedef typename Iterator::value_type ContigType;
        logger.trace() << "Starting to fill edges" << std::endl;
        KmerIndex index(sdbg);
        DbgConstructionHelper helper(sdbg.hasher());
        std::function<void(size_t, ContigType &)> task = [&sdbg, &helper, min_read_size, &index](size_t pos, ContigType &contig) {
            Sequence seq = contig.makeSequence();
            if (seq.size() >= min_read_size) {
                helper.processRead(sdbg, index, seq);
            }
        };
        processRecords(begin, end, logger, threads, task);
        logger.trace() << "Sparse graph edges filled." << std::endl;
    }

    template<class Iterator>
    void RefillSparseDBGEdges(logging::Logger &logger, size_t threads, SparseDBG &sdbg, Iterator begin, Iterator end, KmerIndex &index) {
        logger.trace() << "Starting to fill edges" << std::endl;
        std::function<void(size_t, std::pair<Vertex *, Sequence> &)> task = [&sdbg, &index](size_t pos,
                                                                                    std::pair<Vertex *, Sequence> &contig) {
            Sequence seq = contig.first->getSeq() + contig.second;
            DbgConstructionHelper(sdbg.hasher()).processFullEdgeSequence(sdbg, index, seq);
        };
        processObjects(begin, end, logger, threads, task);
        logger.trace() << "Sparse graph edges filled." << std::endl;
    }

    SparseDBG
    LoadDBGFromEdgeSequences(logging::Logger &logger, size_t threads, const io::Library &lib, hashing::RollingHash &hasher);

    template<class Iterator>
    void fillCoverage(SparseDBG &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                      const hashing::RollingHash &hasher, size_t min_read_size);

    SparseDBG constructSparseDBGFromReads(logging::Logger &logger, const io::Library &reads_file, size_t threads,
                                          const hashing::RollingHash &hasher,
                                          const std::vector<hashing::htype> &hash_list, size_t w);

    void tieTips(logging::Logger &logger, SparseDBG &sdbg, size_t k, size_t w, size_t threads);

    void CalculateCoverage(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex &index,
                           const std::experimental::filesystem::path &dir, const io::Library &lib);

    std::experimental::filesystem::path alignLib(logging::Logger &logger, size_t threads, KmerIndex &index,
                                                 const io::Library &align_lib,
                                                 const std::experimental::filesystem::path &dir);
}
