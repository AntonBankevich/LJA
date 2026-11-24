#pragma once
#include "dbg_graph_aligner.hpp"
#include "assembly_graph/random_access_paths.hpp"
#include "sparse_dbg.hpp"
#include "sequences/read_source.hpp"
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

//        SparseDBG Subgraph(std::vector<Segment<Edge>> &pieces) const;

        void addAllKmers(SparseDBG &dbg, const std::vector<Sequence> &new_seqs, KmerIndex &index) const;

        void checkIndexConsistency(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex& index) const;
    };

//    Fills the edges of sdbg from the sequences of at least min_size nucleotides.
    void FillSparseDBGEdges(SparseDBG &sdbg, const std::vector<Sequence> &seqs, logging::Logger &logger,
                            size_t threads, size_t min_size);
//    Fills the edges of sdbg from the reads of lib, cut into pieces by splitter, of at least min_size nucleotides.
    void FillSparseDBGEdges(SparseDBG &sdbg, const io::Library &lib, logging::Logger &logger, size_t threads,
                            size_t min_size, const io::Splitter &splitter);

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

    SparseDBG constructSparseDBGFromReads(logging::Logger &logger, const io::Library &reads_file, size_t threads,
                                          const hashing::RollingHash &hasher,
                                          const std::vector<hashing::htype> &hash_list, size_t w);

    void CalculateCoverage(logging::Logger &logger, size_t threads, SparseDBG &dbg, KmerIndex &index,
                           const std::experimental::filesystem::path &dir, const io::Library &lib);

    std::experimental::filesystem::path alignLib(logging::Logger &logger, size_t threads, KmerIndex &index,
                                                 const io::Library &align_lib,
                                                 const std::experimental::filesystem::path &dir);
}
