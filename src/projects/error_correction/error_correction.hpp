#pragma once

#include "dbg/sparse_dbg.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"

class AbstractCorrectionAlgorithm {
private:
    std::string name;
public:
    virtual ~AbstractCorrectionAlgorithm() = default;

    explicit AbstractCorrectionAlgorithm(const std::string &name) : name(name) {};
    std::string getName() const {return name;}
    virtual void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads) {};
    virtual std::string correctRead(const std::string &name, ag::GraphPath &) = 0;
};

class ErrorCorrectionEngine {
private:
    AbstractCorrectionAlgorithm &algorithm;
public:
    explicit ErrorCorrectionEngine(AbstractCorrectionAlgorithm &algorithm) : algorithm(algorithm) {}

    size_t run(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &storage) {
        ag::AlignedReadStorage &reads_storage = storage.getReads();
        algorithm.initialize(logger, threads, dbg, storage);
        logger.info() << "Correcting reads using algorithm " << algorithm.getName() << std::endl;
        ParallelCounter cnt(threads);
        omp_set_num_threads(threads);
        logging::ProgressBar progressBar(logger, reads_storage.size(), threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(std::cout, reads_storage, logger, cnt, progressBar)
        for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
            ag::AlignedRead &alignedRead = reads_storage[read_ind];
            progressBar.tick();
            if (!alignedRead.valid())
                continue;
            ag::GraphPath corrected = alignedRead.getPath();
            std::string message = algorithm.correctRead(alignedRead.getId(), corrected);
            if(!message.empty()) {
                if (corrected.truncLen() >= 500) {
                    VERIFY_MSG(alignedRead.getPath() != corrected, message);
                    reads_storage.rerouteRead(alignedRead, corrected, itos(omp_get_thread_num()) + "_" + algorithm.getName() + "_" + message);
                    cnt += 1;
                } else {
#pragma omp critical
                    logger.trace() << "Removed read " << alignedRead.getId() << " due to its size: " << alignedRead.getPath().truncLen() <<"->" << corrected.truncLen() << std::endl;
                    reads_storage.delayedInvalidateRead(alignedRead, itos(omp_get_thread_num()) + "_" + algorithm.getName() + "_" + message + "_" + "invalidated_as_short");
                }
            }
        }
        progressBar.finish();
        reads_storage.applyCorrections(logger, threads);
        for(dbg::Edge &edge: dbg.edges()) edge.is_reliable = false;
        logger.info() << "Corrected " << cnt.get() << " reads using algorithm " << algorithm.getName() << std::endl;
        return cnt.get();
    }
};
