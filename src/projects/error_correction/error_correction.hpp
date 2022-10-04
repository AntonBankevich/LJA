#pragma once

#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_alignment_storage.hpp"

class AbstractCorrectionAlgorithm {
private:
    std::string name;
public:
    AbstractCorrectionAlgorithm(const std::string &name) : name(name) {};
    std::string getName() const {return name;}
    virtual void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads) {};
    virtual std::string correctRead(dbg::GraphAlignment &) = 0;
};

class ErrorCorrectionEngine {
private:
    AbstractCorrectionAlgorithm &algorithm;
public:
    explicit ErrorCorrectionEngine(AbstractCorrectionAlgorithm &algorithm) : algorithm(algorithm) {}

    size_t run(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads_storage) {
        algorithm.initialize(logger, threads, dbg, reads_storage);
        logger.info() << "Correcting reads using algorithm " << algorithm.getName() << std::endl;
        ParallelCounter cnt(threads);
        omp_set_num_threads(threads);
        logging::ProgressBar progressBar(logger, reads_storage.size(), threads);
#pragma omp parallel for default(none) schedule(dynamic, 100) shared(std::cout, reads_storage, logger, cnt, progressBar)
        for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
            AlignedRead &alignedRead = reads_storage[read_ind];
            progressBar.tick();
            if (!alignedRead.valid())
                continue;
            dbg::CompactPath &initial_cpath = alignedRead.path;
            dbg::GraphAlignment corrected = initial_cpath.getAlignment();
            std::string message = algorithm.correctRead(corrected);
            if(!message.empty()) {
                reads_storage.reroute(alignedRead, corrected, itos(omp_get_thread_num()) + "_" + algorithm.getName() + "_" + message);
                cnt += 1;
            }
        }
        progressBar.finish();
        reads_storage.applyCorrections(logger, threads);
        logger.info() << "Corrected " << cnt.get() << " reads using algorithm " << algorithm.getName() << std::endl;
        return cnt.get();
    }
};
