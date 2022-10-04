//
// Created by Andrey Bzikadze on 10/13/21.
//

#pragma once

#include "../dbg/graph_alignment_storage.hpp"
#include "../dbg/sparse_dbg.hpp"
#include "../error_correction/multiplicity_estimation.hpp"
#include "mdbg.hpp"
#include "mdbg_inc.hpp"
#include "paths.hpp"
#include <graphlite/serialize.hpp>
#include <vector>

namespace repeat_resolution {

class RepeatResolver {
    dbg::SparseDBG &dbg;
    RecordStorage *reads_storage;
    std::vector<RecordStorage *> extra_storages{};
    std::uint64_t start_k{1};
    std::uint64_t saturating_k{1};
    std::experimental::filesystem::path dir;
    uint64_t unique_threshold{0};
    bool diploid{false};
    bool debug{false};
    UniqueClassificator classificator;

    [[nodiscard]] std::vector<RecordStorage *> get_storages() const {
        std::vector<RecordStorage *> storages = extra_storages;
        storages.push_back(reads_storage);
        return storages;
    }

 public:
    RepeatResolver(dbg::SparseDBG &dbg,
                   RecordStorage *reads_storage,
                   std::vector<RecordStorage *> extra_storages,
                   uint64_t start_k,
                   uint64_t saturating_k,
                   const std::experimental::filesystem::path &dir,
                   uint64_t unique_threshold,
                   bool diploid,
                   bool debug,
                   logging::Logger &logger)
        : dbg{dbg}, reads_storage{std::move(reads_storage)},
          extra_storages{std::move(extra_storages)}, start_k{start_k},
          saturating_k{saturating_k}, dir{std::move(dir)},
          unique_threshold{unique_threshold}, diploid{diploid}, debug{debug},
          classificator{dbg, *(this->reads_storage), 0, diploid, debug} {
        std::experimental::filesystem::create_directory(this->dir);
        classificator.classify(logger, unique_threshold, dir/"mult_dir");
        // TODO reactivate filtering
//        for (RecordStorage *const storage : get_storages()) {
//            storage->invalidateSubreads(logger, 1);
//        }
    }

    void ResolveRepeats(logging::Logger &logger, size_t threads) {
        logger.info() << "Resolving repeats" << std::endl;
        logger.info() << "Constructing paths" << std::endl;
        RRPaths rr_paths = PathsBuilder::FromDBGStorages(dbg, get_storages());

        logger.info() << "Building graph" << std::endl;
        MultiplexDBG mdbg(dbg, &rr_paths, start_k, classificator);
        if (debug) {
            logger.trace() << "Checking validity of graph" << std::endl;
            mdbg.AssertValidity();
            logger.trace() << "Graph has passed validity check" << std::endl;
        }
        // logger.info() << "Export to dot" << std::endl;
        // mdbg.ExportToDot(dir/"init_graph.dot");
        // logger.info() << "Export to GFA" << std::endl;
        // mdbg.ExportToGFA(dir/"init_graph.gfa");

        logger.info() << "Increasing k" << std::endl;
        MultiplexDBGIncreaser k_increaser{start_k, saturating_k, logger, debug};
        k_increaser.IncreaseUntilSaturation(mdbg, true);
        logger.info() << "Finished increasing k" << std::endl;

        logger.info() << "Exporting remaining active transitions" << std::endl;
        mdbg.ExportActiveTransitions(dir/"mdbg_remaining_trans.txt");

        logger.info() << "Export to Dot" << std::endl;
        mdbg.ExportToDot(dir/"mdbg.hpc.dot");
        logger.info() << "Export of unique edges" << std::endl;
        mdbg.ExportUniqueEdges(dir/"mdbg_is_edge_unique.tsv");
        logger.info() << "Export to GFA and compressed contigs" << std::endl;
        std::vector<Contig> contigs = mdbg.ExportContigsAndGFA(
            dir/"assembly.hpc.fasta", dir/"mdbg.hpc.gfa", threads, logger);
        logger.info() << "Finished repeat resolution" << std::endl;
    }
};

} // End namespace repeat_resolution