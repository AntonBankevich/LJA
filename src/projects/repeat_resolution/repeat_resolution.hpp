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
  std::vector<RecordStorage *> extra_storages;
  std::uint64_t start_k{1};
  std::uint64_t saturating_k{1};
  std::experimental::filesystem::path dir;
  uint64_t unique_threshold{0};
  bool diploid{false};
  bool debug{false};

  [[nodiscard]] std::vector<RecordStorage *> get_storages() const {
    std::vector<RecordStorage *> storages = extra_storages;
    storages.push_back(reads_storage);
    return storages;
  }

public:
  RepeatResolver(dbg::SparseDBG &dbg, RecordStorage *reads_storage,
                 std::vector<RecordStorage *> extra_storages, uint64_t start_k,
                 uint64_t saturating_k, std::experimental::filesystem::path dir,
                 uint64_t unique_threshold, bool diploid, bool debug)
      : dbg{dbg}, reads_storage{reads_storage},
        extra_storages{std::move(extra_storages)}, start_k{start_k},
        saturating_k{saturating_k}, dir{std::move(dir)},
        unique_threshold{unique_threshold}, diploid{diploid}, debug{debug} {
    std::experimental::filesystem::create_directory(this->dir);
  }

  std::vector<Contig> resolve_repeats(logging::Logger &logger) {
    logger.info() << "Resolving repeats" << std::endl;
    for (RecordStorage *const storage : get_storages()) {
      storage->invalidateSubreads(logger, 1);
    }
    RRPaths rr_paths = PathsBuilder::FromDBGStorages(dbg, get_storages());

    UniqueClassificator classificator(dbg, *reads_storage, diploid, debug);
    classificator.classify(logger, unique_threshold, dir / "mult_dir");
    MultiplexDBG mdbg(dbg, &rr_paths, start_k, classificator);
    mdbg.SerializeToDot(dir / "init_graph.dot");
    logger.info() << "Increasing k" << std::endl;
    MultiplexDBGIncreaser k_increaser{start_k, saturating_k, logger, debug};
    k_increaser.IncreaseUntilSaturation(mdbg, true);
    logger.info() << "Finished increasing k" << std::endl;
    mdbg.SerializeToDot(dir / "resolved_graph.dot");
    mdbg.SerializeToGFA(dir / "resolved_graph.gfa");

    std::vector<Contig> edges =
        mdbg.PrintTrimEdges(dir / "compressed.fasta");
    return edges;
  }
};

} // End namespace repeat_resolution