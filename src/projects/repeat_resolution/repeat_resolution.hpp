//
// Created by Andrey Bzikadze on 10/13/21.
//

#pragma once

#include "../dbg/graph_alignment_storage.hpp"
#include "../dbg/sparse_dbg.hpp"
#include "../error_correction/multiplicity_estimation.hpp"
#include "multiplex_dbg.hpp"
#include "paths.hpp"
#include <graphlite/serialize.hpp>
#include <vector>

namespace repeat_resolution {

class RepeatResolver {
  dbg::SparseDBG &dbg;
  RecordStorage *reads_storage;
  std::vector<RecordStorage *> extra_storages;
  std::uint64_t start_k{1};
  std::experimental::filesystem::path dir;
  uint64_t unique_threshold;
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
                 std::experimental::filesystem::path dir,
                 uint64_t unique_threshold, bool diploid, bool debug)
      : dbg{dbg}, reads_storage{reads_storage}, extra_storages{std::move(
                                                    extra_storages)},
        start_k{start_k}, dir{std::move(dir)},
        unique_threshold{unique_threshold}, diploid{diploid}, debug{debug} {
    std::experimental::filesystem::create_directory(this->dir);
  }

  void resolve_repeats(logging::Logger &logger) {
    logger.info() << "Resolving repeats" << std::endl;
    for (RecordStorage *const storage : get_storages()) {
      storage->invalidateSubreads(logger, 1);
    }
    RRPaths rr_paths = PathsBuilder::FromDBGStorages(dbg, get_storages());

    UniqueClassificator classificator(dbg, *reads_storage, diploid, debug);
    classificator.classify(logger, unique_threshold, dir/"mult_dir");
    MultiplexDBG mdbg(dbg, &rr_paths, start_k, classificator, debug, dir,
                      logger);
    mdbg.serialize_to_dot(dir / "init_graph.dot");
    logger.info() << "Increasing k" << std::endl;
    mdbg.incN(100000, debug);
    logger.info() << "Finished increasing k" << std::endl;
    mdbg.serialize_to_dot(dir / "resolved_graph.dot");
  }
};

} // End namespace repeat_resolution