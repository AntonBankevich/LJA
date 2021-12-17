#pragma once

#include "sequences/edit_distance.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/graph_alignment_storage.hpp"
void CorrectTips(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                 const std::vector<RecordStorage *> &storages);

void TipCorrectionPipeline(logging::Logger &logger, dbg::SparseDBG &dbg, RecordStorage &reads,
                           size_t threads,double reliable_threshold);
