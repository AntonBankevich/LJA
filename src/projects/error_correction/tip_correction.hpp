#pragma once

#include "sequences/edit_distance.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
size_t CorrectTips(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                 std::vector<dbg::DBGAlignedReadStorage *> storages);

void TipCorrectionPipeline(logging::Logger &logger, dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads,
                           size_t threads,double reliable_threshold);
