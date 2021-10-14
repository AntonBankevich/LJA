#pragma once

size_t Precorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads_storage,
                    double reliable_threshold);
