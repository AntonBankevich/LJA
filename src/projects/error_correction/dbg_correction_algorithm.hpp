#pragma once

#include "error_correction.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"

namespace dbg {
//    Downcasts ag::AssemblyGraph/ag::AlignedReadStorage to their DBG-specific counterparts exactly
//    once here, so DBG-only correction algorithms can override the typed initialize() below
//    without repeating the cast in every subclass.
    class AbstractDBGCorrectionAlgorithm : public ag::AbstractCorrectionAlgorithm {
    public:
        using ag::AbstractCorrectionAlgorithm::AbstractCorrectionAlgorithm;

        void initialize(logging::Logger &logger, size_t threads, ag::AssemblyGraph &graph, ag::AlignedReadStorage &reads) final {
            initialize(logger, threads, dynamic_cast<dbg::SparseDBG &>(graph), dynamic_cast<dbg::DBGAlignedReadStorage &>(reads));
        }

        virtual void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads) {}
    };
}
