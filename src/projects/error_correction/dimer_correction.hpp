#pragma once
#include "error_correction.hpp"
#include "dbg/graph_modification.hpp"
#include "common/logging.hpp"
using namespace dbg;
class DimerCorrector : public AbstractCorrectionAlgorithm {
private:
    dbg::SparseDBG &sdbg;
    dbg::DBGAlignedReadStorage &reads_storage;
    logging::Logger &logger;
    size_t max_at;
public:
    DimerCorrector(logging::Logger &logger, dbg::SparseDBG &sdbg, dbg::DBGAlignedReadStorage &reads_storage, size_t max_at) :
            AbstractCorrectionAlgorithm("DimerCorrector"), logger(logger), sdbg(sdbg), reads_storage(reads_storage), max_at(max_at) {}

    std::string correctRead(const std::string &name, ag::GraphPath &path) override;
};
