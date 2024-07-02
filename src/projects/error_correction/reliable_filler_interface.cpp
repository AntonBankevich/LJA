#include "reliable_filler_interface.hpp"

using namespace dbg;
size_t AbstractReliableFillingAlgorithm::ReFill(SparseDBG &dbg) {
    for (Edge &edge: dbg.edges()) {
        edge.is_reliable = false;
        edge.mark(ag::common);
    }
    return Fill(dbg);
}

size_t AbstractReliableFillingAlgorithm::LoggedReFill(logging::Logger &logger, SparseDBG &dbg) {
    for (Edge &edge: dbg.edges()) {
        edge.is_reliable = false;
        edge.mark(ag::common);
    }
    return LoggedFill(logger, dbg);
}

size_t AbstractReliableFillingAlgorithm::LoggedFill(logging::Logger &logger, SparseDBG &dbg) {
    logger.info() << "Running reliable marker " << name() << std::endl;
    size_t res = Fill(dbg);
    logger.info() << "Marked " << res << " edges as reliable" << std::endl;
    return res;
}

CompositeReliableFiller::CompositeReliableFiller(std::vector<AbstractReliableFillingAlgorithm *> &&algorithms)
        : algorithms(algorithms) {
    std::vector<std::string> names;
    for (AbstractReliableFillingAlgorithm *alg: algorithms) {
        names.emplace_back(alg->name());
    }
    _name = join(" ", names);
}

size_t CompositeReliableFiller::Fill(SparseDBG &dbg) {
    size_t res = 0;
    for (AbstractReliableFillingAlgorithm *alg: algorithms) {
        res += alg->Fill(dbg);
    }
    return res;
}
