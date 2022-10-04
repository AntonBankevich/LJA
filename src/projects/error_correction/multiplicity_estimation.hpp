#pragma once

#include "dbg/graph_alignment_storage.hpp"
#include "uniqueness.hpp"
#include "ff.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/visualization.hpp"
#include "dbg/compact_path.hpp"
#include <utility>

class MappedNetwork : public Network {
public:
    std::vector<dbg::Edge *> edge_mapping;
    std::unordered_map<dbg::Vertex *, int> vertex_mapping;

    MappedNetwork(const dbg::Component &component, const std::function<bool(const dbg::Edge &)> &unique,
                  double rel_coverage = 1000, double unique_coverage = 0, double double_coverage = 0);
    size_t addTipSinks();
    std::vector<dbg::Edge*> getUnique(logging::Logger &logger);
    std::unordered_map<dbg::Edge *, std::pair<size_t, size_t>> findBounds();
};

std::pair<double, double> minmaxCov(const dbg::Component &subcomponent, const RecordStorage &reads_storage,
                                    const std::function<bool(const dbg::Edge &)> &is_unique);

    class UniqueClassificator : public MultiplicityBounds {
private:
    dbg::SparseDBG &dbg;
    bool diploid;
    bool debug;
    double initial_rel_coverage;

public:
    const RecordStorage &reads_storage;

    void markPseudoHets() const;

    void classify(logging::Logger &logger, size_t unique_len, const std::experimental::filesystem::path &dir);
    explicit UniqueClassificator(dbg::SparseDBG &dbg, const RecordStorage &reads_storage, double initial_rel_coverage, bool diploid, bool debug) :
                    dbg(dbg), reads_storage(reads_storage), initial_rel_coverage(initial_rel_coverage), diploid(diploid), debug(debug) {}
    size_t ProcessUsingCoverage(logging::Logger &logger, const dbg::Component &subcomponent,
                              const std::function<bool(const dbg::Edge &)> &is_unique, double rel_coverage);
    void processSimpleComponent(logging::Logger &logger, const dbg::Component &component) const;
    bool processSimpleRepeat(const dbg::Component &component);
    size_t processComponent(logging::Logger &logger, const dbg::Component &component);
};

RecordStorage ResolveLoops(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads_storage,
                           const AbstractUniquenessStorage &more_unique);

