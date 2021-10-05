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
    std::unordered_map<int, dbg::Edge *> edge_mapping;
    std::unordered_map<dbg::Vertex *, int> vertex_mapping;

    MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique,
                  double rel_coverage = 1000, double unique_coverage = 0, double double_coverage = 0);
    size_t addTipSinks();
    std::vector<dbg::Edge*> getUnique(logging::Logger &logger);
};

class MultiplicityBoundsEstimator {
private:
    SparseDBG &dbg;
    MultiplicityBounds bounds;
public:
    MultiplicityBoundsEstimator(SparseDBG &dbg, const AbstractUniquenessStorage &uniquenessStorage);

    bool updateComponent(logging::Logger &logger, const Component &component, const AbstractUniquenessStorage &uniquenessStorage,
                                double rel_coverage, double unique_coverage = 0);
    void update(logging::Logger &logger, double rel_coverage, const std::experimental::filesystem::path &dir);
};
std::pair<double, double> minmaxCov(const Component &subcomponent, const RecordStorage &reads_storage,
                                    const std::function<bool(const dbg::Edge &)> &is_unique);

class UniqueClassificator : public SetUniquenessStorage{
private:
    SparseDBG &dbg;
    bool diploid;
    bool debug;

public:
    const RecordStorage &reads_storage;

    void markPseudoHets() const;

    void classify(logging::Logger &logger, size_t unique_len, const std::experimental::filesystem::path &dir);
    explicit UniqueClassificator(SparseDBG &dbg, const RecordStorage &reads_storage, bool diploid, bool debug) :
                    dbg(dbg), reads_storage(reads_storage), diploid(diploid), debug(debug) {}
    std::vector<dbg::Edge *> ProcessUsingCoverage(logging::Logger &logger, const Component &subcomponent,
                              const std::function<bool(const dbg::Edge &)> &is_unique, double rel_coverage) const;
    void processSimpleComponent(logging::Logger &logger, const Component &component) const;
    std::vector<const dbg::Edge *> processComponent(logging::Logger &logger, const Component &component) const;
};

RecordStorage ResolveLoops(logging::Logger &logger, size_t threads, SparseDBG &dbg, RecordStorage &reads_storage,
                           const AbstractUniquenessStorage &more_unique);

