#pragma once

#include "dbg/dbg_read_alignment_storage.hpp"
#include "uniqueness.hpp"
#include "ff.hpp"
#include "dbg/sparse_dbg.hpp"
#include <utility>

class MappedNetwork : public Network {
public:
    std::vector<ag::EdgeId> edge_mapping;
    std::unordered_map<ag::VertexId, int> vertex_mapping;

    MappedNetwork(const ag::Component &component, const std::function<bool(const dbg::Edge &)> &unique,
                  double rel_coverage = 1000, double unique_coverage = 0, double double_coverage = 0);
    size_t addTipSinks();
    std::vector<dbg::EdgeId> getUnique(logging::Logger &logger);
    std::unordered_map<dbg::EdgeId, std::pair<size_t, size_t>> findBounds();
};

std::pair<double, double> minmaxCov(const ag::Component &subcomponent, const dbg::DBGAlignedReadStorage &reads_storage,
                                    const std::function<bool(const dbg::Edge &)> &is_unique);

class UniqueClassificator : public MultiplicityBounds {
private:
    dbg::SparseDBG &dbg;
    bool diploid;
    bool debug;
    double initial_rel_coverage;

public:
    const dbg::DBGAlignedReadStorage &reads_storage;

    void markPseudoHets() const;

    void classify(logging::Logger &logger, size_t unique_len, const std::experimental::filesystem::path &dir);
    explicit UniqueClassificator(dbg::SparseDBG &dbg, const dbg::DBGAlignedReadStorage &reads_storage, double initial_rel_coverage, bool diploid, bool debug) :
                    dbg(dbg), reads_storage(reads_storage), initial_rel_coverage(initial_rel_coverage), diploid(diploid), debug(debug) {}
    size_t ProcessUsingCoverage(logging::Logger &logger, const ag::Component &subcomponent,
                              const std::function<bool(const dbg::Edge &)> &is_unique, double rel_coverage);
    void processSimpleComponent(logging::Logger &logger, const ag::Component &component) const;
    bool processSimpleRepeat(const ag::Component &component);
    size_t processComponent(logging::Logger &logger, const ag::Component &component);
};

ag::AlignedReadStorage ResolveLoops(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads_storage,
                           const AbstractUniquenessStorage &more_unique);

