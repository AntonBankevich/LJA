#pragma once
#include "dbg/graph_modification.hpp"
#include "dbg/compact_path.hpp"
#include "error_correction.hpp"
#include "multiplicity_estimation.hpp"
#include "sequences/edit_distance.hpp"
#include "tip_correction.hpp"
#include "correction_utils.hpp"

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump = false);
std::vector<DBGGraphPath> FilterAlternatives(const DBGGraphPath &initial, const std::vector<DBGGraphPath> &als,
                                             size_t max_diff, double threshold);
DBGGraphPath chooseBulgeCandidate(const DBGGraphPath &bulge, const RecordStorage &reads_storage, double threshold,
                                  std::vector<DBGGraphPath> &read_alternatives, std::string &message);
std::pair<DBGGraphPath, size_t> BestAlignmentPrefix(const DBGGraphPath &al, const Sequence & seq);
DBGGraphPath processTip(const DBGGraphPath &tip,
                        const std::vector<DBGGraphPath> & alternatives,
                        double threshold, std::string &message);
size_t collapseBulges(logging::Logger &logger, RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                double threshold, size_t k, size_t threads);
void initialCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                    const std::experimental::filesystem::path &out_file,
                    RecordStorage &reads_storage,
                    RecordStorage &ref_storage,
                    double threshold, double bulge_threshold, double reliable_coverage, bool diploid, size_t unique_threshold, bool dump);

class TournamentPathCorrector : public AbstractCorrectionAlgorithm {
private:
    dbg::SparseDBG &sdbg;
    RecordStorage &reads_storage;
    double threshold;
    double reliable_threshold;
    bool diploid;
    size_t unique_threshold;
    size_t max_size;
public:
    TournamentPathCorrector(dbg::SparseDBG &sdbg, RecordStorage &reads_storage,
                            double threshold, double reliable_threshold, bool diploid, size_t unique_threshold = 60000);
    void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads) override;
    std::string correctRead(DBGGraphPath &path) override;
};

class PrimitiveBulgeCorrector : public AbstractCorrectionAlgorithm {
private:
    double threshold;
public:
    PrimitiveBulgeCorrector(double threshold);
    std::string correctRead(DBGGraphPath &path) override;
};

