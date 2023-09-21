#pragma once
#include "dbg/graph_modification.hpp"
#include "dbg/compact_path.hpp"
#include "error_correction.hpp"
#include "multiplicity_estimation.hpp"
#include "sequences/edit_distance.hpp"
#include "tip_correction.hpp"
#include "correction_utils.hpp"

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump = false);
std::vector<dbg::GraphPath> FilterAlternatives(const dbg::GraphPath &initial, const std::vector<dbg::GraphPath> &als,
                                             size_t max_diff, double threshold);
dbg::GraphPath chooseBulgeCandidate(const dbg::GraphPath &bulge, const RecordStorage &reads_storage, double threshold,
                                  std::vector<dbg::GraphPath> &read_alternatives, std::string &message);
std::pair<dbg::GraphPath, size_t> BestAlignmentPrefix(const dbg::GraphPath &al, const Sequence & seq);
dbg::GraphPath processTip(const dbg::GraphPath &tip,
                        const std::vector<dbg::GraphPath> & alternatives,
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

    bool checkTipSize(const dbg::GraphPath &tip);
public:
    TournamentPathCorrector(dbg::SparseDBG &sdbg, RecordStorage &reads_storage,
                            double threshold, double reliable_threshold, bool diploid, size_t unique_threshold = 60000);
    void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads) override;
    std::string correctRead(dbg::GraphPath &path) override;
};

class PrimitiveBulgeCorrector : public AbstractCorrectionAlgorithm {
private:
    double threshold;
public:
    PrimitiveBulgeCorrector(double threshold);
    std::string correctRead(dbg::GraphPath &path) override;
};

