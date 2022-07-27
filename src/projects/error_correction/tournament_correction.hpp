#pragma once
#include "dbg/graph_modification.hpp"
#include "dbg/compact_path.hpp"
#include "error_correction.hpp"
#include "multiplicity_estimation.hpp"
#include "sequences/edit_distance.hpp"
#include "tip_correction.hpp"
#include "correction_utils.hpp"

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump = false);
std::vector<dbg::GraphAlignment> FilterAlternatives(const dbg::GraphAlignment &initial, const std::vector<dbg::GraphAlignment> &als,
                                            size_t max_diff, double threshold);
dbg::GraphAlignment chooseBulgeCandidate(const dbg::GraphAlignment &bulge, const RecordStorage &reads_storage, double threshold,
                                    std::vector<dbg::GraphAlignment> &read_alternatives, std::string &message);
std::pair<dbg::GraphAlignment, size_t> BestAlignmentPrefix(const dbg::GraphAlignment &al, const Sequence & seq);
dbg::GraphAlignment processTip(const dbg::GraphAlignment &tip,
                         const std::vector<dbg::GraphAlignment> & alternatives,
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
    std::string correctRead(dbg::GraphAlignment &path) override;
};

class PrimitiveBulgeCorrector : public AbstractCorrectionAlgorithm {
private:
    double threshold;
public:
    PrimitiveBulgeCorrector(double threshold);
    std::string correctRead(dbg::GraphAlignment &path) override;
};

