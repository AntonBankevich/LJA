#pragma once
#include "dbg/graph_modification.hpp"
#include "dbg/compact_path.hpp"
#include "multiplicity_estimation.hpp"
#include "sequences/edit_distance.hpp"
#include "tip_correction.hpp"
#include "correction_utils.hpp"

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump = false);
std::vector<dbg::Path> FindBulgeAlternatives(const dbg::Path &path, size_t max_diff);
std::vector<dbg::GraphAlignment> FilterAlternatives(const dbg::GraphAlignment &initial, const std::vector<dbg::GraphAlignment> &als,
                                            size_t max_diff, double threshold);
dbg::GraphAlignment chooseBulgeCandidate(const dbg::GraphAlignment &bulge, const RecordStorage &reads_storage, double threshold,
                                    std::vector<dbg::GraphAlignment> &read_alternatives, std::string &message);
std::pair<dbg::GraphAlignment, size_t> BestAlignmentPrefix(const dbg::GraphAlignment &al, const Sequence & seq);
dbg::GraphAlignment processTip(const dbg::GraphAlignment &tip,
                         const std::vector<dbg::GraphAlignment> & alternatives,
                         double threshold, std::string &message);
size_t correctLowCoveredRegions(logging::Logger &logger, dbg::SparseDBG &sdbg,RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                const std::experimental::filesystem::path &out_file,
                                double threshold, double reliable_threshold, bool diploid, size_t threads, bool dump);
dbg::GraphAlignment findAlternative(logging::Logger &logger, std::ostream &out, const dbg::GraphAlignment &bulge,
                                   const RecordStorage &reads_storage);
size_t collapseBulges(logging::Logger &logger, RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                const std::experimental::filesystem::path &out_file,
                                double threshold, size_t k, size_t threads);
size_t correctAT(logging::Logger &logger, size_t threads, RecordStorage &reads_storage, size_t max_at);
void initialCorrect(dbg::SparseDBG &sdbg, logging::Logger &logger,
                    const std::experimental::filesystem::path &out_file,
                    RecordStorage &reads_storage,
                    RecordStorage &ref_storage,
                    double threshold, double bulge_threshold, double reliable_coverage, bool diploid, size_t threads, bool dump);