#pragma once
#include "dbg/graph_modification.hpp"
#include "dbg/compact_path.hpp"
#include "multiplicity_estimation.hpp"
#include "edit_distance.hpp"
#include "tip_correction.hpp"
#include "correction_utils.hpp"

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump = false);
std::vector<Path> FindBulgeAlternatives(const Path &path, size_t max_diff);
std::vector<GraphAlignment> FilterAlternatives(logging::Logger &logger1, const GraphAlignment &initial, const std::vector<GraphAlignment> &als,
                                            size_t max_diff, double threshold);
GraphAlignment chooseBulgeCandidate(logging::Logger &logger, std::ostream &out, const GraphAlignment &bulge,
                                    const RecordStorage &reads_storage, const RecordStorage &ref_storage, double threshold,
                                    std::vector<GraphAlignment> &read_alternatives, std::string &message, bool dump);
std::pair<GraphAlignment, size_t> BestAlignmentPrefix(const GraphAlignment &al, const Sequence & seq);
GraphAlignment processTip(logging::Logger &logger, std::ostream &out, const GraphAlignment &tip,
                         const std::vector<GraphAlignment> & alternatives,
                         const RecordStorage &ref_storage,
                         double threshold, std::string &message, bool dump = false);
size_t correctLowCoveredRegions(logging::Logger &logger, SparseDBG &sdbg,RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                const std::experimental::filesystem::path &out_file,
                                double threshold, double reliable_threshold, size_t k, size_t threads, bool dump);
GraphAlignment findAlternative(logging::Logger &logger, std::ostream &out, const GraphAlignment &bulge,
                                   const RecordStorage &reads_storage);
size_t collapseBulges(logging::Logger &logger, RecordStorage &reads_storage,
                                RecordStorage &ref_storage,
                                const std::experimental::filesystem::path &out_file,
                                double threshold, size_t k, size_t threads);
size_t correctAT(logging::Logger &logger, RecordStorage &reads_storage, size_t k, size_t threads);
void initialCorrect(SparseDBG &sdbg, logging::Logger &logger,
                    const std::experimental::filesystem::path &out_file,
                    RecordStorage &reads_storage,
                    RecordStorage &ref_storage,
                    double threshold, double bulge_threshold, double reliable_coverage, size_t threads, bool dump);