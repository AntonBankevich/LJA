#pragma once
#include "dbg/graph_modification.hpp"
#include "error_correction.hpp"
#include "multiplicity_estimation.hpp"
#include "sequences/edit_distance.hpp"
#include "tip_correction.hpp"
#include "correction_utils.hpp"

namespace dbg {
    size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates, bool dump = false);

    std::vector<ag::GraphPath>
    FilterAlternatives(const ag::GraphPath &initial, const std::vector<ag::GraphPath> &als,
                       size_t max_diff, double threshold);

    ag::GraphPath
    chooseBulgeCandidate(const ag::GraphPath &bulge, const dbg::DBGAlignedReadStorage &reads_storage, double threshold,
                         std::vector<ag::GraphPath> &read_alternatives, std::string &message);

    std::pair<ag::GraphPath, size_t> BestAlignmentPrefix(const ag::GraphPath &al, const Sequence &seq, size_t max_diff);

    ag::GraphPath processTip(const ag::GraphPath &tip,
                              const std::vector<ag::GraphPath> &alternatives,
                              double threshold, std::string &message);

    void initialCorrect(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                        const std::experimental::filesystem::path &out_file,
                        dbg::DBGAlignedReadStorage &reads_storage,
                        DBGAlignedReadStorage &ref_storage,
                        double threshold, double bulge_threshold, double reliable_coverage, bool diploid,
                        size_t unique_threshold, bool dump);

    class TournamentPathCorrector : public AbstractCorrectionAlgorithm {
    private:
        dbg::SparseDBG &sdbg;
        DBGAlignedReadStorage &reads_storage;
        double threshold;
        double reliable_threshold;
        bool diploid;
        size_t unique_threshold;
        size_t max_size;

        bool checkTipSize(const ag::GraphPath &tip);

    public:
        TournamentPathCorrector(dbg::SparseDBG &sdbg, DBGAlignedReadStorage &reads_storage,
                                double threshold, double reliable_threshold, bool diploid,
                                size_t unique_threshold = 60000);

        void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, DBGAlignedReadStorage &reads) override;

        std::string correctRead(const std::string &name, ag::GraphPath &path) override;
    };

    class PrimitiveBulgeCorrector : public AbstractCorrectionAlgorithm {
    private:
        double threshold;
    public:
        PrimitiveBulgeCorrector(double threshold);

        std::string correctRead(const std::string &name, ag::GraphPath &path) override;
    };

}