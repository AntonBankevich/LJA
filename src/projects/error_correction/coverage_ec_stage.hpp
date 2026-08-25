#pragma once

#include "tournament_correction.hpp"
#include "parameter_estimator.hpp"
#include "partial_rr.hpp"
#include "precorrection.hpp"
#include "dimer_correction.hpp"
#include "manyk_correction.hpp"
#include <dbg/dbg_construction.hpp>
#include <dbg/graph_printing.hpp>
#include <dbg/graph_stats.hpp>
#include <dbg/graph_algorithms.hpp>
#include "dbg/path_dumping.hpp"

namespace dbg {
    std::unordered_map<std::string, std::experimental::filesystem::path>
    CoverageEC(logging::Logger &logger, const std::experimental::filesystem::path &dir,
               const io::Library &graph_lib, const io::Library &read_alignments_lib,
               const io::Library &paths_lib, const io::Library &references_lib,
               size_t threads, size_t k, double threshold, double reliable_coverage,
               bool diploid, bool dump_reads, bool debug) {
        logger.info() << "Performing coverage-based error correction with k = " << k << std::endl;
        if (k % 2 == 0) {
            logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
            k += 1;
        }
        ensure_dir_existance(dir);
        hashing::RollingHash hasher(k);
        dbg::SparseDBG dbg = LoadDBGFromEdgeSequences(logger, threads, graph_lib, hasher);
        ag::Printer gfa_printer;
        ag::Printer dot_printer(ag::VertexPrintStyles::defaultDotInfo(), ag::EdgePrintStyles::defaultDotInfo());
        gfa_printer.setEdgeInfo(ag::EdgeInfo({&ag::GetEdgeNameForSaving},{}, {}));
        size_t extension_size = 800;
        dbg::DBGAlignedReadStorage readStorage = dbg::DBGAlignedReadStorage::Load(logger, threads, read_alignments_lib.front(), dbg, true);
        if(debug) {
            readStorage.logReads(threads, dir/"read_log.txt");
            readStorage.logGraph(dbg, logger.getLoggerStream(logging::LogLevel::trace));
        }
        dbg::DBGAlignedReadStorage refStorage(logger, threads, dbg, std::vector<ag::AlignedRead>(), false);
        coverageStats(logger, dbg);
        if (debug) {
            PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, references_lib, true);
        }
        Precorrector precorrector_early(
                [](const ag::Edge &e) { return e.getCoverage() >= 4 || e.is_reliable; },
                [](const ag::Edge &e) { return e.getCoverage() == 1; });
        Precorrector precorrector_late(
                [](const ag::Edge &e) { return e.getCoverage() >= 1.01 || e.is_reliable; },
                [](const ag::Edge &e) { return e.getCoverage() == 1; });
        DimerCorrector dimerCorrector(logger, dbg, readStorage, StringContig::max_dimer_size);
        TournamentPathCorrector tournamentPathCorrector(dbg, readStorage, threshold, reliable_coverage, diploid, 60000);
        BulgePathCorrector bpCorrector(dbg, readStorage, 80000, 1);
        ag::ErrorCorrectionEngine(precorrector_early).run(logger, threads, dbg, readStorage);
        RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
        readStorage.stopTrackSuffixes();
        dbg.resetEdgeCodes(logger, threads);
        readStorage.trackSuffixes(logger, threads, dbg, 0, extension_size);
        if(debug) readStorage.getReads().checkConsistency();
        ag::ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, readStorage);
        RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
        DatasetParameters params = EstimateDatasetParameters(dbg, readStorage, true);
        params.PrintBasic(logger.getLoggerStream(logging::LogLevel::info));
        params.PrintStatistics(logger.getLoggerStream(logging::LogLevel::trace));
        ManyKCorrect(logger, threads, dbg, readStorage, threshold, reliable_coverage, 800, 4, diploid);
        if (debug)
            PrintPaths(logger, threads, dir / "state_dump", "mk800", dbg, readStorage, paths_lib, references_lib, true);
        readStorage.stopTrackSuffixes();
        RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
        readStorage.trackSuffixes(logger, threads, dbg, 0, 3000);
        ManyKCorrect(logger, threads, dbg, readStorage, threshold, reliable_coverage, 2000, 4, diploid);
        if (debug)
            PrintPaths(logger, threads, dir / "state_dump", "mk2000", dbg, readStorage, paths_lib, references_lib, true);
        readStorage.stopTrackSuffixes();
        RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
        dbg.resetEdgeCodes(logger, threads);
        readStorage.trackSuffixes(logger, threads, dbg, 0, 1000000);
        ag::ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, readStorage);
        ManyKCorrect(logger, threads, dbg, readStorage, threshold, reliable_coverage, 3500, 3, diploid);
        ag::ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, readStorage);
        if (diploid)
            ag::ErrorCorrectionEngine(bpCorrector).run(logger, threads, dbg, readStorage);
        ag::ErrorCorrectionEngine(precorrector_late).run(logger, threads, dbg, readStorage);
        RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
        {
            std::vector<ag::AlignedRead> pseudo_reads = PartialRR(logger, threads, dbg, readStorage.getSuffixes());
            ag::AlignedReadStorage pseudo_reads_storage(dbg, std::move(pseudo_reads));
            pseudo_reads_storage.printSequences(dir / "pseudo_reads.fasta");
        }
        readStorage.stopTrackSuffixes();
        dbg.resetEdgeCodes(logger, threads);
        if(debug) readStorage.getReads().checkConsistency();
        coverageStats(logger, dbg);
        if (debug)
            PrintPaths(logger, threads, dir / "state_dump", "mk3500", dbg, readStorage, paths_lib, references_lib, false);
        if(dump_reads)
            readStorage.getReads().printReadFasta(logger, dir / "corrected_reads.fasta");
        std::experimental::filesystem::path corrected_reads = dir / "corrected_reads.paths";
        readStorage.getReads().printReadPaths(logger, dir / "corrected_reads.aln",
                                   dir / "final_dbg.gfa", corrected_reads, k);
//    dbg.printFastaOld(dir / "final_dbg.fasta"); ???
        gfa_printer.setEdgeInfo(ag::EdgeInfo({&ag::GetEdgeNameForSaving},{}, {}));
        gfa_printer.printGFA(dir / "final_dbg.gfa", ag::Component(dbg), true);
        gfa_printer.setEdgeInfo(ag::EdgeInfo({&ag::GetEdgeNameForSaving}, {}, {}));
        dot_printer.printDot(dir / "final_dbg.dot", ag::Component(dbg));
        logger.info() << "Initial correction results with k = " << k << " printed to " << corrected_reads << std::endl;
        return {{"corrected_reads", corrected_reads},
                {"pseudo_reads",    dir / "pseudo_reads.fasta"},
                {"final_dbg",       dir / "final_dbg.gfa"}};
    }

    class CoverageCorrectionStage : public Stage {
    public:
        CoverageCorrectionStage() : Stage(AlgorithmParameters(
                {"k-mer-size=501", "coverage-threshold=3", "reliable-coverage=10", "diploid", "dump-reads"},
                {}, ""), {"graph", "read_alignments", "paths", "references"}, {"corrected_reads", "pseudo_reads", "final_dbg"}) {
        }

    protected:
        std::unordered_map<std::string, std::experimental::filesystem::path>
        innerRun(logging::Logger &logger, size_t threads,
                 const std::experimental::filesystem::path &dir, bool debug,
                 const AlgorithmParameterValues &parameterValues,
                 const std::unordered_map<std::string, io::Library> &input) override {
            size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
            double reliable_coverage = std::stod(parameterValues.getValue("reliable-coverage"));
            double threshold = std::stod(parameterValues.getValue("coverage-threshold"));
            bool diploid = parameterValues.getCheck("diploid");
            return CoverageEC(logger, dir, input.find("graph")->second, input.find("read_alignments")->second,
                              input.find("paths")->second, input.find("references")->second, threads, k, threshold, reliable_coverage, diploid,
                              parameterValues.getCheck("dump-reads"), debug);
        }
    };
}