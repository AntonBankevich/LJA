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
#include <dbg/aln_reads_reader.hpp>
#include "dbg/path_dumping.hpp"

namespace dbg {
    std::unordered_map<std::string, std::experimental::filesystem::path>
    CoverageEC(logging::Logger &logger, const std::experimental::filesystem::path &dir,
               const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib, const io::Library &references_lib,
               size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
               bool diploid, bool dump_reads, bool debug, bool load) {
        logger.info() << "Performing coverage-based error correction with k = " << k << std::endl;
        if (k % 2 == 0) {
            logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
            k += 1;
        }
        ensure_dir_existance(dir);
        hashing::RollingHash hasher(k);
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        dbg::SparseDBG dbg = load ? LoadDBGFromEdgeSequences(logger, threads, {dir/"initial_dbg.gfa"}, hasher) :
                             DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
        ag::Printer gfa_printer;
        ag::Printer dot_printer(ag::VertexPrintStyles::defaultDotInfo(), ag::EdgePrintStyles::defaultDotInfo());
        gfa_printer.setEdgeInfo(ag::EdgeInfo({&ag::GetEdgeNameForSaving},{}, {}));
        if(debug && !load) {
            logger.info() << "Printing graph to " << (dir/"initial_dbg.gfa") << std::endl;
            gfa_printer.printGFA(dir/"initial_dbg.gfa", dbg);
            logger.info() << "Finished printing graph" << std::endl;
        }
        size_t extension_size = 800;
        std::experimental::filesystem::path al_file = dir / "initial_alignments.aln";
        dbg::SeqReader reader(reads_lib, logger, threads);
        dbg::DBGAlignedReadStorage readStorage = load ?
                dbg::DBGAlignedReadStorage::Load(logger, threads, al_file, dbg, true) :
                dbg::DBGAlignedReadStorage(logger, threads, dbg,
                    AlignReads(logger, threads, reader.begin(), reader.end(), dbg, w),
                                               true);
        if (debug && !load) {
            logger.info() << "Printing read alignments to " << al_file << std::endl;
            readStorage.Save(al_file);
            logger.info() << "Finished printing read alignments to " << al_file << std::endl;
        }
        dot_printer.printDot(dir / "initial_dbg.dot", dbg);
        if(debug) {
            readStorage.logReads(threads, dir/"read_log.txt");
            readStorage.logGraph(dbg, logger.getLoggerStream(logging::LogLevel::trace));
        }
        dbg::DBGAlignedReadStorage refStorage(logger, threads, dbg, std::vector<ag::AlignedRead>(), false);
        coverageStats(logger, dbg);
        if (debug) {
            PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, references_lib, true);
        }
        Precorrector precorrector_early(4);
        Precorrector precorrector_late(1.01);
        DimerCorrector dimerCorrector(logger, dbg, readStorage, StringContig::max_dimer_size);
        TournamentPathCorrector tournamentPathCorrector(dbg, readStorage, threshold, reliable_coverage, diploid, 60000);
        BulgePathCorrector bpCorrector(dbg, readStorage, 80000, 1);
        ErrorCorrectionEngine(precorrector_early).run(logger, threads, dbg, readStorage);
        RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
        readStorage.stopTrackSuffixes();
        dbg.resetEdgeCodes(logger, threads);
        readStorage.trackSuffixes(logger, threads, dbg, 0, extension_size);
        if(debug) readStorage.getReads().checkConsistency();
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, readStorage);
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
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, readStorage);
        ManyKCorrect(logger, threads, dbg, readStorage, threshold, reliable_coverage, 3500, 3, diploid);
        ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, readStorage);
        if (diploid)
            ErrorCorrectionEngine(bpCorrector).run(logger, threads, dbg, readStorage);
        ErrorCorrectionEngine(precorrector_late).run(logger, threads, dbg, readStorage);
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
                {"k-mer-size=501", "window=2000", "coverage-threshold=3", "reliable-coverage=10", "diploid", "load", "dump-reads"},
                {}, ""), {"reads", "pseudo_reads", "paths", "references"}, {"corrected_reads", "pseudo_reads", "final_dbg"}) {
        }

    protected:
        std::unordered_map<std::string, std::experimental::filesystem::path>
        innerRun(logging::Logger &logger, size_t threads,
                 const std::experimental::filesystem::path &dir, bool debug,
                 const AlgorithmParameterValues &parameterValues,
                 const std::unordered_map<std::string, io::Library> &input) override {
            size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
            size_t w = std::stoi(parameterValues.getValue("window"));
            double reliable_coverage = std::stod(parameterValues.getValue("reliable-coverage"));
            double threshold = std::stod(parameterValues.getValue("coverage-threshold"));
            bool diploid = parameterValues.getCheck("diploid");
            bool load = parameterValues.getCheck("load");
            return CoverageEC(logger, dir, input.find("reads")->second, input.find("pseudo_reads")->second,
                              input.find("paths")->second, input.find("references")->second, threads, k, w, threshold, reliable_coverage, diploid,
                              parameterValues.getCheck("dump-reads"), debug, load);
        }
    };
}