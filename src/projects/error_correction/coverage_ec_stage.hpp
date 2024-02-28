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
#include <sequences/seqio.hpp>

namespace dbg {
    std::unordered_map<std::string, std::experimental::filesystem::path>
    CoverageEC(logging::Logger &logger, const std::experimental::filesystem::path &dir,
               const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
               size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
               bool diploid, bool debug, bool load) {
        logger.info() << "Performing coverage-based error correction with k = " << k << std::endl;
        if (k % 2 == 0) {
            logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
            k += 1;
        }
        ensure_dir_existance(dir);
        hashing::RollingHash hasher(k);
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        dbg::SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, construction_lib, dir, threads,
                                                (dir / "disjointigs.fasta").string(), (dir / "vertices.save").string())
                                  :
                             DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
        KmerIndex index(dbg);
        index.fillAnchors(logger, threads, dbg, w);
        size_t extension_size = std::max<size_t>(k * 2, 1000);
        ag::ReadLogger readLogger(threads, dir / "read_log.txt");
        dbg::ReadAlignmentStorage readStorage(dbg, 0, extension_size, true, true, false);
        readStorage.setReadLogger(readLogger);
        dbg::ReadAlignmentStorage refStorage(dbg, 0, extension_size, false, false);
        refStorage.setReadLogger(readLogger);
        io::SeqReader reader(reads_lib);
        readStorage.FillAlignments(logger, threads, reader.begin(), reader.end(), dbg, index);
        printDot(dir / "initial_dbg.dot", Component(dbg), ag::SaveEdgeName<DBGTraits>);
        coverageStats(logger, dbg);
        if (debug) {
            PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
        }
        Precorrector precorrector(4);
        DimerCorrector dimerCorrector(logger, dbg, readStorage, StringContig::max_dimer_size);
        TournamentPathCorrector tournamentPathCorrector(dbg, readStorage, threshold, reliable_coverage, diploid, 60000);
        BulgePathCorrector bpCorrector(dbg, readStorage, 80000, 1);
        ErrorCorrectionEngine(precorrector).run(logger, threads, dbg, readStorage);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, extension_size);
        readStorage.trackSuffixes(logger, threads);
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, readStorage);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, extension_size);
        DatasetParameters params = EstimateDatasetParameters(dbg, readStorage, true);
        params.Print(logger);
        ManyKCorrect(logger, threads, dbg, readStorage, threshold, reliable_coverage, 800, 4, diploid);
        if (debug)
            PrintPaths(logger, threads, dir / "state_dump", "mk800", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 5 / 2, 3000));
        ManyKCorrect(logger, threads, dbg, readStorage, threshold, reliable_coverage, 2000, 4, diploid);
        if (debug)
            PrintPaths(logger, threads, dir / "state_dump", "mk2000", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 7 / 2, 10000000));
        ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, readStorage);
        ManyKCorrect(logger, threads, dbg, readStorage, threshold, reliable_coverage, 3500, 3, diploid);
        ErrorCorrectionEngine(tournamentPathCorrector).run(logger, threads, dbg, readStorage);
        if (diploid)
            ErrorCorrectionEngine(bpCorrector).run(logger, threads, dbg, readStorage);
        std::vector<dbg::GraphPath> pseudo_reads = PartialRR(logger, threads, dbg, readStorage);
        printGraphAlignments(dir / "pseudo_reads.fasta", pseudo_reads);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        coverageStats(logger, dbg);
        if (debug)
            PrintPaths(logger, threads, dir / "state_dump", "mk3500", dbg, readStorage, paths_lib, false);
        readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");

        if (debug)
            DrawSplit(Component(dbg), dir / "split");
//    dbg.printFastaOld(dir / "final_dbg.fasta");
        printGFA(dir / "final_dbg.gfa", Component(dbg), true, &ag::SaveEdgeName<DBGTraits>);
        printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler());
        std::experimental::filesystem::path res;
        res = dir / "corrected_reads.fasta";
        logger.info() << "Initial correction results with k = " << k << " printed to " << res << std::endl;
        return {{"corrected_reads", res},
                {"pseudo_reads",    dir / "pseudo_reads.fasta"},
                {"final_dbg",       dir / "final_dbg.gfa"}};
    }

    class CoverageCorrectionStage : public Stage {
    public:
        CoverageCorrectionStage() : Stage(AlgorithmParameters(
                {"k-mer-size=501", "window=2000", "coverage-threshold=3", "reliable-coverage=10", "diploid", "load"},
                {}, ""), {"reads", "pseudo_reads", "paths"}, {"corrected_reads", "pseudo_reads", "final_dbg"}) {
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
                              input.find("paths")->second, threads, k, w, threshold, reliable_coverage, diploid, debug,
                              load);
        }
    };
}