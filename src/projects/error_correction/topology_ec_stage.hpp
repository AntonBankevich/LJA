#pragma once

#include <dbg/aln_reads_reader.hpp>

#include "gap_closing.hpp"
#include "mult_correction.hpp"
#include "mitochondria_rescue.hpp"

std::unordered_map<std::string, std::experimental::filesystem::path>
TopologyEC(logging::Logger &logger, const std::experimental::filesystem::path &dir,
        const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
        const io::Library &paths_lib, size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
           const std::string &reliability_mode, double ml_threshold,
        size_t unique_threshold, bool diploid, bool debug, bool load) {
    logger.info() << "Performing topology-based error correction using k = " << k << std::endl;
    if (k%2==0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1)
                      << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k);
    io::Library construction_lib = reads_lib + pseudo_reads_lib;
    SparseDBG dbg =
            load ? DBGPipeline(logger, hasher, w, construction_lib, dir, threads,
                               (dir/"disjointigs.fasta").string(),
                               (dir/"vertices.save").string())
                 : DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
    size_t extension_size = 10000000;
    ag::ReadLogger readLogger(threads, dir/"read_log.txt");
    dbg::ReadAlignmentStorage readStorage(dbg, 0, extension_size, true, debug);
    readStorage.setReadLogger(readLogger);
    dbg::ReadAlignmentStorage refStorage(dbg, 0, extension_size, false, false);
    refStorage.setReadLogger(readLogger);
    dbg::SeqReader reader(reads_lib, logger, threads);
    {
        KmerIndex index(dbg);
        index.fillAnchors(logger, threads, dbg, w);
        readStorage.FillAlignments(logger, threads, reader.begin(), reader.end(), dbg, index);
    }
    printDot(dir / "initial_dbg.dot", Component(dbg), ag::SaveEdgeName<DBGTraits>);
    printGFA(dir / "initial_dbg.gfa", Component(dbg), true, &ag::SaveEdgeName<DBGTraits>);
    AbstractReliableFillingAlgorithm * reliableFiller = nullptr;
    if(reliability_mode == "default") {
        reliableFiller = new CompositeReliableFiller(CreateDefaultReliableFiller(dbg, readStorage, reliable_coverage, diploid));
    } else {
        reliableFiller = new MLReliableFiller(reliability_mode, dir, ml_threshold);
    }
    if(debug) {
        DrawSplit(Component(dbg), dir / "before_figs", readStorage.labeler(), 25000);
        PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, false);
    }
    initialCorrect(logger, threads, dbg, dir / "correction.txt", readStorage, refStorage,
                   threshold, 2 * threshold, reliable_coverage, *reliableFiller, diploid, 60000, false);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "low", dbg, readStorage, paths_lib, false);
    GapCloserPipeline(logger, threads, dbg, {&readStorage, &refStorage});
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "gap1", dbg, readStorage, paths_lib, false);
    //    MultCorrect(logger, threads, dbg, dir / "mult1", readStorage, unique_threshold, 40, diploid, debug);
    //    if(debug) PrintPaths(logger, dir/ "state_dump", "mult1", dbg, readStorage, paths_lib, false);
    readStorage.delayedInvalidateBad(logger, threads, 1.01, "after_gap1_1");
    readStorage.applyCorrections(logger, threads);
    RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
    readStorage.delayedInvalidateBad(logger, threads, threshold, "after_gap1_threshold");
    readStorage.applyCorrections(logger, threads);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "bad", dbg, readStorage, paths_lib, false);
    RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "uncovered1", dbg, readStorage, paths_lib, false);
    MultCorrect(logger, threads, dbg, dir / "mult1", readStorage, unique_threshold, 0, diploid, debug);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "mult1", dbg, readStorage, paths_lib, false);
    RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});

    dbg::ReadAlignmentStorage extra_reads = MultCorrect(logger, threads, dbg, dir / "mult2", readStorage, unique_threshold, 0, diploid, debug);
    MRescue(logger, threads, dbg, readStorage, unique_threshold, 0.05);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "mult2", dbg, readStorage, paths_lib, false);
    RemoveUncovered(logger, threads, dbg, {&readStorage, &extra_reads, &refStorage});
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "uncovered2", dbg, readStorage, paths_lib, false);
    GapCloserPipeline(logger, threads, dbg, {&readStorage, &extra_reads, &refStorage});
    if(debug) {
        PrintPaths(logger, threads, dir / "state_dump", "gap2", dbg, readStorage, paths_lib, false);
        DrawSplit(Component(dbg), dir / "split_figs", readStorage.labeler());
    }
    printFasta(dir / "final_dbg.fasta", dbg, &ag::SaveEdgeName<DBGTraits>);
    printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler());
    printGFA(dir / "final_dbg.gfa", Component(dbg), true, &ag::SaveEdgeName<DBGTraits>);
    ag::SaveAllReads<DBGTraits>(dir/"final_dbg.aln", {&readStorage, &extra_reads});
    //readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");
    readStorage.printReadPaths(logger, dir / "corrected_reads.aln",
                                   dir / "final_dbg.gfa", dir / "corrected_reads.paths", k);
    extra_reads.printReadFasta(logger, dir / "pseudo_reads.fasta");
    std::experimental::filesystem::path res;
    res = dir / "corrected_reads.paths";
    logger.info() << "Second phase results with k = " << k << " printed to "
                  << res << std::endl;
    delete reliableFiller;
    return {{"corrected_reads", res}, {"pseudo_reads", dir / "pseudo_reads.fasta"}, {"final_dbg", dir / "final_dbg.gfa"}, {"final_aln", dir / "final_dbg.aln"}};
}


class TopologyCorrectionStage : public Stage {
public:
    TopologyCorrectionStage() : Stage(AlgorithmParameters(
            {"k-mer-size=5001", "window=500", "coverage-threshold=3", "reliable-coverage=10", "unique-threshold=40000", "reliability-mode=default", "ml-threshold=0.5", "diploid", "load"},
            {}, ""), {"reads", "pseudo_reads", "paths"}, {"corrected_reads", "pseudo_reads", "final_dbg", "final_aln"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
        size_t w = std::stoi(parameterValues.getValue("window"));
        double reliable_coverage = std::stod(parameterValues.getValue("reliable-coverage"));
        double threshold = std::stod(parameterValues.getValue("coverage-threshold"));
        size_t unique_threshold = std::stoull(parameterValues.getValue("unique-threshold"));
        bool diploid = parameterValues.getCheck("diploid");
        bool load = parameterValues.getCheck("load");
        double ml_threshold = std::stod(parameterValues.getValue("ml-threshold"));
        return TopologyEC(logger, dir, input.find("reads")->second, input.find("pseudo_reads")->second, input.find("paths")->second,
                          threads, k, w, threshold, reliable_coverage, parameterValues.getValue("reliability-mode"), ml_threshold, unique_threshold, diploid, debug, load);
    }
};
