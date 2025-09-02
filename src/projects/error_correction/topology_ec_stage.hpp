#pragma once

#include <dbg/aln_reads_reader.hpp>
#include "gap_closing.hpp"
#include "mult_correction.hpp"
#include "mitochondria_rescue.hpp"
#include "read_cleaning.hpp"

std::unordered_map<std::string, std::experimental::filesystem::path>
TopologyEC(logging::Logger &logger, const std::experimental::filesystem::path &dir,
        const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
        const io::Library &paths_lib, const io::Library &references_lib, size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
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
                               (dir/"vertices.save").string(), debug)
                 : DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
    size_t extension_size = 10000000;
    dbg::SeqReader reader(reads_lib, logger, threads);
    dbg::DBGAlignedReadStorage readStorage(logger, threads, dbg,
                                           AlignReads(logger, threads, reader.begin(), reader.end(), dbg, w),
                                           true);
    if(debug) readStorage.logReads(threads, dir/"read_log.txt");
    readStorage.trackSuffixes(logger, threads, dbg, 0, extension_size);
    dbg::DBGAlignedReadStorage refStorage(logger, threads, dbg, std::vector<ag::AlignedRead<DBGTraits>>(), false); //0, extension_size, false, false);
    Printer<dbg::DBGTraits> printer;
    printer.setEdgeInfo(ObjInfo<dbg::Edge>({&ag::GetEdgeNameForSaving<DBGTraits>}, {}, {}));
    printer.printDot(dir / "initial_dbg.dot", Component(dbg));
    //printDot(dir / "initial_dbg.dot", Component(dbg), ag::GetEdgeNameForSaving<DBGTraits>);
    if(debug) {
        DrawSplit(Component(dbg), dir / "before_figs", readStorage.getSuffixes().labeler(), 25000);
        PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, references_lib, false);
    }
    initialCorrect(logger, threads, dbg, dir / "correction.txt", readStorage, refStorage,
                   threshold, 2 * threshold, reliable_coverage, diploid, 60000, false);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "low", dbg, readStorage, paths_lib, references_lib, false);
    GapCloserPipeline(logger, threads, dbg);
//    readStorage.checkConsistency();
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "gap1", dbg, readStorage, paths_lib, references_lib, false);
    //    MultCorrect(logger, threads, dbg, dir / "mult1", readStorage, unique_threshold, 40, diploid, debug);
    //    if(debug) PrintPaths(logger, dir/ "state_dump", "mult1", dbg, readStorage, paths_lib, references_lib, false);
    InvalidateLowCovered(logger, threads, readStorage.getReads(), 1.01, 500, "after_gap1_1");
    readStorage.getReads().applyCorrections(logger, threads);
    RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
    InvalidateLowCovered(logger, threads, readStorage.getReads(), threshold, 500, "after_gap1_1");
    readStorage.getReads().applyCorrections(logger, threads);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "bad", dbg, readStorage, paths_lib, references_lib, false);
    RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "uncovered1", dbg, readStorage, paths_lib, references_lib, false);
    MultCorrect(logger, threads, dbg, dir / "mult1", readStorage, unique_threshold, 0, diploid, debug);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "mult1", dbg, readStorage, paths_lib, references_lib, false);
    RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &refStorage.getReads()});

    ag::AlignedReadStorage<DBGTraits> extra_reads = MultCorrect(logger, threads, dbg, dir / "mult2", readStorage, unique_threshold, 0, diploid, debug);
    MRescue(logger, threads, dbg, readStorage, unique_threshold, 0.05);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "mult2", dbg, readStorage, paths_lib, references_lib, false);
    RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &extra_reads, &refStorage.getReads()});
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "uncovered2", dbg, readStorage, paths_lib, references_lib, false);
    GapCloserPipeline(logger, threads, dbg);
    dbg.resetEdgeCodes(logger, threads);
    if(debug) {
        PrintPaths(logger, threads, dir / "state_dump", "gap2", dbg, readStorage, paths_lib, references_lib, false);
        DrawSplit(Component(dbg), dir / "split_figs", readStorage.getSuffixes().labeler());
    }
    printFasta(dir / "final_dbg.fasta", dbg, &ag::GetEdgeNameForSaving<DBGTraits>);
    printer.setEdgeInfo(ObjInfo<Edge>({&ag::GetEdgeNameForSaving<DBGTraits>},{}, {}));
    printer.printGFA(dir / "final_dbg.gfa", Component(dbg), true);
    printer.setEdgeInfo(ObjInfo<Edge>({readStorage.getSuffixes().labeler()},{},{}));
    printer.printDot(dir / "final_dbg.dot", Component(dbg));
    //printDot(dir / "final_dbg.dot", Component(dbg), readStorage.getSuffixes().labeler()); delete if ok
    //printGFA(dir / "final_dbg.gfa", Component(dbg), true, &ag::SaveEdgeName<DBGTraits>); delete if ok
    ag::SaveReads(dir/"final_dbg.aln", readStorage);
    ag::SaveReads(dir / "extra_read.aln", extra_reads);
    //readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");
    readStorage.getReads().printReadPaths(logger, dir / "corrected_reads.aln",
                                   dir / "final_dbg.gfa", dir / "corrected_reads.paths", k);
    extra_reads.printReadFasta(logger, dir / "pseudo_reads.fasta");
    std::experimental::filesystem::path res;
    res = dir / "corrected_reads.paths";
    logger.info() << "Second phase results with k = " << k << " printed to "
                  << res << std::endl;
    return {{"corrected_reads", res}, {"pseudo_reads", dir / "pseudo_reads.fasta"},
            {"final_dbg", dir / "final_dbg.gfa"}, {"final_aln", dir / "final_dbg.aln"},
            {"extra_read_aln", dir / "extra_read.aln"}};
}


class TopologyCorrectionStage : public Stage {
public:
    TopologyCorrectionStage() : Stage(AlgorithmParameters(
            {"k-mer-size=5001", "window=500", "coverage-threshold=3", "reliable-coverage=10", "unique-threshold=40000", "diploid", "load"},
            {}, ""), {"reads", "pseudo_reads", "paths", "references"},
                                      {"corrected_reads", "pseudo_reads", "final_dbg", "final_aln", "extra_read_aln"}) {
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
        return TopologyEC(logger, dir, input.find("reads")->second, input.find("pseudo_reads")->second, input.find("paths")->second,
                          input.find("references")->second, threads, k, w, threshold, reliable_coverage, unique_threshold, diploid, debug, load);
    }
};