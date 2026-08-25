#pragma once

#include <dbg/graph_algorithms.hpp>
#include <dbg/dbg_construction.hpp>
#include <dbg/dbg_read_alignment_storage.hpp>
#include <dbg/graph_printing.hpp>
#include "gap_closing.hpp"
#include "mult_correction.hpp"
#include "mitochondria_rescue.hpp"
#include "read_cleaning.hpp"
#include "dbg/path_dumping.hpp"

std::unordered_map<std::string, std::experimental::filesystem::path>
TopologyEC(logging::Logger &logger, const std::experimental::filesystem::path &dir,
        const io::Library &graph_lib, const io::Library &read_alignments_lib,
        const io::Library &paths_lib, const io::Library &references_lib, size_t threads, size_t k, double threshold, double reliable_coverage,
        size_t unique_threshold, bool diploid, bool debug, bool no_topology_correction) {
    logger.info() << "Performing topology-based error correction using k = " << k << std::endl;
    if (k%2==0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1)
                      << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k);
    dbg::SparseDBG dbg = dbg::LoadDBGFromEdgeSequences(logger, threads, graph_lib, hasher);
    ag::Printer gfa_printer;
    ag::Printer dot_printer(ag::VertexPrintStyles::defaultDotInfo(), ag::EdgePrintStyles::defaultDotInfo());
    dbg::DBGAlignedReadStorage readStorage = dbg::DBGAlignedReadStorage::Load(logger, threads, read_alignments_lib.front(), dbg, true);
    if(debug) {
        readStorage.logReads(threads, dir/"read_log.txt");
        readStorage.logGraph(dbg, logger.getLoggerStream(logging::LogLevel::trace));
    }
    size_t extension_size = 10000000;
    readStorage.trackSuffixes(logger, threads, dbg, 0, extension_size);
    dbg::DBGAlignedReadStorage refStorage(logger, threads, dbg, std::vector<ag::AlignedRead>(), false); //0, extension_size, false, false);
    if(debug) {
        ag::Printer printer(ag::EdgeInfo::Labeler(readStorage.getSuffixes().labeler()));
        printer.DrawSplit(ag::Component(dbg), dir/"before_figs", 25000);
        PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, references_lib, false);
    }
    initialCorrect(logger, threads, dbg, dir / "correction.txt", readStorage, refStorage,
                   threshold, 2 * threshold, reliable_coverage, diploid, 60000, false);
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "low", dbg, readStorage, paths_lib, references_lib, false);
    GapCloserPipeline(logger, threads, dbg);
//    readStorage.checkConsistency();
    if(debug) PrintPaths(logger, threads, dir/ "state_dump", "gap1", dbg, readStorage, paths_lib, references_lib, false);
    ag::AlignedReadStorage extra_reads;
    if (!no_topology_correction) {
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

        extra_reads = MultCorrect(logger, threads, dbg, dir / "mult2", readStorage, unique_threshold, 0, diploid, debug);
        MRescue(logger, threads, dbg, readStorage, unique_threshold, 0.05);
        if(debug) PrintPaths(logger, threads, dir/ "state_dump", "mult2", dbg, readStorage, paths_lib, references_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage.getReads(), &extra_reads, &refStorage.getReads()});
        if(debug) PrintPaths(logger, threads, dir/ "state_dump", "uncovered2", dbg, readStorage, paths_lib, references_lib, false);
        GapCloserPipeline(logger, threads, dbg);
    }
    dbg.resetEdgeCodes(logger, threads);
    if(debug) readStorage.checkConsistency();
    if(debug) {
        PrintPaths(logger, threads, dir / "state_dump", "gap2", dbg, readStorage, paths_lib, references_lib, false);
        ag::Printer printer(ag::EdgeInfo::Labeler(readStorage.getSuffixes().labeler()));
        printer.DrawSplit(ag::Component(dbg), dir/"split_figs", 25000);
    }
    printFasta(dir / "final_dbg.fasta", dbg, &ag::GetEdgeNameForSaving);
    gfa_printer.setEdgeInfo(ag::EdgeInfo({&ag::GetEdgeNameForSaving},{}, {}));
    gfa_printer.printGFA(dir / "final_dbg.gfa", ag::Component(dbg), true);
    gfa_printer.setEdgeInfo(ag::EdgeInfo({readStorage.getSuffixes().labeler()},{},{}));
    dot_printer += ag::EdgeInfo::Tooltiper(readStorage.getSuffixes().labeler());
    dot_printer.printDot(dir / "final_dbg.dot", ag::Component(dbg));
    readStorage.Save(dir/"final_dbg.aln");
    extra_reads.Save(dir / "extra_read.aln");
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
            {"k-mer-size=5001", "coverage-threshold=3", "reliable-coverage=10", "unique-threshold=40000", "diploid", "no-topology-correction"},
            {}, ""), {"graph", "read_alignments", "paths", "references"},
                                      {"corrected_reads", "pseudo_reads", "final_dbg", "final_aln", "extra_read_aln"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
        double reliable_coverage = std::stod(parameterValues.getValue("reliable-coverage"));
        double threshold = std::stod(parameterValues.getValue("coverage-threshold"));
        size_t unique_threshold = std::stoull(parameterValues.getValue("unique-threshold"));
        bool diploid = parameterValues.getCheck("diploid");
        bool no_topology_correction = parameterValues.getCheck("no-topology-correction");
        return TopologyEC(logger, dir, input.find("graph")->second, input.find("read_alignments")->second, input.find("paths")->second,
                          input.find("references")->second, threads, k, threshold, reliable_coverage, unique_threshold, diploid, debug,
                          no_topology_correction);
    }
};