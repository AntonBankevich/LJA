#pragma once

#include "repeat_resolution.hpp"

std::unordered_map<std::string, std::experimental::filesystem::path> MDBGConstruction(
        logging::Logger &logger, size_t threads, size_t k, size_t kmdbg, size_t w, size_t unique_threshold, bool diploid,
        const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &graph_fasta,
        const std::experimental::filesystem::path &read_paths, bool skip, bool debug) {
    logger.info() << "Performing repeat resolution by transforming de Bruijn graph into Multiplex de Bruijn graph" << std::endl;
    hashing::RollingHash hasher(k);
    SparseDBG dbg = dbg::LoadDBGFromFasta({graph_fasta}, hasher, logger, threads);
    size_t extension_size = 10000000;
    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, debug);
    RecordStorage extra_reads(dbg, 0, extension_size, threads, readLogger, false, debug);
    LoadAllReads(read_paths, {&readStorage, &extra_reads}, dbg);
    repeat_resolution::RepeatResolver rr(dbg, &readStorage, {&extra_reads},
                                         k, kmdbg, dir, unique_threshold,
                                         diploid, debug, logger);
    return rr.ResolveRepeats(logger, threads);
}

class MDBGStage : public Stage {
public:
    MDBGStage() : Stage(AlgorithmParameters({"k-mer-size=5001", "max-k=40000", "window=500", "unique-threshold=40000", "diploid"},
                  {}, ""), {"read_aln", "graph"}, {"graph", "assembly"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                                                                  const std::experimental::filesystem::path &dir, bool debug,
                                                                                  const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
        size_t w = std::stoi(parameterValues.getValue("window"));
        double reliable_coverage = std::stod(parameterValues.getValue("reliable-coverage"));
        double threshold = std::stod(parameterValues.getValue("coverage-threshold"));
        bool diploid = parameterValues.getCheck("diploid");
        return CoverageEC(logger, dir, input.find("reads")->second, input.find("pseudo_reads")->second,
                          input.find("paths")->second, threads, k, w, threshold, reliable_coverage, diploid, debug, false);
    }
};