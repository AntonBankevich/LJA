#pragma once

#include "repeat_resolution.hpp"

std::unordered_map<std::string, std::experimental::filesystem::path> MDBGConstruction(
        logging::Logger &logger, size_t threads, size_t k, size_t kmdbg, size_t unique_threshold, bool diploid,
        const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &graph_gfa,
        const std::experimental::filesystem::path &read_paths, bool debug) {
    logger.info() << "Performing repeat resolution by transforming de Bruijn graph into Multiplex de Bruijn graph" << std::endl;
    hashing::RollingHash hasher(k);
    SparseDBG dbg = dbg::LoadDBGFromEdgeSequences({graph_gfa}, hasher, logger, threads);
    IdIndex<Vertex> index(dbg.vertices().begin(), dbg.vertices().end());
    size_t extension_size = 10000000;
    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, debug);
    RecordStorage extra_reads(dbg, 0, extension_size, threads, readLogger, false, debug);
    LoadAllReads(read_paths, {&readStorage, &extra_reads}, index);
    repeat_resolution::RepeatResolver rr(dbg, &readStorage, {&extra_reads},
                                         k, kmdbg, dir, unique_threshold,
                                         diploid, debug, logger);
    return rr.ResolveRepeats(logger, threads);
}

class MDBGStage : public Stage {
public:
    MDBGStage() : Stage(AlgorithmParameters({"k-mer-size=5001", "max-k=40000", "unique-threshold=40000", "diploid"},
                  {}, ""), {"read_aln", "graph"}, {"graph", "assembly"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                                          const std::experimental::filesystem::path &dir, bool debug,
                                          const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
        size_t max_k = std::stoi(parameterValues.getValue("max-k"));
        size_t unique_threshold = std::stoi(parameterValues.getValue("unique-threshold"));
        bool diploid = parameterValues.getCheck("diploid");
        return MDBGConstruction(logger, threads, k, max_k, unique_threshold, diploid, dir,
                          input.find("graph")->second.front(), input.find("read_aln")->second.front(), debug);
    }
};