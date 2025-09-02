#pragma once

std::unordered_map<std::string, std::experimental::filesystem::path>
NoCorrection(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const io::Library &reads_lib,
             const io::Library &pseudo_reads_lib, const io::Library &paths_lib, size_t k, size_t w, bool debug, bool load) {
    logger.info() << "Performing graph construction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    hashing::RollingHash hasher(k);
    ensure_dir_existance(dir);
    io::Library construction_lib = reads_lib + pseudo_reads_lib;
    SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, construction_lib, dir, threads,
                                       (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string(), debug) :
                    DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
    size_t extension_size = std::max<size_t>(k * 2, 1000);
    dbg::SeqReader reader(reads_lib, logger, threads);
    dbg::DBGAlignedReadStorage readStorage(logger, threads, dbg,
                                           AlignReads(logger, threads, reader.begin(), reader.end(), dbg, w),
                                           true);
    readStorage.trackSuffixes(logger, threads, dbg, 0, extension_size);
    dbg::DBGAlignedReadStorage extra_reads(logger, threads, dbg, std::vector<ag::AlignedRead<DBGTraits>>(), false);
    coverageStats(logger, dbg);
    printFasta(dir / "final_dbg.fasta", dbg, &ag::GetEdgeNameForSaving<DBGTraits>);
    Printer<DBGTraits> printer;
    printer.addEdgeInfo(ObjInfo<Edge>({},{}, {}));
    printer.printGFA(dir / "final_dbg.gfa", Component(dbg), true);
    printer.setEdgeInfo(ObjInfo<Edge>({readStorage.getSuffixes().labeler()},{},{}));
    printer.printDot(dir / "final_dbg.dot", Component(dbg));
    //printGFA(dir / "final_dbg.gfa", Component(dbg), true, &ag::GetEdgeNameForSaving<DBGTraits>); delete if ok
    //printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler()); delete if ok
    ag::SaveReads(dir/"final_dbg.aln", readStorage);
    ag::SaveReads(dir / "extra_read.aln", extra_reads);
    //readStorage.getReads().printReadFasta(logger, dir / "corrected_reads.fasta");
    readStorage.getReads().printReadPaths(logger, dir / "corrected_reads.aln",
                                   dir / "final_dbg.gfa", dir / "corrected_reads.paths", k);
    return {{"corrected_reads", dir/"corrected_reads.paths"}, {"final_dbg", dir / "final_dbg.gfa"},
            {"final_aln", dir / "final_dbg.aln"}, {"extra_reads_aln", dir / "extra_reads.aln"}};
}

class NoCorrectionStage : public Stage {
public:
    NoCorrectionStage() : Stage(AlgorithmParameters({"k-mer-size=5001", "window=500", "load"}, {}, ""),
                                {"reads", "pseudo_reads", "paths"}, {"corrected_reads", "final_dbg", "final_aln"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                      const std::experimental::filesystem::path &dir, bool debug,
                      const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
        size_t w = std::stoi(parameterValues.getValue("window"));
        bool load = parameterValues.getCheck("load");
        return std::move(NoCorrection(logger, threads, dir, input.find("reads")->second, input.find("pseudo_reads")->second,
                                      input.find("paths")->second, k, w, debug, load));
    }
};

