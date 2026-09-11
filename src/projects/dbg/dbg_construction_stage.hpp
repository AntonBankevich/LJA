#pragma once

#include <common/pipeline_tools.hpp>
#include "dbg_construction.hpp"
#include "dbg_read_alignment_storage.hpp"
#include "graph_printing.hpp"
#include "aln_reads_reader.hpp"

namespace dbg {
    std::unordered_map<std::string, std::experimental::filesystem::path>
    ConstructDBG(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                 const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
                 size_t threads, size_t k, size_t w, bool debug) {
        logger.info() << "Constructing de Bruijn graph with k = " << k << std::endl;
        if (k % 2 == 0) {
            logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
            k += 1;
        }
        ensure_dir_existance(dir);
        hashing::RollingHash hasher(k);
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        dbg::SparseDBG dbg = DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
        ag::Printer gfa_printer;
        ag::Printer dot_printer(ag::VertexPrintStyles::defaultDotInfo(), ag::EdgePrintStyles::defaultDotInfo());
        gfa_printer.setEdgeInfo(ag::EdgeInfo({&ag::GetEdgeNameForSaving}, {}, {}));
        logger.info() << "Printing graph to " << (dir / "initial_dbg.gfa") << std::endl;
        gfa_printer.printGFA(dir / "initial_dbg.gfa", dbg);
        logger.info() << "Finished printing graph" << std::endl;
        dot_printer.printDot(dir / "initial_dbg.dot", dbg);
        std::experimental::filesystem::path al_file = dir / "initial_alignments.aln";
        dbg::SeqReader reader(reads_lib, logger, threads);
        dbg::DBGAlignedReadStorage readStorage(logger, threads, dbg,
                AlignReads(logger, threads, reader.begin(), reader.end(), dbg, w), true);
        logger.info() << "Printing read alignments to " << al_file << std::endl;
        readStorage.Save(al_file);
        logger.info() << "Finished printing read alignments to " << al_file << std::endl;
        if(debug) {
            readStorage.logReads(threads, dir / "read_log.txt");
            readStorage.logGraph(dbg, logger.getLoggerStream(logging::LogLevel::trace));
        }
//        Pseudo-reads (e.g. draft contigs used only to help build the graph) are aligned separately from real
//        reads and never contribute to edge coverage, since they aren't sequencing observations.
        std::experimental::filesystem::path pseudo_al_file = dir / "initial_pseudo_alignments.aln";
        dbg::SeqReader pseudo_reader(pseudo_reads_lib, logger, threads);
        dbg::DBGAlignedReadStorage pseudoReadStorage(logger, threads, dbg,
                AlignReads(logger, threads, pseudo_reader.begin(), pseudo_reader.end(), dbg, w), false);
        logger.info() << "Printing pseudo-read alignments to " << pseudo_al_file << std::endl;
        pseudoReadStorage.Save(pseudo_al_file);
        logger.info() << "Finished printing pseudo-read alignments to " << pseudo_al_file << std::endl;
//        A well-formed, deliberately empty alignment file. Lets stages with no notion of "extra"
//        (repeat-resolution-derived) reads still satisfy consumers that always expect this input, e.g. Multiplexing.
        dbg::DBGAlignedReadStorage emptyStorage(logger, threads, dbg, std::vector<ag::AlignedRead>(), false);
        emptyStorage.Save(dir / "extra_reads.aln");
        return {{"graph",                  dir / "initial_dbg.gfa"},
                {"read_alignments",        al_file},
                {"pseudo_read_alignments", pseudo_al_file},
                {"extra_reads",            dir / "extra_reads.aln"}};
    }

    class DBGConstructionStage : public Stage {
    public:
        explicit DBGConstructionStage(size_t default_k = 501, size_t default_w = 2000) : Stage(AlgorithmParameters(
                {"k-mer-size=" + std::to_string(default_k), "window=" + std::to_string(default_w)},
                {}, ""), {"reads", "pseudo_reads"}, {"graph", "read_alignments", "pseudo_read_alignments", "extra_reads"}) {
        }

    protected:
        std::unordered_map<std::string, std::experimental::filesystem::path>
        innerRun(logging::Logger &logger, size_t threads,
                 const std::experimental::filesystem::path &dir, bool debug,
                 const AlgorithmParameterValues &parameterValues,
                 const std::unordered_map<std::string, io::Library> &input) override {
            size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
            size_t w = std::stoi(parameterValues.getValue("window"));
            return ConstructDBG(logger, dir, input.find("reads")->second, input.find("pseudo_reads")->second,
                                 threads, k, w, debug);
        }
    };
}
