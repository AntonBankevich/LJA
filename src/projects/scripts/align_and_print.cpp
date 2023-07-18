//
// Created by anton on 24.01.2023.
//

#include <common/cl_parser.hpp>
#include <sequences/contigs.hpp>
#include <experimental/filesystem>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/rolling_hash.hpp>
#include <dbg/dbg_construction.hpp>
#include <dbg/component.hpp>
#include <dbg/graph_alignment_storage.hpp>
#include <dbg/subdatasets.hpp>
#include "dbg_graph_aligner.hpp"

int main(int argc, char **argv) {
    AlgorithmParameters parameters({"vertices=none", "unique=none", "dbg=", "output-dir=",
                        "threads=16", "k-mer-size=", "window=2000", "base=239", "debug", "disjointigs=none",
                        "compress", "dimer-compress=32,32,1"}, {"reads", "paths"}, "");

    CLParser parser(parameters, {"o=output-dir", "t=threads", "k=k-mer-size", "w=window"});
    AlgorithmParameterValues values = parser.parseCL(argc, argv);
    if (!values.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << values.checkMissingValues() << "\n" << std::endl;
        std::cout << values.helpMessage() << std::endl;
        return 1;
    }

    bool debug = values.getCheck("debug");
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(values.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(values.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    size_t k = std::stoi(values.getValue("k-mer-size"));
    const size_t w = std::stoi(values.getValue("window"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(values.getListValue("reads"));
    io::Library paths_lib = oneline::initialize<std::experimental::filesystem::path>(values.getListValue("paths"));
    std::string disjointigs_file = values.getValue("disjointigs");
    std::string vertices_file = values.getValue("vertices");
    std::string dbg_file = values.getValue("dbg");
    hashing::RollingHash hasher(k, std::stoi(values.getValue("base")));
    size_t threads = std::stoi(values.getValue("threads"));
    dbg::SparseDBG dbg = dbg_file == "none" ?
                         DBGPipeline(logger, hasher, w, reads_lib, dir, threads, disjointigs_file, vertices_file) :
                         dbg::LoadDBGFromEdgeSequences({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);
    dbg.fillAnchors(w, logger, threads);
    logger.info() << "Constructing getEdge id mapping" << std::endl;
    std::unordered_map<dbg::Edge*, std::string> edge_mapping;
    dbg::GraphAligner aligner(dbg);
    io::SeqReader graphReader(dbg_file);
    for(StringContig s : graphReader) {
        Contig rseq = s.makeContig();
        DBGGraphPath al = aligner.align(rseq.getSeq());
        VERIFY(al.size() == 1);
        edge_mapping[&al.front().contig()] = s.id;
        edge_mapping[&al.front().contig().rc()] = "-" + s.id;
    }
    std::ofstream os;
    logger.info() << "Performing alignments" << std::endl;
    os.open(dir / "alignments.txt");
    io::SeqReader reader(paths_lib);
    for(StringContig read : reader) {
        Contig rseq = read.makeContig();
        os << ">" << read.id << "\n";
        std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> al = aligner.carefulAlign(rseq);
        for(dbg::PerfectAlignment<Contig, dbg::Edge> &piece : al) {
            os << edge_mapping[&piece.seg_to.contig()] << " ";
        }
        os << "\n";
    }
    os.close();
    logger.info() << "Finished alignments. Resulting paths printed to " << (dir / "alignments.txt") << std::endl;
    return 0;
}
