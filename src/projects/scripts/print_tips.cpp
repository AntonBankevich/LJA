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

int main(int argc, char **argv) {
    AlgorithmParameters params({"vertices=none", "unique=none", "dbg=none", "output-dir=",
                                "threads=16", "k-mer-size=", "window=2000", "debug", "disjointigs=none",
                                "reference=none", "compress", "dimer-compress=1000000000,1000000000,1", "unique-threshold=40000", "radius=6000", "bad-cov=7"},
                               {"paths", "reads"}, "");
    CLParser parser(params, {"o=output-dir", "t=threads", "k=k-mer-size", "w=window"}, {});
    AlgorithmParameterValues parameterValues = parser.parseCL(argc, argv);
    if (!parameterValues.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parameterValues.checkMissingValues() << "\n" << std::endl;
        std::cout << parameterValues.helpMessage() << std::endl;
        return 1;
    }

    bool debug = parameterValues.getCheck("debug");
    StringContig::homopolymer_compressing = parameterValues.getCheck("compress");
    StringContig::SetDimerParameters(parameterValues.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(parameterValues.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
    const size_t w = std::stoi(parameterValues.getValue("window"));
    double bad_cov = std::stod(parameterValues.getValue("bad-cov"));
    size_t unique_threshold = std::stoi(parameterValues.getValue("unique-threshold"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("reads"));
    io::Library paths_lib = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("paths"));
    std::string disjointigs_file = parameterValues.getValue("disjointigs");
    std::string vertices_file = parameterValues.getValue("vertices");
    std::string dbg_file = parameterValues.getValue("dbg");
    hashing::RollingHash hasher(k);
    size_t threads = std::stoi(parameterValues.getValue("threads"));
    dbg::SparseDBG dbg = dbg_file == "none" ?
                    DBGPipeline(logger, hasher, w, reads_lib, dir, threads, disjointigs_file, vertices_file) :
                         dbg::LoadDBGFromEdgeSequences(logger, threads, {std::experimental::filesystem::path(dbg_file)}, hasher);
    size_t radius = std::stoull(parameterValues.getValue("radius"));
    radius -= std::min(radius, k);
    std::ofstream os;
    os.open(dir / "tips.fasta");
    for(dbg::Edge &edge : dbg.edges()) {
        if(edge.getStart().inDeg() == 0) {
            os << ">" << edge.getShortId() << "\n" << edge.getStart().getSeq() << edge.truncSeq().Subseq(0, std::min(
                    edge.truncSize(), radius)) << "\n";
        }
    }
    os.close();
    return 0;
}
