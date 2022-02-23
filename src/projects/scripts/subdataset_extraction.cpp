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
    CLParser parser({"vertices=none", "unique=none", "dbg=none", "output-dir=",
                     "threads=16", "k-mer-size=", "window=2000", "base=239", "debug", "disjointigs=none",
                     "reference=none", "compress", "dimer-compress=1000000000,1000000000,1", "unique-threshold=40000", "bad-cov=7"},
                    {"reads"},
                    {"o=output-dir", "t=threads", "k=k-mer-size", "w=window"},
                    "");
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parser.check() << "\n" << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }

    bool debug = parser.getCheck("debug");
    StringContig::homopolymer_compressing = parser.getCheck("compress");
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    const size_t w = std::stoi(parser.getValue("window"));
    double bad_cov = std::stod(parser.getValue("bad-cov"));
    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    std::string disjointigs_file = parser.getValue("disjointigs");
    std::string vertices_file = parser.getValue("vertices");
    std::string dbg_file = parser.getValue("dbg");
    hashing::RollingHash hasher(k, std::stoi(parser.getValue("base")));
    size_t threads = std::stoi(parser.getValue("threads"));
    dbg::SparseDBG dbg = dbg_file == "none" ?
                    DBGPipeline(logger, hasher, w, reads_lib, dir, threads, disjointigs_file, vertices_file) :
                    dbg::LoadDBGFromFasta({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);
    dbg.fillAnchors(w, logger, threads);
    size_t extension_size = 100000;
    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true, false);
    io::SeqReader reader(reads_lib);
    readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
    std::experimental::filesystem::path subdir = dir / "subdatasets";
    recreate_dir(subdir);
    std::vector<Subdataset> all = SubdatasetSplit(dbg, {&readStorage}, unique_threshold, true);
    size_t cnt = 0;
    for(const Subdataset &subdataset: all) {
        bool ok = false;
        for(dbg::Edge &edge : subdataset.component.edgesInnerUnique()) {
            if(edge.getCoverage() >= 2 && edge.getCoverage() < bad_cov) {
                ok = true;
                break;
            }
        }
        if(!ok)
            continue;
        logger.info() << "Printing subdataset " << cnt << ":";
        for(dbg::Vertex &v : subdataset.component.verticesUnique()) {
            logger << " " << v.getShortId();
        }
        logger << "\n";
        subdataset.Save(subdir / itos(cnt));
        cnt++;
    }
    return 0;
}
