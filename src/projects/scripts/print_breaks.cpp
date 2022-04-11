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
                     "reference=none", "compress", "dimer-compress=1000000000,1000000000,1", "unique-threshold=40000", "radius=6000", "bad-cov=7"},
                    {"paths", "reads"},
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
    io::Library paths_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));
    std::string disjointigs_file = parser.getValue("disjointigs");
    std::string vertices_file = parser.getValue("vertices");
    std::string dbg_file = parser.getValue("dbg");
    hashing::RollingHash hasher(k, std::stoi(parser.getValue("base")));
    size_t threads = std::stoi(parser.getValue("threads"));
    dbg::SparseDBG dbg = dbg_file == "none" ?
                    DBGPipeline(logger, hasher, w, reads_lib, dir, threads, disjointigs_file, vertices_file) :
                    dbg::LoadDBGFromFasta({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);
    size_t radius = std::stoull(parser.getValue("radius"));
    dbg.fillAnchors(100, logger, threads);
    radius -= std::min(radius, k);
    std::ofstream os;
    os.open(dir / "breaks.fasta");
    dbg::GraphAligner aligner(dbg);
    for(StringContig sc : io::SeqReader(paths_lib)) {
        Contig contig = sc.makeContig();
        std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> al = aligner.carefulAlign(contig);
        for(size_t i = 0; i + 1 < al.size(); i++) {
            if(al[i].seg_to.contig() == al[i - 1].seg_to.contig() || al[i].seg_from.right + k + 5000 < al[i + 1].seg_from.left)
                continue;
            if(al[i].seg_to.contig().end()->outDeg() == 0 || al[i + 1].seg_to.contig().start()->inDeg() == 0) {
                logger.trace() << "Possible break " << al[i] << " " << al[i + 1] << std::endl;
                Segment<Contig> seg = al[i].seg_from.unite(al[i + 1].seg_from).extendRight(k).extendBy(radius);
                logger.trace() << "Printing segment " << seg << std::endl;
                os << ">" << seg.contig().id << "_" << seg.left << "_" << seg.right << "\n" << seg.seq() << "\n";
            }
        }
    }
    os.close();
    return 0;
}
