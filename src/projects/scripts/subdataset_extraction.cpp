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
                     "reference=none", "compress", "dimer-compress=1000000000,1000000000,1", "unique-threshold=40000", "radius=1000", "bad-cov=7", "track-paths"},
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
    const std::experimental::filesystem::path dir(parser.getValue("output-dir")); //initialization of dir
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
    bool track_paths = parser.getCheck("track-paths");
    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library paths_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));
    io::Library ref_lib;
    if(parser.getValue("reference") != "none")
        ref_lib =  oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reference"));
    std::string disjointigs_file = parser.getValue("disjointigs");
    std::string vertices_file = parser.getValue("vertices");
    std::string dbg_file = parser.getValue("dbg");
    hashing::RollingHash hasher(k, std::stoi(parser.getValue("base")));
    size_t threads = std::stoi(parser.getValue("threads"));
    dbg::SparseDBG dbg = dbg_file == "none" ?
                    DBGPipeline(logger, hasher, w, reads_lib, dir, threads, disjointigs_file, vertices_file) :
                    dbg::LoadDBGFromFasta({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads); //Create dbg
    dbg.fillAnchors(w, logger, threads);// Creates index of edge k-mers
    size_t extension_size = 100000;
    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, false, track_paths);//Structure for read alignments
    io::SeqReader reader(reads_lib);//Reader that can read reads from file
    readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);//Align reads to the graph
    std::experimental::filesystem::path subdir = dir / "subdatasets";
    recreate_dir(subdir);
    std::vector<Subdataset> subdatasets;
    GraphAlignmentStorage storage(dbg);
    for(StringContig stringContig : io::SeqReader(ref_lib)) {
        storage.addContig(stringContig.makeContig());
    }
    if(paths_lib.empty()) {
        logger.info() << "No paths provided. Splitting the whole graph." << std::endl;
//        std::function<bool(const dbg::Component&)> f = [bad_cov](const dbg::Component &component) {
//            for(dbg::Edge &edge : component.edgesInnerUnique()) {
//                if(edge.getCoverage() >= 2 && edge.getCoverage() < bad_cov) {
//                    return false;
//                    break;
//                }
//            }
//            return true;
//        };
//        std::vector<dbg::Component> components = oneline::filter(dbg::LengthSplitter(unique_threshold).splitGraph(dbg), f);
        std::vector<dbg::Component> components = dbg::LengthSplitter(unique_threshold).splitGraph(dbg); //Split graph into components
        subdatasets = oneline::initialize<Subdataset>(components);//Create subdatasets corresponding to components
    } else {
        logger.info() << "Extracting subdatasets around contigs" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        io::SeqReader reader(paths_lib);
        size_t radius = std::stoull(parser.getValue("radius"));
        for(StringContig scontig : io::SeqReader(paths_lib)) {
            Contig contig = scontig.makeContig();
            std::cout << contig.id << " " << contig.size() << " " << dbg::GraphAligner(dbg).carefulAlign(contig).size() << std::endl;
            storage.addContig(contig);
            subdatasets.emplace_back(dbg::Component::neighbourhood(dbg, contig, dbg.hasher().getK() + radius));
            subdatasets.back().id = contig.id;
        }
    }
    storage.Fill(threads);
    FillSubdatasets(subdatasets, {&readStorage}, true);//Assign reads to datasets
    size_t cnt = 0;
    for(const Subdataset &subdataset: subdatasets) {//Print subdatasets to disk
        logger.info() << "Printing subdataset " << cnt << " " << subdataset.id << ":";
        for(dbg::Vertex &v : subdataset.component.verticesUnique()) {
            logger << " " << v.getShortId();
        }
        logger << "\n";
        std::string name = itos(cnt);
        if(!subdataset.id.empty())
            name += "_" + name;
        subdataset.Save(subdir / name, storage.labeler() + readStorage.labeler());
        cnt++;
    }
    return 0;
}
