#include <common/cl_parser.hpp>
#include <sequences/contigs.hpp>
#include <experimental/filesystem>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/rolling_hash.hpp>
#include <dbg/dbg_construction.hpp>
#include <assembly_graph/data_structures/component.hpp>
#include <dbg/dbg_read_alignment_storage.hpp>
#include <dbg/subdatasets.hpp>
#include <dbg/aln_reads_reader.hpp>
#include "dbg/dbg_graph_aligner.hpp"
#include "dbg/path_dumping.hpp"

using namespace dbg;
int main(int argc, char **argv) {
    AlgorithmParameters params({"unique=none", "dbg=none", "output-dir=",
                               "threads=16", "k-mer-size=", "window=2000", "debug",
                               "reference=none", "compress", "dimer-compress=1000000000,1000000000,1",
                               "unique-threshold=40000", "radius=1000", "bad-cov=7", "track-paths", "add-paths"},
                               {"paths", "reads", "pseudo-reads", "contigs"}, "");
    CLParser parser(params,
                    {"o=output-dir", "t=threads", "k=k-mer-size", "w=window"},
                    {});
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
    const std::experimental::filesystem::path dir(parameterValues.getValue("output-dir")); //initialization of dir
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
    bool track_paths = parameterValues.getCheck("track-paths");
    bool add_paths = parameterValues.getCheck("add-paths");
    size_t unique_threshold = std::stoi(parameterValues.getValue("unique-threshold"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("reads"));
    io::Library pseudo_reads_lib = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("pseudo-reads"));
    io::Library contig_lib = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("contigs"));
    io::Library construction_lib = reads_lib + pseudo_reads_lib;
    io::Library paths_lib = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("paths"));
    if(add_paths)
        construction_lib = construction_lib + paths_lib;
    io::Library ref_lib;
    if(parameterValues.getValue("reference") != "none")
        ref_lib =  oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("reference"));
    std::string dbg_file = parameterValues.getValue("dbg");
    hashing::RollingHash hasher(k);
    size_t threads = std::stoi(parameterValues.getValue("threads"));
    dbg::SparseDBG dbg = dbg_file == "none" ?
                    DBGPipeline(logger, hasher, w, construction_lib, dir, threads) :
                         dbg::LoadDBGFromEdgeSequences(logger, threads, {std::experimental::filesystem::path(dbg_file)}, hasher); //Create dbg
    size_t extension_size = 100000;
    io::SeqReader reader(reads_lib);//Reader that can read reads from file
    dbg::DBGAlignedReadStorage readStorage(logger, threads, dbg,
                                           AlignReads(logger, threads, reader.begin(), reader.end(), dbg, w),
                                           false);
    dbg::KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, w);
    readStorage.trackSuffixes(logger, threads, dbg, 0, 1000000);
    std::experimental::filesystem::path subdir = dir / "subdatasets";
    recreate_dir(subdir);
    std::vector<Subdataset> subdatasets;
    AlignedContigStorage storage(dbg);
    for(StringContig stringContig : dbg::SeqReader(ref_lib, logger, threads)) {
        storage.addContig(stringContig.makeContig());
    }
    if(paths_lib.empty()) {
        logger.info() << "No paths provided. Splitting the whole graph." << std::endl;
//        std::function<bool(const ag::Component&)> f = [bad_cov](const ag::Component &component) {
//            for(dbg::Edge &edge : component.edgesInnerUnique()) {
//                if(edge.getCoverage() >= 2 && edge.getCoverage() < bad_cov) {
//                    return false;
//                    break;
//                }
//            }
//            return true;
//        };
//        std::vector<ag::Component> components = oneline::filter(ag::LengthSplitter(unique_threshold).splitGraph(dbg), f);
        std::vector<ag::Component> components = ag::LengthSplitter(unique_threshold).splitGraph(dbg); //Split graph into components
        subdatasets = oneline::initialize<Subdataset>(std::move(components));//Create subdatasets corresponding to components
    } else {
        logger.info() << "Extracting subdatasets around contigs" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        size_t radius = std::stoull(parameterValues.getValue("radius"));
        for(StringContig scontig : io::SeqReader(paths_lib)) {
            Contig contig = scontig.makeContig();
            std::cout << contig.getInnerId() << " " << contig.truncSize() << " " << index.carefulAlign(contig).size() << std::endl;
            storage.addContig(contig);
            std::vector<ag::AlignmentChain<Contig, dbg::Edge>> contig_al = index.carefulAlign(contig);
            subdatasets.emplace_back(ag::Component::neighbourhood(dbg, contig_al, k + radius));
            subdatasets.back().id = contig.getInnerId();
        }
    }
    for(StringContig stringContig: io::SeqReader(contig_lib)) {
        storage.addContig(stringContig.makeContig());
    }
    logger.info() << "Filling path storage" << std::endl;
    storage.Fill(threads, index);
    FillSubdatasets(subdatasets, {&readStorage}, true);//Assign reads to datasets
    size_t cnt = 0;
    ag::Printer printer(ag::EdgePrintStyles::defaultDotInfo() + storage.edgeInfo() + ag::EdgeInfo::Tooltiper(readStorage.getSuffixes().labeler()));
    printer.printDot(dir / "graph.dot", ag::Component(dbg));
    //printDot(dir / "graph.dot", ag::Component(dbg), storage.labeler() + readStorage.getSuffixes().labeler());
    for(const Subdataset &subdataset: subdatasets) {//Print subdatasets to disk
        logger.info() << "Printing subdataset " << cnt << " " << subdataset.id << ":";
        for(dbg::Vertex &v : subdataset.component.verticesUnique()) {
            logger << " " << v.getId();
        }
        logger << "\n";
        std::string name = itos(cnt);
        if(!subdataset.id.empty())
            name += "_" + name;
        subdataset.Save(subdir / name, printer);
        cnt++;
    }
    logger.info() << "Finished extracting subdatasets" << std::endl;
    return 0;
}
