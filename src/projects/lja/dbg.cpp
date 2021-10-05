//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include "error_correction/diploidy_analysis.hpp"
#include "error_correction/mult_correction.hpp"
#include "error_correction/parameter_estimator.hpp"
#include "dbg/visualization.hpp"
#include "error_correction/error_correction.hpp"
#include "dbg/graph_algorithms.hpp"
#include "dbg/dbg_construction.hpp"
#include "dbg/dbg_disjointigs.hpp"
#include "dbg/minimizer_selection.hpp"
#include "dbg/sparse_dbg.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include "common/hash_utils.hpp"
#include "error_correction/initial_correction.hpp"
#include "sequences/seqio.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include "../dbg/graph_printing.hpp"
#include "subdataset_processing.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <wait.h>
using namespace dbg;

void analyseGenome(SparseDBG &dbg, const std::string &ref_file, size_t min_len,
                   const std::experimental::filesystem::path &path_dump,
                   const std::experimental::filesystem::path &mult_dump, logging::Logger &logger) {
    logger.info() << "Reading reference" << std::endl;
    std::vector<StringContig> ref = io::SeqReader(ref_file).readAll();
    logger.info() << "Finished reading reference. Starting alignment" << std::endl;
    std::vector<Segment<Edge>> path;
    std::ofstream os;
    os.open(path_dump);
    size_t cur = 0;
    for(StringContig & contig : ref) {
        Sequence seq = contig.makeSequence();
        os << "New chromosome " << contig.id << "(" << contig.size() << ")" << std::endl;
        if(seq.size() < min_len) {
            continue;
        }
        auto tmp = GraphAligner(dbg).align(seq);
        for(size_t i = 0; i < tmp.size(); i++) {
            const Segment<Edge> &seg = tmp[i];
            os << "[" << cur << ", " << cur + seg.size() << "] -> [" << seg.left << ", " << seg.right <<"] ";
            os << tmp.getVertex(i).hash() << tmp.getVertex(i).isCanonical() << " "
               << tmp.getVertex(i + 1).hash() << tmp.getVertex(i + 1).isCanonical() << " "
               << seg.size() << " " << seg.contig().getCoverage() << std::endl;
            cur += seg.size();
        }
        logger.info() << "Aligned chromosome " << contig.id << " . Path length " << tmp.size() << std::endl;
        path.insert(path.end(), tmp.begin(), tmp.end());
    }
    os.close();
    logger.info() << "Reference path consists of " << path.size() << " edges" << std::endl;
    size_t max_cov = 50;
    std::vector<size_t> cov(max_cov + 1);
    std::vector<size_t> cov_len(max_cov + 1);
    std::vector<size_t> cov_bad(max_cov + 1);
    std::vector<size_t> cov_bad_len(max_cov + 1);
    std::vector<size_t> cov_good(max_cov + 1);
    std::vector<size_t> cov_good_len(max_cov + 1);
    std::unordered_map<Edge const *, size_t> eset;
    for(size_t i = 0; i < path.size(); i++) {
        const Edge &edge = path[i].contig();
        eset[&edge] += 1;
    }
    std::ofstream os_mult;
    os_mult.open(mult_dump);
    for(auto & it : eset) {
        os_mult << it.second << " " << it.first->getCoverage() << " " << it.first->size() << std::endl;
    }
    os_mult.close();
    for(auto & pair : dbg) {
        Vertex &vert = pair.second;
        for (Edge &edge : vert) {
            size_t cov_val = std::min(max_cov, size_t(edge.getCoverage()));
            if (eset.find(&edge) == eset.end() && eset.find(&edge.rc()) == eset.end()) {
                cov_bad[cov_val] += 1;
                cov_bad_len[cov_val] += edge.size();
            } else {
                cov_good[cov_val] += 1;
                cov_good_len[cov_val] += edge.size();
            }
            cov[cov_val] += 1;
            cov_len[cov_val] += edge.size();
        }
    }
    logger.info() << "All coverages" << std::endl;
    logger << cov << std::endl << cov_len << std::endl;
    logger.info() << "Coverages of edges in genome path" << std::endl;
    logger << cov_good << std::endl << cov_good_len << std::endl;
    logger.info() << "Coverages of edges outside genome path" << std::endl;
    logger << cov_bad << std::endl << cov_bad_len << std::endl;
}

void LoadCoverage(const std::experimental::filesystem::path &fname, logging::Logger &logger, SparseDBG &dbg) {
    logger.info() << "Loading edge coverages." << std::endl;
    std::ifstream is;
    is.open(fname);
    size_t n;
    is >> n;
    for (size_t i = 0; i < n; i++) {
        hashing::htype vid;
        is >> vid;
        Vertex *v = &dbg.getVertex(vid);
        size_t inDeg, outDeg;
        is >> outDeg >> inDeg;
        for (size_t j = 0; j < inDeg + outDeg; j++) {
            if (j == outDeg)
                v = &v->rc();
            size_t next;
            is >> next;
            Edge &edge = v->getOutgoing(char(next));
            size_t cov;
            is >> cov;
            edge.incCov(cov);
        }
    }
    is.close();
    logger.info() << "Finished loading edge coverages." << std::endl;
}

std::string constructMessage() {
    std::stringstream ss;
    ss << "JumboDB: a tool for constructing de Bruijn graph for arbitrarily large value of k\n";
    ss << "Usage: dbg [options] -o <output-dir> -k <int> --reads <reads_file> [--reads <reads_file2> ...]\n\n";
    ss << "Basic options:\n";
    ss << "  -o <file_name> (or --output-dir <file_name>)  Name of output folder. Resulting graph will be stored there.\n";
    ss << "  -k <int>                                      Value of k (vertex size) to be used for de Bruijn graph construction. k should be odd (otherwise k + 1 is used instead).\n";
    ss << "  --reads <file_name>                           Name of file that contains reads in fasta or fastq format. This option can be used any number of times in the same command line. In this case reads from all specified files will be used as an input.\n";
    ss << "  -h (or --help)                                Print this help message.\n";
    ss << "\nAdvanced options:\n";
    ss << "  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.\n";
    ss << "  -w <int> (or --window <int>`)                 The window size to be used for sparse de Bruijn graph construction. The default value is 2000. Note that all reads of length less than k + w are ignored during graph construction.\n";
    ss << "  --compress                                    Compress all homolopymers in reads.\n";
    ss << "  --coverage                                    Calculate edge coverage of edges in the constructed de Bruijn graph.\n";
    return ss.str();
}

int main(int argc, char **argv) {
    CLParser parser({"vertices=none", "unique=none", "coverages=none", "dbg=none", "output-dir=",
                     "threads=16", "k-mer-size=", "window=2000", "base=239", "debug", "disjointigs=none", "reference=none",
                     "correct", "simplify", "coverage", "cov-threshold=2", "rel-threshold=10", "tip-correct",
                     "initial-correct", "mult-correct", "mult-analyse", "compress", "dimer-compress=1000000000,1000000000,1", "help", "genome-path",
                     "dump", "extension-size=none", "print-all", "extract-subdatasets", "print-alignments", "subdataset-radius=10000",
                     "split", "diploid"},
                    {"reads", "pseudo-reads", "align", "paths", "print-segment"},
                    {"h=help", "o=output-dir", "t=threads", "k=k-mer-size","w=window"},
                    constructMessage());
    parser.parseCL(argc, argv);
    if (parser.getCheck("help")) {
        std::cout << parser.message() << std::endl;
        return 0;
    }
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
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    const size_t w = std::stoi(parser.getValue("window"));
    logger << std::endl;
    logger.info() << "Hello. You are running jumboDBG, a tool for construction of de Bruijn graphs for arbitrarily large values of k\n";
    logger.info() << "Note that jumboDB does not perform any error correction and ignores all reads shorter than k + w = " << k + w << std::endl;
    if(parser.getCheck("extract-subdatasets")) {
        logger.info() << "Enabled subdataset extraction" << std::endl;
    }
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    hashing::RollingHash hasher(k, std::stoi(parser.getValue("base")));
    io::Library pseudo_reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("pseudo-reads"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library paths_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));
    io::Library genome_lib = {};
    if (parser.getValue("reference") != "none") {
        logger.info() << "Added reference to graph construction. Careful, some edges may have coverage 0" << std::endl;
        genome_lib = {std::experimental::filesystem::path(parser.getValue("reference"))};
    }
    io::Library construction_lib = reads_lib + pseudo_reads_lib + genome_lib;
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);

    std::string disjointigs_file = parser.getValue("disjointigs");
    std::string vertices_file = parser.getValue("vertices");
    std::string dbg_file = parser.getValue("dbg");
    SparseDBG dbg = dbg_file == "none" ?
                    DBGPipeline(logger, hasher, w, construction_lib, dir, threads, disjointigs_file, vertices_file) :
                    LoadDBGFromFasta({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);

    bool calculate_alignments = parser.getCheck("initial-correct") ||
            parser.getCheck("mult-correct") || parser.getCheck("print-alignments") || parser.getCheck("split");
    bool calculate_coverage = parser.getCheck("coverage") || parser.getCheck("simplify") ||
            parser.getCheck("correct") || parser.getValue("reference") != "none" ||
            parser.getCheck("tip-correct") ||
            parser.getCheck("initial-correct") || parser.getCheck("mult-correct") || !paths_lib.empty();
    calculate_coverage = calculate_coverage && !calculate_alignments;
    if (!parser.getListValue("align").empty() || parser.getCheck("print-alignments") ||
                parser.getCheck("mult-correct") || parser.getCheck("mult-analyse") ||
                parser.getCheck("split")|| calculate_coverage || calculate_alignments) {
        dbg.fillAnchors(w, logger, threads);
    }

    if (calculate_coverage) {
        if (parser.getValue("coverages") == "none") {
            CalculateCoverage(dir, hasher, w, reads_lib, threads, logger, dbg);
        } else {
            LoadCoverage(parser.getValue("coverages"), logger, dbg);
        }
    }
    size_t extension_size = std::max<size_t>(k * 5 / 2, 3000);
    if (k > 1500)
        extension_size = 1000000;
    if(parser.getValue("extension-size") != "none")
        extension_size = std::stoull(parser.getValue("extension-size"));

    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true);
    RecordStorage refStorage(dbg, 0, extension_size, threads, readLogger, false, false);

    if(calculate_alignments) {
        logger.info() << "Collecting read alignments" << std::endl;
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        logger.info() << "Collecting reference alignments" << std::endl;
        io::SeqReader refReader(genome_lib);
        refStorage.fill(refReader.begin(), refReader.end(), dbg, w + k - 1, logger, threads);
    }

    if(parser.getCheck("mult-analyse")) {
        NewMultCorrect(dbg, logger, dir, readStorage, 70000, threads, parser.getCheck("dump"));
    }

    if(parser.getCheck("mult-correct")) {
        MultCorrect(dbg, logger, dir, readStorage, 50000,threads, parser.getCheck("diploid"), debug);
    }

    if(parser.getCheck("initial-correct")) {
        size_t threshold = std::stoull(parser.getValue("cov-threshold"));
        size_t reliable = std::stoull(parser.getValue("rel-threshold"));
        std::vector<StringContig> ref_vector;
        if (parser.getValue("reference") != "none") {
            ref_vector = io::SeqReader(parser.getValue("reference")).readAll();
        }
        initialCorrect(dbg, logger, dir / "correction.txt", readStorage, refStorage, threshold, 2 * threshold, reliable,
                       threads, parser.getCheck("dump"));
        Component comp(dbg);
        DrawSplit(comp, dir / "split");
    }

    if(!paths_lib.empty()) {
        logger << "Printing additional figures with paths" << std::endl;
        logger << "Aligning contigs from paths files" << std::endl;
        GraphAlignmentStorage storage(dbg);
        io::SeqReader reader(paths_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.fill(contig);
        }
        {
            logger.info() << "Printing graph with paths to dot file " << (dir / "paths.dot") << std::endl;
            std::ofstream coordinates_dot;
            coordinates_dot.open(dir / "paths.dot");
            printDot(coordinates_dot, Component(dbg), storage.labeler());
            coordinates_dot.close();
        }
        {
            std::ofstream coordinates;
            coordinates.open(dir / "paths.txt");
            storage.print(coordinates);
            coordinates.close();
        }
    }

    if(parser.getCheck("print-alignments") || parser.getCheck("mult-correct")) {
        readStorage.printAlignments(logger, dir/"alignments.txt");
    }

    if(parser.getCheck("mult-correct") || parser.getCheck("initial-correct")) {
        readStorage.printFasta(logger, dir/"corrected.fasta");
    }

    if(parser.getCheck("print-all")) {
        logger.info() << "Printing segments of paths" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        GraphAlignmentStorage storage(dbg);
        io::SeqReader reader(paths_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.fill(contig);
        }
        reader.reset();
        size_t cnt = 0;
        for(StringContig scontig : reader) {
            cnt += 1;
            if(cnt > 200)
                break;
            Contig contig = scontig.makeContig();
            std::string fname = "contig_" + std::to_string(cnt) + "_" + mask(contig.id) + ".dot";
            const std::experimental::filesystem::path seg_file = dir / fname;
            logger.info() << "Printing contig " << contig.id << " to dot file " << (seg_file) << std::endl;
            std::ofstream coordinates_dot;
            Component comp = Component::neighbourhood(dbg, contig, dbg.hasher().getK() + 100);
            coordinates_dot.open(seg_file);
            printDot(coordinates_dot, Component(comp), storage.labeler());
            coordinates_dot.close();
        }
    }

    if(parser.getCheck("split")) {
        std::experimental::filesystem::path subdatasets_dir = dir / "subdatasets";
        ensure_dir_existance(subdatasets_dir);
        std::experimental::filesystem::path executable(argv[0]);
        std::experimental::filesystem::path py_path = executable.parent_path() / "run_rr.py";
        logger.trace() << "py_path set to " << py_path.string() << std::endl;
        RepeatResolver rr(dbg, {&readStorage}, subdatasets_dir, py_path, true);
        logger.info() << "Extracting subdatasets for connected components" << std::endl;
        std::function<bool(const Edge &)> is_unique = [](const Edge &){return false;};
        std::vector<RepeatResolver::Subdataset> subdatasets = rr.SplitDataset(is_unique);
#pragma omp parallel for schedule(dynamic, 1) default(none) shared(subdatasets, logger, rr)
        for(size_t snum = 0; snum < subdatasets.size(); snum++) {
            RepeatResolver::Subdataset &subdataset = subdatasets[snum];
            std::experimental::filesystem::path outdir = subdataset.dir / "mltik";
            ensure_dir_existance(outdir);
            rr.prepareDataset(subdataset);
        }
        logger.info() << "Finished extracting subdatasets for connected components" << std::endl;
    }

    if(parser.getCheck("extract-subdatasets")) {
        std::experimental::filesystem::path subdatasets_dir = dir / "subdatasets";
        ensure_dir_existance(subdatasets_dir);
        logger.info() << "Extracting subdatasets around contigs" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        std::vector<Component> comps;
        std::vector<std::ofstream *> os;
        GraphAlignmentStorage storage(dbg);
        io::SeqReader reader(paths_lib);
        size_t cnt = 0;
        size_t radius = std::stoull(parser.getValue("subdataset-radius"));
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            comps.emplace_back(Component::neighbourhood(dbg, contig, dbg.hasher().getK() + radius));
            os.emplace_back(new std::ofstream());
            os.back()->open(subdatasets_dir / (std::to_string(cnt) + ".fasta"));
            cnt += 1;
        }
        io::SeqReader read_reader(reads_lib);
        std::function<void(size_t, StringContig &)> task = [&dbg, &comps, &os, &hasher, w, &logger](size_t pos, StringContig & scontig) {
            string initial_seq = scontig.seq;
            Contig contig = scontig.makeContig();
            if(contig.size() < hasher.getK() + w - 1)
                return;
            GraphAlignment al = GraphAligner(dbg).align(contig.seq);
            for(size_t j = 0; j < comps.size(); j++) {
                for(size_t i = 0; i <= al.size(); i++) {
                    if(comps[j].contains(al.getVertex(i))) {
#pragma omp critical
                        *os[j] << ">" << contig.id << "\n" <<initial_seq << "\n";
                        break;
                    }
                }
            }
        };
        processRecords(read_reader.begin(), read_reader.end(), logger, threads, task);
        for(size_t j = 0; j < comps.size(); j++) {
            os[j]->close();
            delete os[j];
            os[j] = nullptr;
        }
    }

    std::vector<std::string> path_segments = parser.getListValue("print-segment");
    if(!path_segments.empty()) {
        logger.info() << "Printing segments of paths" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        GraphAlignmentStorage storage(dbg);
        std::vector<std::tuple<std::string, size_t, size_t, std::string>> seg_recs;
        for(const std::string& s : path_segments) {
            std::vector<std::string> parsed = split(s, "[,]");
            seg_recs.emplace_back(parsed[0], std::stoull(parsed[1]), std::stoull(parsed[2]), s);
        }
        std::vector<Contig> segs;
        io::SeqReader reader(paths_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.fill(contig);
            for(auto & seg_rec : seg_recs) {
                if(std::get<0>(seg_rec) == contig.id) {
                    segs.emplace_back(contig.seq.Subseq(std::get<1>(seg_rec), std::min(contig.size(), std::get<2>(seg_rec))), std::get<3>(seg_rec));
                } else if (std::get<0>(seg_rec) == "-" + contig.id) {
                    segs.emplace_back((!contig.seq).Subseq(std::get<1>(seg_rec), std::min(contig.size(), std::get<2>(seg_rec))), std::get<3>(seg_rec));
                }
            }
        }
        for(Contig &seg : segs) {
            const std::experimental::filesystem::path seg_file = dir / ("seg_" + mask(seg.id) + ".dot");
            logger.info() << "Printing segment " << seg.id << " to dot file " << (seg_file) << std::endl;
            std::ofstream coordinates_dot;
            Component comp = Component::neighbourhood(dbg, seg, dbg.hasher().getK() + 100);
            coordinates_dot.open(seg_file);
            printDot(coordinates_dot, Component(dbg), storage.labeler());
            coordinates_dot.close();
        }
    }

    if(parser.getCheck("genome-path")) {
        logger.info() << "Printing additional figures with reference alignments" << std::endl;
        logger.info() << "Aligning genome" << std::endl;
        GraphAlignmentStorage storage(dbg);
        io::SeqReader reader(genome_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.fill(contig);
        }
        {
            logger.info() << "Printing graph to dot file " << (dir / "genome_path.dot") << std::endl;
            std::ofstream coordinates_dot;
            coordinates_dot.open(dir / "genome_path.dot");
            printDot(coordinates_dot, Component(dbg), storage.labeler());
            coordinates_dot.close();
        }
        {
            std::ofstream coordinates;
            coordinates.open(dir / "genome_path.txt");
            storage.print(coordinates);
            coordinates.close();
        }
    }

    if(parser.getValue("dbg") == "none") {
        if(debug) {
            logger.info() << "Printing graph to fasta file " << (dir / "graph.fasta") << std::endl;
            printFasta(dir / "graph.fasta", Component(dbg));
            logger.info() << "Printing assembly to fasta file " << (dir / "assembly.fasta") << std::endl;
            printAssembly(dir / "assembly.fasta", Component(dbg));
        } else {
            logger.info() << "Printing graph to fasta file " << (dir / "graph.fasta") << std::endl;
            printAssembly(dir / "graph.fasta", Component(dbg));
        }
        logger.info() << "Printing graph to gfa file " << (dir / "graph.gfa") << std::endl;
        printGFA(dir / "graph.gfa", Component(dbg), calculate_coverage);
        logger.info() << "Printing graph to dot file " << (dir / "graph.dot") << std::endl;
        printDot(dir / "graph.dot", Component(dbg));
    }

    if (parser.getCheck("tip-correct")) {
        logger.info() << "Removing tips from reads" << std::endl;
        std::experimental::filesystem::path out = dir / "tip_correct.fasta";
        ParallelRecordCollector<Contig> alignment_results(threads);

        std::function<void(size_t, StringContig &)> task = [&dbg, &alignment_results, &hasher, w, &logger](size_t pos, StringContig & contig) {
            Contig read = contig.makeContig();
            if(read.size() < w + hasher.getK() - 1)
                return;
            GraphAlignment gal = GraphAligner(dbg).align(read.seq);
            if (gal.size() > 0 && gal.front().contig().getCoverage() < 2 && gal.start().inDeg() == 0 && gal.start().outDeg() == 1) {
                gal = gal.subalignment(1, gal.size());
            }
            if (gal.size() > 0 && gal.back().contig().getCoverage() < 2 && gal.finish().outDeg() == 0 && gal.finish().inDeg() == 1) {
                gal = gal.subalignment(0, gal.size() - 1);
            }
            for(Segment<Edge> seg : gal) {
                if (seg.contig().getCoverage() < 2)
                    return;
            }
            if (gal.size() > 0) {
                alignment_results.emplace_back(gal.Seq(), read.id);
            }
        };
        std::ofstream os(dir / "tip_correct.fasta");
        io::SeqReader reader(reads_lib);
        processRecords(reader.begin(), reader.end(), logger, threads, task);
        for(Contig & rec : alignment_results) {
            os << ">" << rec.id << "\n" << rec.seq << "\n";
        }
        os.close();
    }

    if (!parser.getListValue("align").empty()) {
        io::Library align_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("align"));
        alignLib(logger, dbg, align_lib, hasher, w, dir, threads);
    }

//    findTips(logger, dbg, threads);
    if(parser.getCheck("correct")) {
        io::SeqReader reader(reads_lib);
        error_correction::correctSequences(dbg, logger, reader.begin(), reader.end(),
                                           dir / "corrected.fasta", dir / "bad.fasta", threads, w + hasher.getK() - 1);
    }
    if (parser.getValue("reference") != "none") {
        analyseGenome(dbg, parser.getValue("reference"), k + w - 1, dir / "ref.info", dir / "mult.info", logger);
    }
    if (parser.getCheck("simplify")) {
        logger.info() << "Removing low covered edges" << std::endl;
        size_t threshold = std::stoull(parser.getValue("cov-threshold"));
        std::vector<Sequence> edges;
        std::vector<hashing::htype> vertices_again;
        for(auto & it : dbg) {
            Vertex &vert = it.second;
            bool add = false;
            for(Edge & edge : vert) {
                if (edge.getCoverage() >= threshold) {
                    edges.push_back(vert.seq + edge.seq);
                    add = true;
                }
            }
            for(Edge & edge : vert.rc()) {
                if (edge.getCoverage() >= threshold){
                    edges.push_back(vert.rc().seq + edge.seq);
                    add = true;
                }
            }
            if (add)
                vertices_again.push_back(vert.hash());
        }
        SparseDBG simp_dbg(vertices_again.begin(), vertices_again.end(), hasher);
        FillSparseDBGEdges(simp_dbg, edges.begin(), edges.end(), logger, threads, 0);
        for(auto & it : simp_dbg) {
            Vertex &vert = it.second;
            Vertex &other = dbg.getVertex(vert.hash());
            for(Edge & edge : vert) {
                edge.incCov(other.getOutgoing(edge.seq[0]).intCov());
            }
            for(Edge & edge : vert.rc()) {
                edge.incCov(other.rc().getOutgoing(edge.seq[0]).intCov());
            }
        }
        printDot(dir / "simp_graph1.dot", Component(simp_dbg));
        mergeAll(logger, simp_dbg, threads);
        printFasta(dir / "simp_graph.fasta", Component(simp_dbg));
        printDot(dir / "simp_graph.dot", Component(simp_dbg));
    }
    logger.info() << "DBG construction finished" << std::endl;
    logger.info() << "Please cite our paper if you use jumboDBG in your research: https://www.biorxiv.org/content/10.1101/2020.12.10.420448" << std::endl;
    return 0;
}
