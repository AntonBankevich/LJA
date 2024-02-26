//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include "error_correction/diploidy_analysis.hpp"
#include "error_correction/mult_correction.hpp"
#include "error_correction/parameter_estimator.hpp"
#include "dbg/visualization.hpp"
#include "dbg/graph_algorithms.hpp"
#include "dbg/dbg_construction.hpp"
#include "dbg/dbg_disjointigs.hpp"
#include "dbg/minimizer_selection.hpp"
#include "dbg/sparse_dbg.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include "common/hash_utils.hpp"
#include "error_correction/tournament_correction.hpp"
#include "sequences/seqio.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include "../dbg/graph_printing.hpp"
#include "subdataset_processing.hpp"
#include "dbg/dbg_graph_aligner.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <wait.h>
#include <dbg/id_index.hpp>

using namespace dbg;

void analyseGenome(SparseDBG &dbg, KmerIndex &index, const std::string &ref_file,
                   const std::experimental::filesystem::path &path_dump,
                   const std::experimental::filesystem::path &cov_dump,
                   const std::experimental::filesystem::path &mult_dump, logging::Logger &logger) {
    logger.info() << "Reading reference" << std::endl;
    std::vector<StringContig> ref = io::SeqReader(ref_file).readAll();
    logger.info() << "Finished reading reference. Starting alignment" << std::endl;
    std::vector<dbg::GraphPath> paths;
    std::ofstream os;
    os.open(path_dump);
    size_t cur = 0;
    std::unordered_map<Edge *, size_t> mult;
    size_t num = 0;
    for(StringContig & contig : ref) {
        Sequence seq = contig.makeSequence();
        os << "New chromosome " << contig.id << "(" << contig.size() << ")" << std::endl;
        if(seq.size() < index.minReadLen()) {
            continue;
        }
        auto tmp = index.align(seq);
        for(size_t i = 0; i < tmp.size(); i++) {
            const Segment<Edge> &seg = tmp[i];
            mult[&seg.contig()]++;
            mult[&seg.contig().rc()]++;
            os << "[" << cur << ", " << cur + seg.size() << "] -> " << tmp[i].contig().getInnerId() << " [" << seg.left << ", " << seg.right << "]\n";
            cur += seg.size();
        }
        logger.info() << "Aligned chromosome " << contig.id << " . Path length " << tmp.size() << std::endl;
        num += tmp.size();
        paths.emplace_back(std::move(tmp));
    }
    os.close();
    std::ofstream mos;
    mos.open(mult_dump);
    for(Edge &edge: dbg.edges()) {
        mos << edge.getInnerId() << " " << mult[&edge] << "\n";
    }
    mos.close();
    logger.info() << "Reference path consists of " << num << " edges" << std::endl;
    size_t max_cov = 50;
    std::vector<size_t> cov(max_cov + 1);
    std::vector<size_t> cov_len(max_cov + 1);
    std::vector<size_t> cov_bad(max_cov + 1);
    std::vector<size_t> cov_bad_len(max_cov + 1);
    std::vector<size_t> cov_good(max_cov + 1);
    std::vector<size_t> cov_good_len(max_cov + 1);
    std::unordered_map<Edge const *, size_t> eset;
    for(dbg::GraphPath &path: paths)
        for(Edge &edge : path.edges())
            eset[&edge] += 1;
    std::ofstream os_mult;
    os_mult.open(cov_dump);
    for(auto & it : eset) {
        os_mult << it.second << " " << it.first->getCoverage() << " " << it.first->truncSize() << std::endl;
    }
    os_mult.close();
    for(auto & vert : dbg.verticesUnique()) {
        for (Edge &edge : vert) {
            size_t cov_val = std::min(max_cov, size_t(edge.getCoverage()));
            if (eset.find(&edge) == eset.end() && eset.find(&edge.rc()) == eset.end()) {
                cov_bad[cov_val] += 1;
                cov_bad_len[cov_val] += edge.truncSize();
            } else {
                cov_good[cov_val] += 1;
                cov_good_len[cov_val] += edge.truncSize();
            }
            cov[cov_val] += 1;
            cov_len[cov_val] += edge.truncSize();
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
    IdIndex<Vertex> id_index(dbg.vertices().begin(), dbg.vertices().end());
    for (size_t i = 0; i < n; i++) {
        Vertex::id_type vid;
        is >> vid;
        Vertex *v = &id_index.getById(vid);
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
    ss << "JumboDBG: a tool for constructing de Bruijn graph for arbitrarily large value of k\n";
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
    AlgorithmParameters parameters({"vertices=none", "unique=none", "coverages=none", "dbg=none", "output-dir=",
                                   "threads=16", "k-mer-size=", "window=2000", "base=239", "debug", "disjointigs=none", "reference=none",
                                   "simplify", "coverage", "cov-threshold=2", "rel-threshold=10", "tip-correct",
                                   "initial-correct", "mult-correct", "mult-analyse", "compress", "dimer-compress=1000000000,1000000000,1", "help", "genome-path",
                                   "dump", "extension-size=none", "print-all", "extract-subdatasets", "print-alignments", "subdataset-radius=10000",
    "split", "diploid"}, {"reads", "pseudo-reads", "align", "paths", "print-segment"}, constructMessage());
    CLParser parser(parameters, {"h=help", "o=output-dir", "t=threads", "k=k-mer-size","w=window"});
    AlgorithmParameterValues params = parser.parseCL(argc, argv);
    if (params.getCheck("help")) {
        std::cout << params.helpMessage() << std::endl;
        return 0;
    }
    if (!params.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << params.checkMissingValues() << "\n" << std::endl;
        std::cout << params.helpMessage() << std::endl;
        return 1;
    }

    bool debug = params.getCheck("debug");
    StringContig::homopolymer_compressing = params.getCheck("compress");
    StringContig::SetDimerParameters(params.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(params.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    size_t k = std::stoi(params.getValue("k-mer-size"));
    const size_t w = std::stoi(params.getValue("window"));
    logger << std::endl;
    logger.info() << "Hello! You are running jumboDBG, a tool for construction of de Bruijn graphs for arbitrarily large values of k\n";
    logger.info() << "Note that jumboDBG does not perform any error correction and ignores all reads shorter than k + w = " << k + w << std::endl;
    if(params.getCheck("extract-subdatasets")) {
        logger.info() << "Enabled subdataset extraction" << std::endl;
    }
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    hashing::RollingHash hasher(k);
    io::Library pseudo_reads_lib = oneline::initialize<std::experimental::filesystem::path>(params.getListValue("pseudo-reads"));
    io::Library reads_lib = oneline::initialize<std::experimental::filesystem::path>(params.getListValue("reads"));
    io::Library paths_lib = oneline::initialize<std::experimental::filesystem::path>(params.getListValue("paths"));
    io::Library genome_lib = {};
    if (params.getValue("reference") != "none") {
        logger.info() << "Added reference to graph construction. Careful, some edges may have coverage 0" << std::endl;
        genome_lib = {std::experimental::filesystem::path(params.getValue("reference"))};
    }
    io::Library construction_lib = reads_lib + pseudo_reads_lib + genome_lib;
    size_t threads = std::stoi(params.getValue("threads"));
    omp_set_num_threads(threads);

    std::string disjointigs_file = params.getValue("disjointigs");
    std::string vertices_file = params.getValue("vertices");
    std::string dbg_file = params.getValue("dbg");
    SparseDBG dbg = dbg_file == "none" ?
                    DBGPipeline(logger, hasher, w, construction_lib, dir, threads, disjointigs_file, vertices_file) :
                    LoadDBGFromEdgeSequences({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);

    bool calculate_alignments = params.getCheck("initial-correct") ||
            params.getCheck("mult-correct") || params.getCheck("print-alignments") || params.getCheck("split");
    bool calculate_coverage = params.getCheck("coverage") || params.getCheck("simplify") ||
            params.getValue("reference") != "none" || params.getCheck("tip-correct") ||
            params.getCheck("initial-correct") || params.getCheck("mult-correct") || !paths_lib.empty();
    calculate_coverage = calculate_coverage && !calculate_alignments;
    KmerIndex index(dbg);
    if (!params.getListValue("align").empty() || params.getCheck("print-alignments") ||
                params.getCheck("mult-correct") || params.getCheck("mult-analyse") ||
                params.getCheck("split")|| calculate_coverage || calculate_alignments) {
        index.fillAnchors(logger, threads, dbg, w);
    }

    if (calculate_coverage) {
        if (params.getValue("coverages") == "none") {
            CalculateCoverage(logger, threads, dbg, index, dir, reads_lib);
        } else {
            LoadCoverage(params.getValue("coverages"), logger, dbg);
        }
    }
    size_t extension_size = std::max<size_t>(k * 5 / 2, 3000);
    if (k > 1500)
        extension_size = 1000000;
    if(params.getValue("extension-size") != "none")
        extension_size = std::stoull(params.getValue("extension-size"));

    ReadLogger readLogger(threads, dir/"read_log.txt");
    RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true);
    RecordStorage refStorage(dbg, 0, extension_size, threads, readLogger, false, false);

    if(calculate_alignments) {
        logger.info() << "Collecting read alignments" << std::endl;
        io::SeqReader reader(reads_lib);
        readStorage.fill(logger, threads, reader.begin(), reader.end(), dbg, index);
        logger.info() << "Collecting reference alignments" << std::endl;
        io::SeqReader refReader(genome_lib);
        refStorage.fill(logger, threads, refReader.begin(), refReader.end(), dbg, index);
    }

    if(params.getCheck("mult-correct")) {
        MultCorrect(logger, threads, dbg, dir, readStorage, 50000, 0, params.getCheck("diploid"), debug);
    }

    if(params.getCheck("initial-correct")) {
        size_t threshold = std::stoull(params.getValue("cov-threshold"));
        size_t reliable = std::stoull(params.getValue("rel-threshold"));
        std::vector<StringContig> ref_vector;
        if (params.getValue("reference") != "none") {
            ref_vector = io::SeqReader(params.getValue("reference")).readAll();
        }
        initialCorrect(logger, threads, dbg, dir / "correction.txt", readStorage, refStorage,
                       threshold, 2 * threshold, reliable, false, 60000, params.getCheck("dump"));
        Component comp(dbg);
        DrawSplit(comp, dir / "split");
    }

    if(!paths_lib.empty()) {
        logger << "Printing additional figures with paths" << std::endl;
        logger << "Aligning contigs from paths files" << std::endl;
        GraphPathStorage storage(dbg);
        io::SeqReader reader(paths_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.addContig(contig);
        }
        storage.Fill(threads, index);
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

    if(params.getCheck("print-alignments") || params.getCheck("mult-correct")) {
        readStorage.printReadAlignments(logger, dir / "alignments.txt");
    }

    if(params.getCheck("mult-correct") || params.getCheck("initial-correct")) {
        readStorage.printReadFasta(logger, dir / "corrected.fasta");
    }

    if(params.getCheck("print-all")) {
        logger.info() << "Printing segments of paths" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        GraphPathStorage storage(dbg);
        io::SeqReader reader(paths_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.addContig(contig);
        }
        storage.Fill(threads, index);
        reader.reset();
        size_t cnt = 0;
        for(StringContig scontig : reader) {
            cnt += 1;
            if(cnt > 200)
                break;
            Contig contig = scontig.makeContig();
            std::string fname = "contig_" + std::to_string(cnt) + "_" + mask(contig.getInnerId()) + ".dot";
            const std::experimental::filesystem::path seg_file = dir / fname;
            logger.info() << "Printing contig " << contig.getInnerId() << " to dot file " << (seg_file) << std::endl;
            std::ofstream coordinates_dot;
            std::vector<PerfectAlignment<Contig, Edge>> contig_al = index.carefulAlign(contig);
            Component comp = Component::neighbourhood(dbg, contig_al, k + 100);
            coordinates_dot.open(seg_file);
            printDot(coordinates_dot, Component(comp), storage.labeler());
            coordinates_dot.close();
        }
    }

    if(params.getCheck("extract-subdatasets")) {
        std::experimental::filesystem::path subdatasets_dir = dir / "subdatasets";
        ensure_dir_existance(subdatasets_dir);
        logger.info() << "Extracting subdatasets around contigs" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        std::vector<Component> comps;
        std::vector<std::ofstream *> os;
        GraphPathStorage storage(dbg);
        io::SeqReader reader(paths_lib);
        size_t cnt = 0;
        size_t radius = std::stoull(params.getValue("subdataset-radius"));
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            std::vector<PerfectAlignment<Contig, Edge>> contig_al = index.carefulAlign(contig);
            comps.emplace_back(Component::neighbourhood(dbg, contig_al, k + radius));
            os.emplace_back(new std::ofstream());
            os.back()->open(subdatasets_dir / (std::to_string(cnt) + ".fasta"));
            cnt += 1;
        }
        io::SeqReader read_reader(reads_lib);
        std::function<void(size_t, StringContig &)> task = [&dbg, &comps, &os, &hasher, w, &index](size_t pos, StringContig & scontig) {
            string initial_seq = scontig.seq;
            Contig contig = scontig.makeContig();
            if(contig.truncSize() < hasher.getK() + w - 1)
                return;
            dbg::GraphPath al = index.align(contig.getSeq());
            for(size_t j = 0; j < comps.size(); j++) {
                for(size_t i = 0; i <= al.size(); i++) {
                    if(comps[j].contains(al.getVertex(i))) {
#pragma omp critical
                        *os[j] << ">" << contig.getInnerId() << "\n" << initial_seq << "\n";
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

    std::vector<std::string> path_segments = params.getListValue("print-segment");
    if(!path_segments.empty()) {
        logger.info() << "Printing segments of paths" << std::endl;
        logger.info() << "Aligning paths" << std::endl;
        GraphPathStorage storage(dbg);
        std::vector<std::tuple<std::string, size_t, size_t, std::string>> seg_recs;
        for(const std::string& s : path_segments) {
            std::vector<std::string> parsed = split(s, "[,]");
            seg_recs.emplace_back(parsed[0], std::stoull(parsed[1]), std::stoull(parsed[2]), s);
        }
        std::vector<Contig> segs;
        io::SeqReader reader(paths_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.addContig(contig);
            for(auto & seg_rec : seg_recs) {
                if(std::get<0>(seg_rec) == contig.getInnerId()) {
                    segs.emplace_back(contig.getSeq().Subseq(std::get<1>(seg_rec), std::min(contig.truncSize(), std::get<2>(seg_rec))), std::get<3>(seg_rec));
                } else if (std::get<0>(seg_rec) == "-" + contig.getInnerId()) {
                    segs.emplace_back((!contig.getSeq()).Subseq(std::get<1>(seg_rec), std::min(contig.truncSize(), std::get<2>(seg_rec))), std::get<3>(seg_rec));
                }
            }
        }
        storage.Fill(threads, index);
        for(Contig &seg : segs) {
            const std::experimental::filesystem::path seg_file = dir / ("seg_" + mask(seg.getInnerId()) + ".dot");
            logger.info() << "Printing segment " << seg.getInnerId() << " to dot file " << (seg_file) << std::endl;
            std::ofstream coordinates_dot;
            std::vector<PerfectAlignment<Contig, Edge>> contig_al = index.carefulAlign(seg);
            Component comp = Component::neighbourhood(dbg, contig_al, k + 100);
            coordinates_dot.open(seg_file);
            printDot(coordinates_dot, Component(dbg), storage.labeler());
            coordinates_dot.close();
        }
    }

    if(params.getCheck("genome-path")) {
        logger.info() << "Printing additional figures with reference alignments" << std::endl;
        logger.info() << "Aligning genome" << std::endl;
        GraphPathStorage storage(dbg);
        io::SeqReader reader(genome_lib);
        for(StringContig scontig : reader) {
            Contig contig = scontig.makeContig();
            storage.addContig(contig);
        }
        storage.Fill(threads, index);
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

    if(params.getValue("dbg") == "none") {
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

    if (params.getCheck("tip-correct")) {
        logger.info() << "Removing tips from reads" << std::endl;
        std::experimental::filesystem::path out = dir / "tip_correct.fasta";
        ParallelRecordCollector<Contig> alignment_results(threads);

        std::function<void(size_t, StringContig &)> task = [&dbg, &alignment_results, &hasher, w, &index](size_t pos, StringContig & contig) {
            Contig read = contig.makeContig();
            if(read.truncSize() < w + hasher.getK() - 1)
                return;
            dbg::GraphPath gal = index.align(read.getSeq());
            if (gal.size() > 0 && gal.front().contig().getCoverage() < 2 && gal.start().inDeg() == 0 && gal.start().outDeg() == 1) {
                gal = gal.subPath(1, gal.size());
            }
            if (gal.size() > 0 && gal.back().contig().getCoverage() < 2 && gal.finish().outDeg() == 0 && gal.finish().inDeg() == 1) {
                gal = gal.subPath(0, gal.size() - 1);
            }
            for(Segment<Edge> seg : gal) {
                if (seg.contig().getCoverage() < 2)
                    return;
            }
            if (gal.size() > 0) {
                alignment_results.emplace_back(gal.Seq(), read.getInnerId());
            }
        };
        std::ofstream os(dir / "tip_correct.fasta");
        io::SeqReader reader(reads_lib);
        processRecords(reader.begin(), reader.end(), logger, threads, task);
        for(Contig & rec : alignment_results) {
            os << ">" << rec.getInnerId() << "\n" << rec.getSeq() << "\n";
        }
        os.close();
    }

    if (!params.getListValue("align").empty()) {
        io::Library align_lib = oneline::initialize<std::experimental::filesystem::path>(params.getListValue("align"));
        alignLib(logger, threads, index, align_lib, dir);
    }

//    findTips(logger, dbg, threads);
    if (params.getValue("reference") != "none") {
        analyseGenome(dbg, index, params.getValue("reference"), dir / "ref.info", dir / "cov.info", dir / "mult.info", logger);
    }
    if (params.getCheck("simplify")) {
        logger.info() << "Removing low covered edges" << std::endl;
        size_t threshold = std::stoull(params.getValue("cov-threshold"));
        std::vector<Sequence> edges;
        std::vector<Sequence> vertices_again;
        for(auto & vert : dbg.verticesUnique()) {
            bool add = false;
            for(Edge & edge : vert) {
                if (edge.getCoverage() >= threshold) {
                    edges.push_back(vert.getSeq() + edge.truncSeq());
                    add = true;
                }
            }
            for(Edge & edge : vert.rc()) {
                if (edge.getCoverage() >= threshold){
                    edges.push_back(vert.rc().getSeq() + edge.truncSeq());
                    add = true;
                }
            }
            if (add)
                vertices_again.push_back(vert.getSeq());
        }
        SparseDBG simp_dbg(hasher);
        for(Sequence &seq : vertices_again) {
            simp_dbg.addKmerVertex(seq);
        }
        FillSparseDBGEdges(simp_dbg, edges.begin(), edges.end(), logger, threads, 0);
        for(auto & vert : simp_dbg.verticesUnique()) {
            Vertex &other = index.getVertex(hashing::KWH(hasher, vert.getSeq(), 0).hash());
            for(Edge & edge : vert) {
                edge.incCov(other.getOutgoing(edge.truncSeq()[0]).intCov());
            }
            for(Edge & edge : vert.rc()) {
                edge.incCov(other.rc().getOutgoing(edge.truncSeq()[0]).intCov());
            }
        }
        printDot(dir / "simp_graph1.dot", Component(simp_dbg));
        ag::MergeAll<DBGTraits>(logger, threads, simp_dbg);
        printFasta(dir / "simp_graph.fasta", Component(simp_dbg));
        printDot(dir / "simp_graph.dot", Component(simp_dbg));
    }
    logger.info() << "DBG construction finished" << std::endl;
    logger.info() << "Please cite our paper if you use jumboDBG in your research: https://www.biorxiv.org/content/10.1101/2020.12.10.420448" << std::endl;
    return 0;
}
