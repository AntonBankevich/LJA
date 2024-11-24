#pragma once

#include "diploidy_analysis.hpp"
#include "parameter_estimator.hpp"
#include "precorrection.hpp"
#include "dimer_correction.hpp"
#include <dbg/dbg_construction.hpp>
#include <dbg/graph_printing.hpp>
#include <dbg/graph_stats.hpp>
#include <dbg/graph_alignment_storage.hpp>
#include <dbg/visualization.hpp>
#include <sequences/seqio.hpp>

#include "assembly_graph/visualization.hpp"

void analyseGenome(SparseDBG &dbg, KmerIndex &index, const std::string &ref_file,
                   const std::experimental::filesystem::path &path_dump,
                   const std::experimental::filesystem::path &cov_dump,
                   const std::experimental::filesystem::path &mult_dump, logging::Logger &logger) {
    logger.info() << "Reading reference " << ref_file << std::endl;
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
        logger.info() << seq.size() << " : " << index.minReadLen() << "\n";
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

std::unordered_map<std::string, std::experimental::filesystem::path>
MLGraphEC(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                             const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
                                             size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
                                             bool diploid, bool debug, bool load, std::string reference) {

    logger.info() << "Performing coverage-based error correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k);
    io::Library genome_lib = {};
    if (reference != "none") {
        logger.info() << "Added reference to graph construction. Careful, some edges may have coverage 0" << std::endl;
        genome_lib = {std::experimental::filesystem::path(reference)};
    }
    io::Library construction_lib = reads_lib + pseudo_reads_lib + genome_lib;
    dbg::SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, construction_lib, dir, threads,
                                            (dir / "disjointigs.fasta").string(), (dir / "vertices.save").string())
                              :
                         DBGPipeline(logger, hasher, w, construction_lib, dir, threads);
    KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, w);
    CalculateCoverage(logger, threads, dbg, index, dir, reads_lib);
    size_t extension_size = std::max<size_t>(k * 2, 1000);
    ag::ReadLogger readLogger(threads, dir / "read_log.txt");
    dbg::ReadAlignmentStorage readStorage(dbg, 0, extension_size, true, true, false);
    readStorage.setReadLogger(readLogger);
    dbg::ReadAlignmentStorage refStorage(dbg, 0, extension_size, false, false);
    refStorage.setReadLogger(readLogger);
    io::SeqReader reader(reads_lib);
    readStorage.FillAlignments(logger, threads, reader.begin(), reader.end(), dbg, index);
    if(reference != "none") {
        io::SeqReader refReader(genome_lib);
        refStorage.FillAlignments(logger, threads, refReader.begin(), refReader.end(), dbg, index);
    }
    ObjInfo<dbg::Vertex> vertexInfo = VertexPrintStyles<dbg::DBGTraits>::defaultDotInfo();
    ObjInfo<dbg::Edge> edgeInfo = EdgePrintStyles<dbg::DBGTraits>::defaultDotInfo();
    Printer<DBGTraits> printer(vertexInfo, edgeInfo);
    printer.printDot(dir / "initial_dbg.dot", Component(dbg));
    coverageStats(logger, dbg);
    if (debug) {
        PrintPaths(logger, threads, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
    }
    Precorrector precorrector(4);
    DimerCorrector dimerCorrector(logger, dbg, readStorage, StringContig::max_dimer_size);
    BulgePathCorrector bpCorrector(dbg, readStorage, 80000, 1);
    ErrorCorrectionEngine(precorrector).run(logger, threads, dbg, readStorage);
    RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, extension_size);
    readStorage.trackSuffixes(logger, threads);
    ErrorCorrectionEngine(dimerCorrector).run(logger, threads, dbg, readStorage);
    RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, extension_size);
    DatasetParameters params = EstimateDatasetParameters(dbg, readStorage, true);
    params.Print(logger);
    logger.info() << "Saving to dot file\n";
    printer.printDot(dir / "graph.dot", Component(dbg));
    ObjInfo<dbg::Edge> edgeGFAInfo = EdgePrintStyles<dbg::DBGTraits>::defaultGFAInfo();
    Printer<DBGTraits> printer2(vertexInfo, edgeGFAInfo);
    printer2.printGFA(dir / "graph.gfa", Component(dbg), true);
    KmerIndex index2(dbg);
    index2.fillAnchors(logger, threads, dbg, w);
    if(reference != "none") {
        analyseGenome(dbg, index2, reference, dir / "ref.info", dir / "cov.info", dir / "mult.info", logger);
    }
    return {{"graph_dot", dir/ "graph.dot"}, {"graph_gfa", dir/ "graph.gfa"}, {"mult_info", dir / "mult.info"}, {"ref_info", dir / "ref.info" }};
}

class MLGraphCorrectionStage : public Stage {
public:
    MLGraphCorrectionStage() : Stage(AlgorithmParameters(
            {"k-mer-size=501", "window=2000", "reference=none", "coverage-threshold=3", "reliable-coverage=10", "diploid", "load", "dump-reads"},
            {}, ""), {"reads", "pseudo_reads", "paths"}, {"graph_dot", "graph_gfa", "mult_info", "ref_info"}) {
    }
protected:
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                          const std::experimental::filesystem::path &dir, bool debug,
                          const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override {
        size_t k = std::stoi(parameterValues.getValue("k-mer-size"));
        size_t w = std::stoi(parameterValues.getValue("window"));
        double reliable_coverage = std::stod(parameterValues.getValue("reliable-coverage"));
        double threshold = std::stod(parameterValues.getValue("coverage-threshold"));
        bool diploid = parameterValues.getCheck("diploid");
        bool load = parameterValues.getCheck("load");
        std::string reference = parameterValues.getValue("reference");
        return MLGraphEC(logger, dir, input.find("reads")->second, input.find("pseudo_reads")->second,
                          input.find("paths")->second, threads, k, w, threshold, reliable_coverage, diploid, debug, load, reference);
    }
};
