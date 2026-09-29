#include "multiplex_correction_stage.hpp"
#include "assembly_graph/data_structures/aligned_read_statistics_tracker.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "path_tracker.hpp"
#include <limits>

using namespace spg;

std::unordered_map<std::string, std::experimental::filesystem::path>
spg::RunMultiplexingAndCorrection(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                     size_t k, size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                     const io::Library &extra_reads_files, const io::Library &paths,
                     size_t initial_core_length, size_t max_core_length, size_t core_length_step,
                     double vertex_coverage_threshold, bool debug) {
    size_t unique_threshold = 40000;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    hashing::RollingHash hasher(k, 239);
    logger.info() << "Loading graph" << std::endl;
    dbg::SparseDBG spg = graph_gfa.empty() ? DBGPipeline(logger, hasher, w, reads_files, dir, threads) :
                         dbg::LoadDBGFromEdgeSequences(logger, threads, graph_gfa, hasher);
    logger.info() << "Loading reads" << std::endl;
    dbg::DBGAlignedReadStorage dbg_storage = dbg::DBGAlignedReadStorage::Load(logger, threads,
                                                                              reads_files + extra_reads_files, spg,
                                                                              true);
    dbg_storage.trackSuffixes(logger, threads, spg, 0, 10000000);
    dbg_storage.stopTrackCoverage();
//    SPGCoverage stays attached as a listener for the whole run: it keeps outgoing_read_count
//    correct through both the multiplexing (ResolutionListener) and correction (AlignedReadStorageListener)
//    steps below, which is what lets the two be interleaved.
    ag::AlignedReadStatisticsTracker SPGCoverage(logger, threads, spg, dbg_storage, dbg_storage.getSuffixes());
//    Separate listener computing per-vertex SPG coverage (see CoverageSamplingTracker); must be attached
//    before the edgeToSupreVertex loop below so it can seed samples from the original DBG edges.
    ag::CoverageSamplingTracker spgVertexCoverage(spg, dbg_storage, k, threads);
    spg.disableHashing();
    ag::EdgeInfo edge_info = ag::EdgePrintStyles::defaultDotInfo() + ag::EdgeInfo::Labeler(SPGCoverage.getEdgeLabeler());
    ag::DLLAlignmentStorage path_storage(spg);
    ag::Printer printer(ag::VertexPrintStyles::spgLabeler() + ag::VertexInfo::Labeler(SPGCoverage.getVertexLabeler()) +
                         ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()), edge_info);
    if (debug) {
        PrepareDLLPathTracker(logger, threads, spg, w, paths, path_storage);
        path_storage.startLogging(logger.trace());
        path_storage.startDrawing(dir / "path_tracking", printer);
    } else {
        path_storage.detach();
    }
    std::experimental::filesystem::path figs = dir / "figs";
    if (debug) {
        recreate_dir(figs);
        ag::Printer dbg_printer(ag::VertexPrintStyles::defaultDotInfo() + ag::VertexInfo::Labeler(SPGCoverage.getVertexLabeler()) +
                         ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()), edge_info);
        dbg_printer.printDot(figs / "dbg_initial.dot", spg);
    }
    std::vector<ag::EdgeId> eids = oneline::map(spg.edgesUnique().begin(), spg.edgesUnique().end(), IdTransformer<Edge>());
    UniqueVertexStorage unique_storage(spg, 40000);
    //TODO: create convert function that would return a new supregraph package an dthis runs in parallel.
    for (ag::EdgeId eid : eids) {
        Vertex &new_vertex = spg.edgeToSupreVertex(*eid);
    }

    ag::VertexCoverageReliableFiller vertexReliableFiller(vertex_coverage_threshold);
    ag::ReliablePathCorrector reliablePathCorrector;

    logger.info() << "Multiplexing and correcting" << std::endl;
    AndreyRule rule(dbg_storage.getSuffixes(), unique_storage);
    std::vector<size_t> thresholds;
    for (size_t threshold = initial_core_length; threshold < max_core_length; threshold += core_length_step)
        thresholds.push_back(threshold);
    thresholds.push_back(max_core_length);
    thresholds.push_back(std::numeric_limits<size_t>::max());
    spg::Multiplexer multiplexer(spg, dbg_storage, rule, thresholds.front());
    if (debug) {
        printer.setEdgeInfo(ag::EdgePrintStyles::spgLabeler() + ag::EdgeInfo::Labeler(SPGCoverage.getEdgeLabeler()) +
                             ag::EdgePrintStyles::simpleColorer("black"));
        printer.setVertexInfo(ag::VertexPrintStyles::spgLabeler() + ag::VertexInfo::Labeler(SPGCoverage.getVertexLabeler()) +
                               ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()) +
                               ag::VertexPrintStyles::defaultDotColorer() + ag::VertexPrintStyles::defaultTooltiper() +
                               ag::VertexInfo::Colorer(unique_storage.getColorer("white", "green")));
        printer.printDot(dir / "initial.dot", spg);
    }
    size_t cnt = 1;
    if (debug)
        dbg_storage.logReads(threads, dir/"read_log.txt");
    for (size_t threshold : thresholds) {
        logger.info() << "Multiplexing with core-length threshold "
                       << (threshold == std::numeric_limits<size_t>::max() ? std::string("inf") : itos(threshold)) << std::endl;
        multiplexer.setMaxCoreLength(threshold);
        size_t mult_cnt = 0;
        while (multiplexer.hasReadyCore() || multiplexer.hasPendingMerge()) {
            auto res = multiplexer.process(logger, threads);
            if (debug && !res.empty()) {
                logger.trace() << "Operation " << cnt << ": " << res << std::endl;
                // printer.printDot(figs / ("supregraph_" + itos(cnt, 5) + ".dot"), ag::Component::neighbourhood(spg, res, 100000, 20));
                cnt++;
            }
            if (res.size() > 1)
                mult_cnt++;
        }
        // for (Vertex & v : spg.vertices()) {
        //     if (v.isCore() && (v.inDeg() == 1 || v.outDeg() == 1)) {
        //         VERIFY(v.coverage_info.weight == 0);
        //     }
        // }
        logger.info() << "Multiplexing performed " << mult_cnt << " times" << std::endl;
        spg.removeMarked();
        logger.info() << "Filling reliable vertices" << std::endl;
        vertexReliableFiller.loggedRefill(logger, spg);
        logger.info() << "Correcting read paths" << std::endl;
        ag::ErrorCorrectionEngine(reliablePathCorrector).run(logger, threads, spg, dbg_storage);
        logger.info() << "Removing uncovered vertices" << std::endl;
        size_t removed_vertices = ag::SimpleRemoveUncovered(logger, threads, spg);
        logger.info() << "Removed " << removed_vertices << " uncovered vertices" << std::endl;
        // TODO: make core queue in multiplexing listen to merging and simplify/remove reset
        multiplexer.reset();
        SPGCoverage.printStatistics(logger.trace(), spg);
        printer.printDot(dir / ("supregraph_" + itos(threshold, 5) + ".dot"), spg);
    }
    for (Vertex &v: spg.vertices())
        v.reliability = ag::VertexReliability::unknown;
    ag::MergeAllSPG(logger, debug ? 1 : threads, spg);
    CleanSupregraph(spg);
    dbg_storage.stopTrackSuffixes();
    spg.resetEdgeCodes(logger, threads);
    SPGCoverage.printStatistics(logger.trace(), spg);
    logger.info() << "Printing final graph" << std::endl;
    printer.printDot(dir / "supregraph_final.dot", spg);
    dbg_storage.getReads().Save(dir/"reads.aln");
    ag::Printer(ag::VertexPrintStyles::defaultLabeler()).printDirectGFA(dir / "supregraph_final.gfa", spg);
    printer.DrawSplit(ag::Component(spg), dir/"split", 30000);
    return {{"supregraph_final", dir / "supregraph_final.gfa"}};
}

spg::MultiplexAndCorrectionPhase::MultiplexAndCorrectionPhase() : Stage(AlgorithmParameters(
        {"k-mer-size=501", "window=2000", "initial-core-length=600", "max-core-length=5000", "core-length-step=300",
         "vertex-coverage-threshold=4"},
        {}, ""), {"graph", "reads", "extra_reads", "paths"}, {"supregraph_final"}) {
}

std::unordered_map<std::string, std::experimental::filesystem::path>
spg::MultiplexAndCorrectionPhase::innerRun(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                               bool debug, const AlgorithmParameterValues &parameterValues,
                               const std::unordered_map<std::string, io::Library> &input) {
    size_t k = std::stoull(parameterValues.getValue("k-mer-size"));
    size_t w = std::stoull(parameterValues.getValue("window"));
    size_t initial_core_length = std::stoull(parameterValues.getValue("initial-core-length"));
    size_t max_core_length = std::stoull(parameterValues.getValue("max-core-length"));
    size_t core_length_step = std::stoull(parameterValues.getValue("core-length-step"));
    double vertex_coverage_threshold = std::stod(parameterValues.getValue("vertex-coverage-threshold"));
    return RunMultiplexingAndCorrection(logger, threads, dir, k, w, input.at("graph"), input.at("reads"),
        input.at("extra_reads"), input.at("paths"), initial_core_length, max_core_length, core_length_step,
        vertex_coverage_threshold, debug);
}
