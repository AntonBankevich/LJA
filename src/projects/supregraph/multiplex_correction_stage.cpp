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
                     size_t reliable_read_count, size_t suspicious_read_count, bool debug) {
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
                                                                              false);
    dbg_storage.trackSuffixes(logger, threads, spg, 0, 10000000);
//    SPGCoverage stays attached as a listener for the whole run: it keeps outgoing_read_count
//    correct through both the multiplexing (ResolutionListener) and correction (AlignedReadStorageListener)
//    steps below, which is what lets the two be interleaved.
    ag::AlignedReadStatisticsTracker SPGCoverage(logger, threads, spg, dbg_storage, dbg_storage.getSuffixes());
//    Separate listener computing per-vertex SPG coverage (see CoverageSamplingTracker); must be attached
//    before the edgeToSupreVertex loop below so it can seed samples from the original DBG edges.
    ag::CoverageSamplingTracker spgVertexCoverage(spg, dbg_storage, k, threads);
    spg.disableHashing();
    ag::EdgeInfo edge_info = ag::EdgePrintStyles::defaultDotInfo() + ag::EdgeInfo::Labeler(SPGCoverage.getEdgeLabeler());
    ag::Printer printer(ag::VertexPrintStyles::defaultDotInfo() + ag::VertexInfo::Labeler(SPGCoverage.getVertexLabeler()) +
                         ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()), edge_info);
    ag::DLLAlignmentStorage path_storage(spg);
    std::unique_ptr<spg::PathTracker> path_tracker;
    if (debug) {
        PrepareDLLPathTracker(logger, threads, spg, w, paths, path_storage);
        path_tracker = std::make_unique<spg::PathTracker>(spg, path_storage, printer, dir / "path_tracking");
    } else {
        path_storage.detach();
    }
    std::experimental::filesystem::path figs = dir / "figs";
    if (debug) {
        recreate_dir(figs);
        printer.printDot(figs / "supregraph_initial.dot", spg);
    }
    std::vector<ag::EdgeId> eids = oneline::map(spg.edgesUnique().begin(), spg.edgesUnique().end(), IdTransformer<Edge>());
    UniqueVertexStorage unique_storage(spg, 40000);
    for (ag::EdgeId eid : eids) {
        Vertex &new_vertex = spg.edgeToSupreVertex(*eid);
    }
    dbg_storage.stopTrackCoverage();

    Precorrector precorrector(
            [reliable_read_count](const ag::Edge &e) { return e.outgoing_read_count + e.rc().outgoing_read_count >= reliable_read_count; },
            [suspicious_read_count](const ag::Edge &e) { return e.outgoing_read_count + e.rc().outgoing_read_count <= suspicious_read_count; });

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
        multiplexer.setMaxCoreLength(threshold);
        while (multiplexer.hasReadyCore() || multiplexer.hasPendingMerge()) {
            auto res = multiplexer.process(logger, threads);
            if (debug && !res.empty()) {
                logger.trace() << "Operation " << cnt << ": " << res << std::endl;
                printer.printDot(figs / ("supregraph_" + itos(cnt, 5) + ".dot"), ag::Component::neighbourhood(spg, res, 100000, 20));
                cnt++;
            }
        }
        spg.removeMarked();
        logger.info() << "Correcting reads with core-length threshold "
                       << (threshold == std::numeric_limits<size_t>::max() ? std::string("inf") : itos(threshold)) << std::endl;
        ag::ErrorCorrectionEngine(precorrector).run(logger, threads, spg, dbg_storage);
        ag::SimpleRemoveUncovered(logger, threads, spg);
    }
    ag::MergeAllSPG(logger, debug ? 1 : threads, spg);
    CleanSupregraph(spg);
    spg.resetEdgeCodes(logger, threads);
    logger.info() << "Printing final graph" << std::endl;
    printer.printDot(dir / "supregraph_final.dot", spg);
    ag::Printer(ag::VertexPrintStyles::defaultLabeler()).printDirectGFA(dir / "supregraph_final.gfa", spg);
    return {{"supregraph_final", dir / "supregraph_final.gfa"}};
}

spg::MultiplexAndCorrectionPhase::MultiplexAndCorrectionPhase() : Stage(AlgorithmParameters(
        {"k-mer-size=500", "window=500", "initial-core-length=500", "max-core-length=5000", "core-length-step=200",
         "reliable-read-count=4", "suspicious-read-count=1"},
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
    size_t reliable_read_count = std::stoull(parameterValues.getValue("reliable-read-count"));
    size_t suspicious_read_count = std::stoull(parameterValues.getValue("suspicious-read-count"));
    return RunMultiplexingAndCorrection(logger, threads, dir, k, w, input.at("graph"), input.at("reads"),
        input.at("extra_reads"), input.at("paths"), initial_core_length, max_core_length, core_length_step,
        reliable_read_count, suspicious_read_count, debug);
}
