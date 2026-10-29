#include "multiplex_correction_stage.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "path_tracker.hpp"
#include "sample_segment_tracker.hpp"
#include "coverage_tracker.hpp"
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
    auto [spg, dbg_storage, path_storage] = ConstructSupregraphFromRawInput(logger, threads, dir, k, w, graph_gfa,
                                                                             reads_files, extra_reads_files, paths, debug);
    if (debug) {
        logger.info() << "[DIAG] checking consistency right after ConstructSupregraphFromRawInput" << std::endl;
        VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed right after ConstructSupregraphFromRawInput");
    }
    UniqueVertexStorage unique_storage(spg, unique_threshold);
    spg::CoverageTracker spgVertexCoverage(spg, dbg_storage, k, logger, threads);
    if (debug) {
        logger.info() << "[DIAG] checking consistency right after UniqueVertexStorage/CoverageTracker construction" << std::endl;
        VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed right after UniqueVertexStorage/CoverageTracker construction");
    }
    ag::Printer printer(ag::VertexPrintStyles::spgLabeler() +
                         ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()), ag::EdgePrintStyles::defaultDotInfo());
    if (debug)
        path_storage.startDrawing(dir / "path_tracking", printer);
    std::experimental::filesystem::path figs = dir / "figs";
    if (debug) {
        recreate_dir(figs);
        ag::Printer dbg_printer(ag::VertexPrintStyles::defaultDotInfo() +
                         ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()), ag::EdgePrintStyles::defaultDotInfo());
        dbg_printer.printDot(figs / "dbg_initial.dot", spg);
    }
    dbg_storage.trackSuffixes(logger, threads, spg, 0, 10000000);

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
        printer.setEdgeInfo(ag::EdgePrintStyles::spgLabeler() +
                             ag::EdgePrintStyles::simpleColorer("black"));
        printer.setVertexInfo(ag::VertexPrintStyles::spgLabeler() +
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
        size_t iter_cnt = 0;
        while (multiplexer.hasReadyCore() || multiplexer.hasPendingMerge()) {
            auto res = multiplexer.process(logger, threads);
            if (debug && !res.empty()) {
                logger.trace() << "Operation " << cnt << ": " << res << std::endl;
                // printer.printDot(figs / ("supregraph_" + itos(cnt, 5) + ".dot"), ag::Component::neighbourhood(spg, res, 100000, 20));
                cnt++;
            }
            iter_cnt++;
//            DIAGNOSTIC: every iteration, to localize exactly which operation breaks the index.
            if (debug) {
                logger.info() << "[DIAG] checking after multiplex iteration " << iter_cnt << " (threshold " << threshold << "), last op: " << res << std::endl;
                VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed after multiplex iteration " << iter_cnt << " (threshold " << threshold << "), last op: " << res);
            }
            if (res.size() > 1)
                mult_cnt++;
        }
        logger.info() << "Multiplexing performed " << mult_cnt << " times (" << iter_cnt << " loop iterations)" << std::endl;
        spg.removeMarked();
        if (debug) {
            logger.info() << "[DIAG] checking after removeMarked (threshold " << threshold << ")" << std::endl;
            VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed after removeMarked (threshold " << threshold << ")");
        }
        logger.info() << "Filling reliable vertices" << std::endl;
        vertexReliableFiller.loggedRefill(logger, spg);
        if (debug) {
            logger.info() << "[DIAG] checking after loggedRefill (threshold " << threshold << ")" << std::endl;
            VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed after loggedRefill (threshold " << threshold << ")");
        }
        logger.info() << "Correcting read paths" << std::endl;
        ag::ErrorCorrectionEngine(reliablePathCorrector).run(logger, threads, spg, dbg_storage);
        if (debug) {
            logger.info() << "[DIAG] checking after ErrorCorrectionEngine (threshold " << threshold << ")" << std::endl;
            VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed after ErrorCorrectionEngine (threshold " << threshold << ")");
        }
        logger.info() << "Removing uncovered vertices" << std::endl;
        size_t removed_vertices = ag::SimpleRemoveUncovered(logger, threads, spg);
        logger.info() << "Removed " << removed_vertices << " uncovered vertices" << std::endl;
        if (debug) {
            logger.info() << "[DIAG] checking after SimpleRemoveUncovered (threshold " << threshold << ")" << std::endl;
            VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed after SimpleRemoveUncovered (threshold " << threshold << ")");
        }
        // TODO: make core queue in multiplexing listen to merging and simplify/remove reset
        multiplexer.reset();
        if (debug) {
            logger.info() << "[DIAG] checking after multiplexer.reset (threshold " << threshold << ")" << std::endl;
            VERIFY_MSG(dbg_storage.checkConsistency(), "[DIAG] consistency check failed after multiplexer.reset (threshold " << threshold << ")");
        }
        printer.printDot(dir / ("supregraph_" + itos(threshold, 5) + ".dot"), spg);
    }
    for (Vertex &v: spg.vertices())
        v.reliability = ag::VertexReliability::unknown;
    ag::MergeAllSPG(logger, debug ? 1 : threads, spg);
    CleanSupregraph(spg);
    dbg_storage.stopTrackSuffixes();
    spg.resetEdgeCodes(logger, threads);
    logger.info() << "Printing final graph" << std::endl;
    printer.printDot(dir / "supregraph_final.dot", spg);
    dbg_storage.getReads().Save(dir/"reads.aln");
    ag::Printer(ag::VertexPrintStyles::defaultLabeler()).printDirectGFA(dir / "supregraph_final.gfa", spg);
    printer.DrawSplit(ag::Component(spg), dir/"split", 30000);
//    path_storage detaches itself automatically before spg is destroyed (see SupregraphFromRawInput's
//    doc comment on struct member ordering).
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
