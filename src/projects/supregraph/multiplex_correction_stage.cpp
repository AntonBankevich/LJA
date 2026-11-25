#include "multiplex_correction_stage.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "path_tracker.hpp"
#include "path_spy.hpp"
#include "coverage_tracker.hpp"
#include <cmath>
#include <limits>

using namespace spg;

std::unordered_map<std::string, std::experimental::filesystem::path>
spg::RunMultiplexingAndCorrection(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                     size_t k, size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                     const io::Library &extra_reads_files, const io::Library &paths,
                     size_t initial_core_length, size_t max_core_length, size_t core_length_step,
                     double vertex_coverage_threshold, size_t divergence_base, double initial_divergence_rate,
                     double final_divergence_rate, const ag::MultiplicityParams &multiplicity_params, bool debug) {
    size_t unique_threshold = 40000;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    auto [spg, dbg_storage, path_storage] = ConstructSupregraphFromRawInput(logger, threads, dir, k, w, graph_gfa,
                                                                             reads_files, extra_reads_files, paths, debug);
    UniqueVertexStorage unique_storage(spg, unique_threshold);
    OldVertexTracker vertex_tracker(spg, debug);
    spg::CoverageTracker spgVertexCoverage(spg, dbg_storage, logger, threads);
    ag::Printer printer(ag::VertexPrintStyles::spgLabeler() +
                         ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()), ag::EdgePrintStyles::defaultDotInfo());
    printer += path_storage.getPrinter();
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
        printer += path_storage.getPrinter();
    }
    size_t cnt = 1;
    if (debug)
        dbg_storage.logReads(threads, dir/"read_log.txt");
    for (size_t pass = 0; pass < thresholds.size(); pass++) {
        const size_t threshold = thresholds[pass];
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
            if (res.size() > 1)
                mult_cnt++;
        }
        logger.info() << "Multiplexing performed " << mult_cnt << " times (" << iter_cnt << " loop iterations)" << std::endl;
        spg.removeMarked();
        logger.info() << "Filling reliable vertices" << std::endl;
        vertexReliableFiller.loggedRefill(logger, spg);
        const double progress = thresholds.size() == 1 ? 1.0 : double(pass) / double(thresholds.size() - 1);
        const double rate = initial_divergence_rate * std::pow(final_divergence_rate / initial_divergence_rate, progress);
        ag::AlignedSearchParams search_params;
        search_params.divergence = {divergence_base, rate};
        ag::ReliablePathCorrector reliablePathCorrector(search_params);
        logger.info() << "Correcting read paths, allowed divergence " << divergence_base << " + " << rate << " * length"
                      << std::endl;
        ag::ErrorCorrectionEngine(reliablePathCorrector).run(logger, threads, spg, dbg_storage);
        logger.info() << "Removing uncovered vertices" << std::endl;
        size_t removed_vertices = spg::SimpleRemoveUncovered(logger, threads, spg);
        logger.info() << "Removed " << removed_vertices << " uncovered vertices" << std::endl;
        // TODO: make core queue in multiplexing listen to merging and simplify/remove reset
        multiplexer.reset();
        printer.printDot(dir / ("supregraph_" + itos(threshold, 5) + ".dot"), spg);
    }
    for (Vertex &v: spg.vertices())
        v.reliability = ag::VertexReliability::unknown;
    spg::MergeAllSPG(logger, debug ? 1 : threads, spg);
    CleanSupregraph(spg);
    dbg_storage.stopTrackSuffixes();
    spg.resetEdgeCodes(logger, threads);
    logger.info() << "Estimating multiplicities of long vertices" << std::endl;
    ag::LongVertexMultiplicities multiplicities =
            ag::LongVertexMultiplicities::Estimate(logger, threads, spg, multiplicity_params);
//    Paths are tracked only in debug mode (see ConstructSupregraphFromRawInput).
    ag::LongVertexMultiplicities::ReferenceCounter reference;
    if (debug && !paths.empty()) {
        reference = ag::LongVertexMultiplicities::FromPaths(path_storage);
        multiplicities.compareWithReference(logger, reference);
    }
    multiplicities.printReport(dir / "multiplicities.tsv", reference);
    multiplicities.printDeadEnds(dir / "dead_ends.tsv", reference);
    printer += multiplicities.vertexInfo();
    logger.info() << "Printing final graph" << std::endl;
    printer.printDot(dir / "supregraph_final.dot", spg);
    dbg_storage.getReads().Save(dir/"reads.aln");
    auto alignments = PrintStorageAlignments(logger, threads, dbg_storage.getReads(), vertex_tracker, dir / "alignments");
    ag::Printer(ag::VertexPrintStyles::defaultLabeler()).printDirectGFA(dir / "supregraph_final.gfa", spg);
    printer.DrawSplit(ag::Component(spg), dir/"split", 10000);
//    path_storage detaches itself automatically before spg is destroyed (see SupregraphFromRawInput's
//    doc comment on struct member ordering).
    return {{"supregraph_final", dir / "supregraph_final.gfa"}, {"read_alignments", alignments.first}};
}

spg::MultiplexAndCorrectionPhase::MultiplexAndCorrectionPhase() : Stage(AlgorithmParameters(
        {"k-mer-size=501", "window=2000", "initial-core-length=600", "max-core-length=5000", "core-length-step=300",
         "vertex-coverage-threshold=4", "divergence-base=30", "initial-divergence-rate=0.001",
         "final-divergence-rate=0.05", "multiplicity-min-long=40000", "multiplicity-max-long=2000",
         "unit-coverage=0"},
        {}, ""), {"graph", "reads", "extra_reads", "paths"}, {"supregraph_final", "read_alignments"}) {
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
    size_t divergence_base = std::stoull(parameterValues.getValue("divergence-base"));
    double initial_divergence_rate = std::stod(parameterValues.getValue("initial-divergence-rate"));
    double final_divergence_rate = std::stod(parameterValues.getValue("final-divergence-rate"));
    ag::MultiplicityParams multiplicity_params;
    multiplicity_params.min_long = std::stoull(parameterValues.getValue("multiplicity-min-long"));
    multiplicity_params.max_long_per_problem = std::stoull(parameterValues.getValue("multiplicity-max-long"));
    multiplicity_params.unit_coverage = std::stod(parameterValues.getValue("unit-coverage"));
    return RunMultiplexingAndCorrection(logger, threads, dir, k, w, input.at("graph"), input.at("reads"),
        input.at("extra_reads"), input.at("paths"), initial_core_length, max_core_length, core_length_step,
        vertex_coverage_threshold, divergence_base, initial_divergence_rate, final_divergence_rate,
        multiplicity_params, debug);
}
