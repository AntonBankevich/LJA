#include "multiplexing_stage.hpp"
#include "assembly_graph/dll_path.hpp"

#include "path_spy.hpp"
#include "path_tracker.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "sample_segment_tracker.hpp"
#include "coverage_tracker.hpp"

void spg::CleanSupregraph(ag::AssemblyGraph &dbg) {
    std::vector<VertexId> ends;
    for(Vertex &v : dbg.vertices())
        if(v.inDeg() == 0 && v.outDeg() == 1 && v.front().isPrefix())
            ends.emplace_back(v.getId());
    for(VertexId vid : ends) {
        while(vid->inDeg() == 0 && vid->outDeg() == 1 && vid->front().isPrefix()) {
            VertexId next = vid->frontVertex().getId();
            dbg.isolateAndMark(*vid);
            vid = next;
        }
    }
    dbg.removeMarked();
}

UniqueClassificator
spg::ConstructUnique(logging::Logger &logger, size_t threads, const io::Library &reads_files, dbg::SparseDBG &dbg,
                     const std::experimental::filesystem::path &dir) {
    dbg::DBGAlignedReadStorage dbg_storage = dbg::DBGAlignedReadStorage::Load(logger, threads, reads_files, dbg, true);
    logger.info() << "Constructing path index" << std::endl;
    dbg_storage.trackSuffixes(logger, threads, dbg, 0, 1000000);
    logger.info() << "Printing initial graph" << std::endl;
    logger.info() << "Reconstructing uniqueness" << std::endl;
    UniqueClassificator classificator(dbg, dbg_storage, 0, false, false);
    classificator.classify(logger, 40000, dir / "mult");
    return std::move(classificator);
}

void PrintConnectedComponents(ag::Printer &printer, const std::experimental::filesystem::path &split, ag::AssemblyGraph &spg) {
    recreate_dir(split);
    size_t ccnt = 1;
    for(const ag::Component &component : ag::CCSplitter().split(ag::Component(spg))) {
        printer.printDot(split / (itos(ccnt) + ".dot"), component);
        ccnt++;
    }
}

spg::SupregraphFromDBG
spg::ConstructSupregraphFromDBG(logging::Logger &logger, size_t threads, size_t k, ag::AssemblyGraph &&graph,
                                dbg::DBGAlignedReadStorage &&dbg_storage) {
    Supregraph spg = Supregraph::FromDBG(logger, threads, std::move(graph), k);
    return {std::move(spg), std::move(dbg_storage)};
}

spg::SupregraphFromRawInput
spg::ConstructSupregraphFromRawInput(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                                     size_t k, size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                                     const io::Library &extra_reads_files, const io::Library &paths, bool debug) {
    hashing::RollingHash hasher(k, 239);
    logger.info() << "Loading graph" << std::endl;
    dbg::SparseDBG dbg = graph_gfa.empty() ? DBGPipeline(logger, hasher, w, reads_files, dir, threads) :
                         dbg::LoadDBGFromEdgeSequences(logger, threads, graph_gfa, hasher);
    logger.info() << "Loading reads" << std::endl;
    dbg::DBGAlignedReadStorage dbg_storage = dbg::DBGAlignedReadStorage::Load(logger, threads,
                                                                              reads_files + extra_reads_files, dbg,
                                                                              false);
//    Attached to dbg (not the Supregraph it becomes) so it observes the edgeToSupreVertex conversions too.
    ag::DLLAlignmentStorage path_storage(dbg);
    if (debug) {
//        Must run before disableHashing(): KmerIndex needs dbg's own k-mer hashing, which only exists
//        while the graph is still a plain DBG (uniform k-mer vertices), not after SPG conversion.
        PrepareDLLPathTracker(logger, threads, dbg, w, paths, path_storage);
        path_storage.startLogging(logger.trace());
    } else {
        path_storage.detach();
    }
    dbg.disableHashing();
    Supregraph spg = Supregraph::FromDBG(logger, threads, std::move(dbg), k);
    return {std::move(spg), std::move(dbg_storage), std::move(path_storage)};
}

std::unordered_map<std::string, std::experimental::filesystem::path>
spg::RunMultiplexing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, size_t k,
                     size_t w, const io::Library &graph_gfa, const io::Library &reads_files,
                     const io::Library &extra_reads_files, const io::Library &paths, bool debug) {
    if (k%2==0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1)
                      << " to make it odd" << std::endl;
        k += 1;
    }
    auto [spg, dbg_storage, dll_tracker] = ConstructSupregraphFromRawInput(logger, threads, dir, k, w, graph_gfa,
                                                                            reads_files, extra_reads_files, paths, debug);
//    Built only now, after every edge has been converted: see ConstructSupregraphFromDBG's doc comment
//    for why CoverageTracker/UniqueVertexStorage can't be built any earlier. vertex_tracker/
//    modificationLogger attach here too (rather than to the pre-conversion dbg, like path_storage does
//    inside ConstructSupregraphFromRawInput) -- missing the edgeToSupreVertex conversions is an acceptable
//    gap for these two debug-only trackers, and avoids them outliving `built` in local destruction order.
    UniqueVertexStorage unique_storage(spg);
    spg::CoverageTracker spgVertexCoverage(spg, dbg_storage, k, logger, threads);
    dbg_storage.stopTrackCoverage();
    OldVertexTracker vertex_tracker(spg, debug);
    ag::LoggingListener modificationLogger(spg, logger.getLoggerStream(logging::LogLevel::trace));
    if (!debug) {
        vertex_tracker.detach();
        modificationLogger.detach();
    }
    ag::EdgeInfo edge_info = ag::EdgePrintStyles::defaultDotInfo();
                                  ag::EdgeInfo::Tooltiper(dbg_storage.getSuffixes().labeler());
    ag::Printer printer(ag::VertexPrintStyles::defaultDotInfo() + ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()), edge_info);
    std::experimental::filesystem::path figs = dir/ "figs";
    recreate_dir(figs);
    printer.printDot(figs / "supregraph_initial.dot", spg);
    if (debug) {
        dll_tracker.print(logger.debug());
        dll_tracker.startDrawing(dir / "path_tracking", printer);
    }
    dbg_storage.trackSuffixes(logger, threads, spg, 0, 10000000);
    logger.info() << "Multiplexing" << std::endl;
//    ChainRule rule(path_index, 4000);
    AndreyRule rule(dbg_storage.getSuffixes(), unique_storage);
    spg::Multiplexer multiplexer(spg, dbg_storage, rule, 200000);
//    multiplexer.fullMultiplex(logger, threads);
    printer.setEdgeInfo(ag::EdgePrintStyles::spgLabeler()
        + ag::EdgePrintStyles::simpleColorer("black")
        + ag::EdgeInfo::Tooltiper(dbg_storage.getSuffixes().labeler()));
    printer.setVertexInfo(ag::VertexPrintStyles::spgLabeler() +
        ag::VertexInfo::Labeler(spgVertexCoverage.getVertexLabeler()) +
        ag::VertexPrintStyles::defaultDotColorer() + ag::VertexPrintStyles::defaultTooltiper() +
        ag::VertexInfo::Colorer(unique_storage.getColorer("white", "green")));
    printer += dll_tracker.getPrinter();
    printer.printDot(dir/"initial.dot", spg);
    size_t cnt = 1;
    size_t iter_cnt = 0;
    while (!multiplexer.finished()) {
        auto res = multiplexer.process(logger, threads);
        if(debug && !res.empty()) {
            logger.trace() << "Operation " << cnt << ": " << res << std::endl;
            printer.printDot(figs / ("supregraph_" + itos(cnt, 5) + ".dot"), ag::Component::neighbourhood(spg, res, 100000, 20));
            cnt++;
        }
        iter_cnt++;
//        Tied to actual loop iterations, not to debug/res.empty() -- those gate dot-file printing, not
//        whether a graph edit happened (merge_queue draining, trivial tip removal, etc. all still edit
//        the graph with an empty `res`). Periodic, not every iteration -- checkConsistency() is O(reads)
//        (it walks every read twice over, once to clear, once to restore) and this is strictly a debug aid.
        if (debug && iter_cnt % 20 == 0)
            VERIFY_MSG(dbg_storage.checkConsistency(), "Read storage consistency check failed after iteration " << iter_cnt);
    }
    if (debug)
        logger.info() << "Multiplexing loop ran " << iter_cnt << " iterations, periodic consistency check ran "
                       << (iter_cnt / 20) << " times" << std::endl;
    spg.removeMarked();
    if(debug) PrintConnectedComponents(printer, dir/"final_graph_figs", spg);
    ag::MergeAllSPG(logger, debug ? 1 : threads, spg);
    CleanSupregraph(spg);
    if (debug) {
        for (Vertex &v : spg.vertices()) {
            logger.trace() << v.getId() << " " << v.size() << std::endl;
            for (OldVertexTracker::DirectEmbedding embedding : vertex_tracker.getAllSubvertices(v.getId())) {
                logger.trace() << embedding.inner_vertex.innerId() << "[" << embedding.from << "," << embedding.to << "]" << std::endl;
            }
        }
    }
    spg.resetEdgeCodes(logger, threads);
    logger.info() << "Printing final graph" << std::endl;
    printer.printDot(dir / "supregraph_final.dot", spg);
    if(debug) PrintConnectedComponents(printer, dir / "split", spg);
    ag::Printer(ag::VertexPrintStyles::defaultLabeler()).printDirectGFA(dir / "supregraph_final.gfa", spg);
//    vertex_tracker/modificationLogger are declared after the ConstructSupregraphFromRawInput call (attached
//    to spg directly, not to the pre-conversion dbg), so they already destruct before spg in local
//    destruction order -- no manual detach needed. dll_tracker (bound to the SupregraphFromRawInput
//    result's path_storage member) detaches itself via the same struct-member-ordering guarantee (see
//    SupregraphFromRawInput's doc comment).
    return {{"supregraph_final", dir / "supregraph_final.gfa"}};
}

spg::SupreGraphPhase::SupreGraphPhase() : Stage(AlgorithmParameters(
        {"k-mer-size=5001", "window=500"},
        {}, ""), {"graph", "reads", "extra_reads", "paths"}, {"supregraph_final"}) {
}

std::unordered_map<std::string, std::experimental::filesystem::path>
spg::SupreGraphPhase::innerRun(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                               bool debug, const AlgorithmParameterValues &parameterValues,
                               const std::unordered_map<std::string, io::Library> &input) {
    size_t k = std::stoull(parameterValues.getValue("k-mer-size"));
    size_t w = std::stoull(parameterValues.getValue("window"));
    return RunMultiplexing(logger, threads, dir, k, w, input.at("graph"), input.at("reads"),
        input.at("extra_reads"), input.at("paths"), debug);
}
