#include "multiplexing_stage.hpp"
#include "assembly_graph/dll_path.hpp"

#include "path_spy.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "assembly_graph/data_structures/aligned_read_statistics_tracker.hpp"

void spg::CleanSupregraph(ag::AssemblyGraph &dbg) {
    std::vector<VertexId> ends;
    for(Vertex &v : dbg.vertices())
        if(v.inDeg() == 0 && v.outDeg() == 1 && v.front().isPrefix())
            ends.emplace_back(v.getId());
    for(VertexId vid : ends) {
        while(vid->inDeg() == 0 && vid->outDeg() == 1 && vid->front().isPrefix()) {
            VertexId next = vid->front().getFinish().getId();
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

void PrintConnectedComponents(ag::Printer &printer, const std::experimental::filesystem::path &split, dbg::SparseDBG &spg) {
    recreate_dir(split);
    size_t ccnt = 1;
    for(const ag::Component &component : ag::CCSplitter().split(ag::Component(spg))) {
        printer.printDot(split / (itos(ccnt) + ".dot"), component);
        ccnt++;
    }
}

void PreparePathTracker(logging::Logger &logger, size_t threads, dbg::SparseDBG &spg, size_t w,
            const io::Library &paths, spg::OldPathTracker &path_tracker) {
    std::vector<Contig> contigs = io::SeqReader(paths).readAllAsContigs();
    dbg::KmerIndex index(spg);
    index.fillAnchors(logger, threads, spg, w);
    for (Contig &contig : contigs) {
        std::vector<ag::AlignmentChain<Contig, ag::Edge>> al = index.carefulAlign(contig);
        path_tracker.addPath(contig.getInnerId(), al);
    }
}

void PrepareDLLPathTracker(logging::Logger &logger, size_t threads, dbg::SparseDBG &spg, size_t w,
            const io::Library &paths, ag::DLLAlignmentStorage &path_tracker) {
    std::vector<Contig> contigs = io::SeqReader(paths).readAllAsContigs();
    dbg::KmerIndex index(spg);
    index.fillAnchors(logger, threads, spg, w);
    for (Contig &contig : contigs) {
        std::vector<ag::AlignmentChain<Contig, ag::Edge>> al = index.carefulAlign(contig);
        path_tracker.addContig(contig, al);
    }
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
    hashing::RollingHash hasher(k, 239);
    logger.info() << "Loading graph" << std::endl;
    // TODO: get rid of DBD here completely
    dbg::SparseDBG spg = graph_gfa.empty() ? DBGPipeline(logger, hasher, w, reads_files, dir, threads) :
                         dbg::LoadDBGFromEdgeSequences(logger, threads, graph_gfa, hasher);
    logger.info() << "Loading reads" << std::endl;
    UniqueClassificator classificator = ConstructUnique(logger, threads, reads_files, spg, dir);
    dbg::DBGAlignedReadStorage dbg_storage = dbg::DBGAlignedReadStorage::Load(logger, threads,
                                                                              reads_files + extra_reads_files, spg,
                                                                              false);
    dbg_storage.trackSuffixes(logger, threads, spg, 0, 10000000);
    ag::AlignedReadStatisticsTracker SPGCoverage(logger, threads ,spg, dbg_storage, dbg_storage.getSuffixes());
    spg.disableHashing();
    ag::EdgeInfo edge_info = ag::EdgePrintStyles::defaultDotInfo() + ag::EdgeInfo::Labeler(SPGCoverage.getEdgeLabeler());
                                  ag::EdgeInfo::Tooltiper(dbg_storage.getSuffixes().labeler());
    ag::Printer printer(ag::VertexPrintStyles::defaultDotInfo() + ag::VertexInfo::Labeler(SPGCoverage.getVertexLabeler()), edge_info);
    std::experimental::filesystem::path figs = dir/ "figs";
    recreate_dir(figs);
    printer.printDot(figs / "supregraph_initial.dot", spg);
    std::vector<ag::EdgeId> eids = oneline::map(spg.edgesUnique().begin(), spg.edgesUnique().end(), IdTransformer<Edge>());
    UniqueVertexStorage unique_storage(spg);
    OldVertexTracker vertex_tracker(spg, debug);
    OldPathTracker path_tracker(spg, vertex_tracker, printer, dir/"state_dump");
    ag::DLLAlignmentStorage dll_tracker(spg);

    if (debug) {
        // PreparePathTracker(logger, threads, spg, w, paths, path_tracker);
        PrepareDLLPathTracker(logger, threads, spg, w, paths, dll_tracker);
        dll_tracker.print(logger.debug());
    } else {
        vertex_tracker.detach();
        path_tracker.detach();
        dll_tracker.detach();
    }
    ag::LoggingListener modificationLogger(spg, logger.getLoggerStream(logging::LogLevel::trace));
    if (!debug)
        modificationLogger.detach();
    // TODO: make multiplexing work with any graph. Then somehow avoid the code below.
    for(ag::EdgeId eid : eids) {
        Vertex &new_vertex = spg.edgeToSupreVertex(*eid);
        if(classificator.isUnique(*eid))
            unique_storage.add(new_vertex);
    }
    dbg_storage.stopTrackCoverage();
    logger.info() << "Multiplexing" << std::endl;
//    ChainRule rule(path_index, 4000);
    AndreyRule rule(dbg_storage.getSuffixes(), unique_storage);
    spg::Multiplexer multiplexer(spg, dbg_storage, rule, 200000);
//    multiplexer.fullMultiplex(logger, threads);
    printer.setEdgeInfo(ag::EdgePrintStyles::spgLabeler()
        + ag::EdgeInfo::Labeler(SPGCoverage.getEdgeLabeler())
        + ag::EdgePrintStyles::simpleColorer("black")
        + ag::EdgeInfo::Tooltiper(dbg_storage.getSuffixes().labeler()));
    printer.setVertexInfo(ag::VertexPrintStyles::spgLabeler() + ag::VertexInfo::Labeler(SPGCoverage.getVertexLabeler())+
        ag::VertexPrintStyles::defaultDotColorer() + ag::VertexPrintStyles::defaultTooltiper() +
        ag::VertexInfo::Colorer(unique_storage.getColorer("white", "green")));
    printer += dll_tracker.getPrinter();
    printer.printDot(dir/"initial.dot", spg);
    size_t cnt = 1;
    while (!multiplexer.finished()) {
        auto res = multiplexer.process(logger, threads);
        if(debug && !res.empty()) {
            logger.trace() << "Operation " << cnt << ": " << res << std::endl;
            printer.printDot(figs / ("supregraph_" + itos(cnt, 5) + ".dot"), ag::Component::neighbourhood(spg, res, 100000, 20));
            cnt++;
        }
    }
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
