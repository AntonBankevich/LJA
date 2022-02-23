#include "pipeline.hpp"
#include "repeat_resolution/repeat_resolution.hpp"

using namespace dbg;
void multigraph::LJAPipeline::PrintPaths(logging::Logger &logger, const std::experimental::filesystem::path &dir, const std::string &stage,
                SparseDBG &dbg, RecordStorage &readStorage, const io::Library &paths_lib, bool small) {
    stage_num += 1;
    std::string stage_name = itos(stage_num) + "_" + stage;
    logger.info() << "Dumping current state. Stage id: " << stage_name << std::endl;
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
    printDot(dir / (stage_name + ".dot"), Component(dbg), readStorage.labeler());
    dbg.printFastaOld(dir / (stage_name + ".fasta"));
    if(!small)
        readStorage.printFullAlignments(logger, dir / (stage_name + ".als"));
    std::vector<Contig> paths;
    for(StringContig sc : io::SeqReader(paths_lib)) {
        Contig contig = sc.makeContig();
        if(contig.size() > 100000) {
            paths.emplace_back(contig.seq.Subseq(0, 50000), contig.id + "_start");
            paths.emplace_back(contig.seq.Subseq(contig.size() - 50000), contig.id + "_end");
        } else {
            paths.emplace_back(std::move(contig));
        }
    }
    GraphAlignmentStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.fill(contig);
    }
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getId());
        Component comp = small ? Component::neighbourhood(dbg, contig, dbg.hasher().getK() + 500) :
                         Component::longEdgeNeighbourhood(dbg, contig, 20000);
        std::function<std::string(dbg::Edge &)> labeler = readStorage.labeler() + storage.labeler();
        printDot(dir / "paths" / contig.getId() / (stage_name + ".dot"), comp, labeler);
    }
    std::ofstream ref_os;
    ref_os.open(dir / (stage_name + ".ref"));
    for(Contig &contig : ref){
        ref_os << contig.getId() << std::endl;
        for(const PerfectAlignment<Contig, dbg::Edge> &al : GraphAligner(dbg).carefulAlign(contig)) {
            ref_os << al.seg_from << " " << al.seg_to << std::endl;
        }
    }
    ref_os.close();
}

std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
multigraph::LJAPipeline::AlternativeCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                      const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
                      size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
                      bool close_gaps, bool remove_bad, bool skip, bool debug, bool load) {
    logger.info() << "Performing initial correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, close_gaps, load, remove_bad, k, w, &reads_lib,
            &pseudo_reads_lib, &paths_lib, threads, threshold, reliable_coverage, debug, this] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                        DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 2, 1000);
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true, false);
        RecordStorage refStorage(dbg, 0, extension_size, threads, readLogger, false, false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        coverageStats(logger, dbg);
        if(debug) {
            PrintPaths(logger, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
        }
        Precorrect(logger, threads, dbg, readStorage, reliable_coverage);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, extension_size);
        readStorage.trackSuffixes(logger, threads);
//        CorrectDimers(logger, readStorage, k, threads, reliable_coverage);
        correctAT(logger, readStorage, k, threads);
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 800, 4, threads);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mk800", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 5 / 2, 3000));
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 2000, 4, threads);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mk2000", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 7 / 2, 5000));
//        CorrectDimers(logger, readStorage, k, threads, reliable_coverage);
        correctAT(logger, readStorage, k, threads);
        correctLowCoveredRegions(logger, dbg, readStorage, refStorage, "/dev/null", threshold, reliable_coverage, k, threads, false);
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 3500, 4, threads);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        coverageStats(logger, dbg);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mk3500", dbg, readStorage, paths_lib, false);
        readStorage.printReadFasta(logger, dir / "corrected.fasta");
        if(debug)
            DrawSplit(Component(dbg), dir / "split");
        dbg.printFastaOld(dir / "graph.fasta");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected.fasta";
    logger.info() << "Initial correction results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "graph.fasta"};
}

std::vector<std::experimental::filesystem::path> multigraph::LJAPipeline::NoCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                                                              const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
                                                              size_t threads, size_t k, size_t w, bool skip, bool debug, bool load) {
    logger.info() << "Performing initial correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, load, k, w, &reads_lib,
            &pseudo_reads_lib, &paths_lib, threads, debug, this] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                        DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 2, 1000);
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true, false);
        RecordStorage extra_reads(dbg, 0, extension_size, threads, readLogger, false, true, false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        coverageStats(logger, dbg);
        if(debug) {
            PrintPaths(logger, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
        }
        dbg.printFastaOld(dir / "final_dbg.fasta");
        printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler());
        printGFA(dir / "final_dbg.gfa", Component(dbg), true);
        SaveAllReads(dir/"final_dbg.aln", {&readStorage, &extra_reads});
        readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");
    };
    if(!skip)
        runInFork(ic_task);

    return {dir/"corrected_reads.fasta", dir / "final_dbg.fasta", dir / "final_dbg.aln"};
}

std::vector<std::experimental::filesystem::path> multigraph::LJAPipeline::SecondPhase(
        logging::Logger &logger, const std::experimental::filesystem::path &dir,
        const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
        const io::Library &paths_lib, size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
        size_t unique_threshold, bool diploid, bool skip, bool debug, bool load) {
    logger.info() << "Performing second phase of error correction using k = " << k << std::endl;
    if (k%2==0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1)
                      << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, load, k, w,
            &reads_lib, &pseudo_reads_lib, &paths_lib,
            threads, threshold, reliable_coverage,
            debug, unique_threshold, diploid, this]
    {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg =
                load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads,
                                   (dir/"disjointigs.fasta").string(),
                                   (dir/"vertices.save").string())
                     : DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = 10000000;
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, debug);
        RecordStorage refStorage(dbg, 0, extension_size, threads, readLogger, false, false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        if(debug) {
            DrawSplit(Component(dbg), dir / "before_figs", readStorage.labeler(), 25000);
            PrintPaths(logger, dir / "state_dump", "initial", dbg, readStorage, paths_lib, false);
        }
        initialCorrect(dbg, logger, dir / "correction.txt", readStorage, refStorage,
                       threshold, 2 * threshold, reliable_coverage, threads, false);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "low", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, threads, dbg, {&readStorage, &refStorage});
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "gap1", dbg, readStorage, paths_lib, false);
        readStorage.invalidateBad(logger, threads, threshold, "after_gap1");
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "bad", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "uncovered1", dbg, readStorage, paths_lib, false);
        RecordStorage extra_reads = MultCorrect(dbg, logger, dir, readStorage, unique_threshold, threads, diploid, debug);
        MRescue(logger, threads, dbg, readStorage, unique_threshold, 0.05);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mult", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &extra_reads, &refStorage});
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "uncovered2", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, threads, dbg, {&readStorage, &extra_reads, &refStorage});
        if(debug) {
            PrintPaths(logger, dir / "state_dump", "gap2", dbg, readStorage, paths_lib, false);
            DrawSplit(Component(dbg), dir / "split_figs", readStorage.labeler());
        }
        dbg.printFastaOld(dir / "final_dbg.fasta");
        printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler());
        printGFA(dir / "final_dbg.gfa", Component(dbg), true);
        SaveAllReads(dir/"final_dbg.aln", {&readStorage, &extra_reads});
        readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected_reads.fasta";
    logger.info() << "Second phase results with k = " << k << " printed to "
                  << res << std::endl;
    return {res, dir / "final_dbg.fasta", dir / "final_dbg.aln"};
}

std::vector<std::experimental::filesystem::path> multigraph::LJAPipeline::MDBGPhase(
        logging::Logger &logger, size_t threads, size_t k, size_t kmdbg, size_t w, size_t unique_threshold, bool diploid,
        const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &graph_fasta,
        const std::experimental::filesystem::path &read_paths, bool skip, bool debug) {
    logger.info() << "Performing repeat resolution by transforming de Bruijn graph into Multiplex de Bruijn graph" << std::endl;
    std::function<void()> ic_task = [&logger, threads, debug, k, kmdbg, &graph_fasta, unique_threshold, diploid, &read_paths, &dir] {
        hashing::RollingHash hasher(k, 239);
        SparseDBG dbg = dbg::LoadDBGFromFasta({graph_fasta}, hasher, logger, threads);
        size_t extension_size = 10000000;
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, debug);
        RecordStorage extra_reads(dbg, 0, extension_size, threads, readLogger, false, debug);
        LoadAllReads(read_paths, {&readStorage, &extra_reads}, dbg);
        repeat_resolution::RepeatResolver rr(dbg, &readStorage, {&extra_reads},
                                             k, kmdbg, dir, unique_threshold,
                                             diploid, debug, logger);
        rr.ResolveRepeats(logger, threads);
    };
    if(!skip)
        runInFork(ic_task);
    return {dir / "assembly.hpc.fasta", dir / "mdbg.hpc.gfa"};
}

std::vector<std::experimental::filesystem::path> multigraph::LJAPipeline::PolishingPhase(
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &output_dir,
        const std::experimental::filesystem::path &gfa_file,
        const std::experimental::filesystem::path &corrected_reads,
        const io::Library &reads, size_t dicompress, size_t min_alignment, bool skip, bool debug) {
    logger.info() << "Performing polishing and homopolymer uncompression" << std::endl;
    std::function<void()> ic_task = [&logger, threads, &output_dir, debug, &gfa_file, &corrected_reads, &reads, dicompress, min_alignment, &dir] {
        io::SeqReader reader(corrected_reads);
        multigraph::MultiGraph vertex_graph;
        vertex_graph.LoadGFA(gfa_file, true);
        multigraph::MultiGraph edge_graph = vertex_graph.DBG();
        std::vector<Contig> contigs = edge_graph.getEdges(false);
        auto res = PrintAlignments(logger, threads, contigs, reader.begin(), reader.end(), min_alignment, dir);
        std::vector<Contig> uncompressed = Polish(logger, threads, contigs, res.first, reads, dicompress);
        std::vector<Contig> assembly = printUncompressedResults(logger, threads, edge_graph, uncompressed, output_dir, debug);
        logger.info() << "Printing final assembly to " << (output_dir / "assembly.fasta") << std::endl;
        std::ofstream os_cut;
        os_cut.open(output_dir / "assembly.fasta");
        for(Contig &contig : assembly) {
            if(contig.size() > 1500)
                os_cut << ">" << contig.id << "\n" << contig.seq << "\n";
        }
        os_cut.close();
    };
    if(!skip)
        runInFork(ic_task);
    return {output_dir / "assembly.fasta", output_dir / "mdbg.gfa"};
}

