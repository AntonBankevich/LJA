#include "visualization.hpp"
#include "graph_alignment_storage.hpp"

size_t stage_num = 0;

void PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const string &stage,
           dbg::SparseDBG &dbg, RecordStorage &readStorage, const io::Library &paths_lib,
           bool small) {
    stage_num += 1;
    std::string stage_name = itos(stage_num) + "_" + stage;
    logger.info() << "Dumping current state. Stage id: " << stage_name << std::endl;
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
    printDot(dir / (stage_name + ".dot"), dbg::Component(dbg), readStorage.labeler());
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
        storage.addContig(contig);
    }
    storage.Fill(threads);
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getId());
        dbg::Component comp = small ? dbg::Component::neighbourhood(dbg, contig, dbg.hasher().getK() + 500) :
                              dbg::Component::longEdgeNeighbourhood(dbg, contig, 20000);
        std::function<std::string(dbg::Edge &)> labeler = readStorage.labeler() + storage.labeler();
        printDot(dir / "paths" / contig.getId() / (stage_name + ".dot"), comp, labeler);
    }
}