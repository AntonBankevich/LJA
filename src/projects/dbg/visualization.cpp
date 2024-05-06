#include "visualization.hpp"
#include "graph_alignment_storage.hpp"
#include "graph_printing.hpp"

size_t stage_num = 0;

void PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const string &stage,
           dbg::SparseDBG &dbg, dbg::ReadAlignmentStorage &readStorage, const io::Library &paths_lib,
           bool small) {
    stage_num += 1;
    std::string stage_name = itos(stage_num) + "_" + stage;
    logger.info() << "Dumping current state. Stage id: " << stage_name << std::endl;
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
    printDot(dir / (stage_name + ".dot"), dbg::Component(dbg), readStorage.labeler());
    dbg::printFasta(dir / (stage_name + ".fasta"), dbg);
    if(!small)
        readStorage.printFullAlignments(logger, dir / (stage_name + ".als"));
    std::vector<Contig> paths;
    for(StringContig sc : io::SeqReader(paths_lib)) {
        Contig contig = sc.makeContig();
        if(contig.truncSize() > 100000) {
            paths.emplace_back(contig.getSeq().Subseq(0, 50000), contig.getInnerId() + "_start");
            paths.emplace_back(contig.getSeq().Subseq(contig.truncSize() - 50000), contig.getInnerId() + "_end");
        } else {
            paths.emplace_back(std::move(contig));
        }
    }
    GraphPathStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.addContig(contig);
    }
    if(paths.empty())
        return;
    dbg::KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, 500);
    storage.Fill(threads, index);
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getInnerId());
        const std::vector<ag::AlignmentChain<Contig, dbg::Edge>> contig_al = index.carefulAlign(contig);
        dbg::Component comp = small ? dbg::Component::neighbourhood(dbg, contig_al, 1000) :
                              dbg::Component::longEdgeNeighbourhood(dbg, contig_al, 20000);
        std::function<std::string(dbg::Edge &)> labeler = readStorage.labeler() + storage.labeler();
        printDot(dir / "paths" / contig.getInnerId() / (stage_name + ".dot"), comp, labeler);
    }
}
