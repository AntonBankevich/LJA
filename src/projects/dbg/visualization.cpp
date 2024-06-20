#include "visualization.hpp"
#include "dbg_read_alignment_storage.hpp"
#include "graph_printing.hpp"

size_t stage_num = 0;

void PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const string &stage,
                dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &readStorage, const io::Library &paths_lib, const io::Library &references_lib,
                bool small) {
    stage_num += 1;
    Printer<dbg::DBGTraits> printer;
    ObjInfo<dbg::Edge> edge_printing_style = EdgePrintStyles<dbg::DBGTraits>::defaultDotLabeler();
    if(readStorage.tracksSuffixes())
        edge_printing_style = edge_printing_style + ObjInfo<dbg::Edge>::Tooltiper(readStorage.getSuffixes().labeler());
    printer.setEdgeInfo(edge_printing_style);
    printer.setVertexInfo(VertexPrintStyles<dbg::DBGTraits>::defaultDotInfo());
    std::string stage_name = itos(stage_num) + "_" + stage;
    logger.info() << "Dumping current state. Stage id: " << stage_name << std::endl;
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
    printer.printDot(dir / (stage_name + ".dot"), dbg);
    dbg::printFasta(dir / (stage_name + ".fasta"), dbg);
    if(!small)
        readStorage.getReads().printFullAlignments(logger, dir / (stage_name + ".als"));
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
    GraphAlignedReadStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.addContig(contig);
    }
    for(StringContig sc : io::SeqReader(references_lib)) {
        Contig tmp = sc.makeContig();
        storage.addContig(tmp);
    }
    if(paths.empty())
        return;
    dbg::KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, 500);
    storage.Fill(threads, index);
    edge_printing_style = ObjInfo<dbg::Edge>::Tooltiper(storage.labeler()) + edge_printing_style;
    printer.setEdgeInfo(edge_printing_style);
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getInnerId());
        const std::vector<ag::AlignmentChain<Contig, dbg::Edge>> contig_al = index.carefulAlign(contig);
        dbg::Component comp = small ? dbg::Component::neighbourhood(dbg, contig_al, 1000) :
                              dbg::Component::longEdgeNeighbourhood(dbg, contig_al, 20000);
        printer.printDot(dir / "paths" / contig.getInnerId() / (stage_name + ".dot"), comp);
    }
}
