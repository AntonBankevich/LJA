#include "path_dumping.hpp"
#include "dbg_read_alignment_storage.hpp"
#include "graph_printing.hpp"
using namespace ag;
size_t stage_num = 0;

void dbg::AddContigs(ag::DLLAlignmentStorage &storage, const KmerIndex &index, const std::vector<Contig> &contigs) {
    for(const Contig &c : contigs) {
        Contig contig = c;
        std::vector<ag::AlignmentChain<Contig, ag::Edge>> als = index.carefulAlign(contig);
        storage.addContig(contig, als);
    }
}

void dbg::PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const string &stage,
                     dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &readStorage, const io::Library &paths_lib, const io::Library &references_lib,
                     bool small) {
    stage_num += 1;
    Printer printer;
    ObjInfo<ag::Edge> edge_printing_style = EdgePrintStyles::defaultDotLabeler();
    if(readStorage.tracksSuffixes())
        edge_printing_style = edge_printing_style + ObjInfo<ag::Edge>::Tooltiper(readStorage.getSuffixes().labeler());
    printer.setEdgeInfo(edge_printing_style);
    printer.setVertexInfo(VertexPrintStyles::defaultDotInfo());
    std::string stage_name = itos(stage_num) + "_" + stage;
    logger.info() << "Dumping current state. Stage id: " << stage_name << std::endl;
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
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
    std::vector<Contig> references = io::SeqReader(references_lib).readAllAsContigs();
    if(paths.empty() && references.empty()) {
        printer.printDot(dir / (stage_name + ".dot"), dbg);
        return;
    }
    dbg::KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, 500);
    ag::DLLAlignmentStorage storage(dbg);
    AddContigs(storage, index, paths);
    AddContigs(storage, index, references);
    edge_printing_style = storage.getEdgeInfo() + edge_printing_style;
    printer.setEdgeInfo(edge_printing_style);
    printer.printDot(dir / (stage_name + ".dot"), dbg);
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getInnerId());
        const std::vector<ag::AlignmentChain<Contig, ag::Edge>> contig_al = index.carefulAlign(contig);
        Component comp = small ? Component::neighbourhood(dbg, contig_al, 1000, 300) :
                              Component::longEdgeNeighbourhood(dbg, contig_al, 20000, 300);
        printer.printDot(dir / "paths" / contig.getInnerId() / (stage_name + ".dot"), comp);
    }
}
