#include "path_dumping.hpp"
#include "dbg_read_alignment_storage.hpp"
#include "graph_printing.hpp"
using namespace ag;
size_t stage_num = 0;

void dbg::AlignedContigStorage::addContig(const Contig &contig) {
    stored_contigs.emplace_back(new Contig(contig));
    stored_contigs.emplace_back(new Contig(contig.RC()));
}

void dbg::AlignedContigStorage::Fill(size_t threads, dbg::KmerIndex &index) {
    omp_set_num_threads(threads);
    ParallelRecordCollector<ag::AlignmentChain<Contig, dbg::Edge>> records(threads);
    // #pragma omp parallel for schedule(dynamic, 10) default(none) shared(stored_contigs, records, index)
    for(size_t i = 0; i < stored_contigs.size(); i++) {
        Contig &contig = *stored_contigs[i];
        std::vector<ag::AlignmentChain<Contig, dbg::Edge>> path = index.carefulAlign(contig);
        for(ag::AlignmentChain<Contig, dbg::Edge> &al : path) {
            records.emplace_back(al);
        }
    }
    std::vector<ag::AlignmentChain<Contig, dbg::Edge>> rec_list = records.collect();
    __gnu_parallel::sort(rec_list.begin(), rec_list.end());
    std::vector<std::pair<ag::ConstEdgeId , std::vector<ag::AlignmentChain<Contig, dbg::Edge>>>> res;
    std::vector<ag::AlignmentChain<Contig, dbg::Edge>> next;
    for(ag::AlignmentChain<Contig, dbg::Edge> rec : rec_list) {
        if(!next.empty() && next[0].seg_to.contig() != rec.seg_to.contig()) {
            res.emplace_back(next[0].seg_to.contig().getId(), std::move(next));
            next.clear();
        }
        next.emplace_back(rec);
    }
    if(!next.empty()) {
        res.emplace_back(next[0].seg_to.contig().getId(), std::move(next));
    }
    alignments = {res.begin(), res.end()};
}

void dbg::AlignedContigStorage::print(std::ostream &os) {
    for(auto &it : alignments) {
        const dbg::Edge &edge = *it.first;
        os << edge.getInnerId() << "\n";
        if (alignments.find(edge.getId()) == alignments.end())
            return;
        const std::vector<ag::AlignmentChain<Contig, dbg::Edge>> &als = alignments.find(edge.getId())->second;
        if (als.empty()) {
            return;
        }
        os << als[0].seg_from << "->" << als[0].seg_to;
        for (size_t i = 1; i < als.size(); i++) {
            const ag::AlignmentChain<Contig, dbg::Edge> &al = als[i];
            os << "\n" << al.seg_from << "->" << al.seg_to;
        }
    }
}

std::function<std::string(const dbg::Edge &edge)> dbg::AlignedContigStorage::pathInfo() const {
    std::function<std::string(const dbg::Edge &edge)> res = [this](const dbg::Edge &edge) {
        if (alignments.find(edge.getId()) == alignments.end())
            return std::string("");
        std::stringstream ss;
        const std::vector<ag::AlignmentChain<Contig, dbg::Edge>> &als = alignments.find(edge.getId())->second;
        if (als.empty()) {
            return std::string("");
        }
        size_t num = std::min<size_t>(10, als.size());
        ss << als[0].seg_from << "->" << als[0].seg_to.coordinaresStr();
        for (size_t i = 1; i < num; i++) {
            const ag::AlignmentChain<Contig, dbg::Edge> &al = als[i];
            ss << "\\n" << al.seg_from << "->" << al.seg_to.coordinaresStr();
        }
        return ss.str();
    };
    return res;
}

std::function<std::string(const dbg::Edge &edge)> dbg::AlignedContigStorage::colorer(const std::string &color) const {
    std::function<std::string(const dbg::Edge &edge)> res = [this, color](const dbg::Edge &edge) {
        return alignments.find(edge.getId()) == alignments.end() || alignments.find(edge.getId())->second.empty() ?
                   "black" : color;
    };
    return res;
}

void dbg::PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const string &stage,
                     dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &readStorage, const io::Library &paths_lib, const io::Library &references_lib,
                     bool small) {
    stage_num += 1;
    Printer printer;
    ObjInfo<dbg::Edge> edge_printing_style = EdgePrintStyles::defaultDotLabeler();
    if(readStorage.tracksSuffixes())
        edge_printing_style = edge_printing_style + ObjInfo<dbg::Edge>::Tooltiper(readStorage.getSuffixes().labeler());
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
    dbg::AlignedContigStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.addContig(contig);
    }
    for(StringContig sc : io::SeqReader(references_lib)) {
        Contig tmp = sc.makeContig();
        storage.addContig(tmp);
    }
    if(paths.empty() && references_lib.empty()) {
        printer.printDot(dir / (stage_name + ".dot"), dbg);
        return;
    }
    dbg::KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, 500);
    storage.Fill(threads, index);
    edge_printing_style = storage.edgeInfo() + edge_printing_style;
    printer.setEdgeInfo(edge_printing_style);
    printer.printDot(dir / (stage_name + ".dot"), dbg);
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getInnerId());
        const std::vector<ag::AlignmentChain<Contig, dbg::Edge>> contig_al = index.carefulAlign(contig);
        Component comp = small ? Component::neighbourhood(dbg, contig_al, 1000, 300) :
                              Component::longEdgeNeighbourhood(dbg, contig_al, 20000, 300);
        printer.printDot(dir / "paths" / contig.getInnerId() / (stage_name + ".dot"), comp);
    }
}
