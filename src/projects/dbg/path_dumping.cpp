#include "path_dumping.hpp"
#include "dbg_read_alignment_storage.hpp"
#include "graph_printing.hpp"
using namespace ag;
size_t stage_num = 0;

void ag::AlignedContigStorage::addContig(Contig contig) {
    contigs.emplace_back(std::move(contig));
    contigs.emplace_back(contigs.back().RC());
}
void ag::AlignedContigStorage::print(std::ostream &os) {
    for(auto &it : edge_alignments) {
        const Edge &edge = *it.first;
        os << edge.getInnerId() << "\n";
        if (edge_alignments.find(edge.getId()) == edge_alignments.end())
            return;
        const std::vector<ag::AlignmentChain<Contig, Edge>> &als = edge_alignments.find(edge.getId())->second;
        if (als.empty()) {
            return;
        }
        os << als[0].seg_from << "->" << als[0].seg_to;
        for (size_t i = 1; i < als.size(); i++) {
            const ag::AlignmentChain<Contig, Edge> &al = als[i];
            os << "\n" << al.seg_from << "->" << al.seg_to;
        }
    }
}

std::function<std::string(const ag::Edge &edge)> ag::AlignedContigStorage::pathInfo() const {
    std::function<std::string(const ag::Edge &edge)> res = [this](const ag::Edge &edge) {
        if (edge_alignments.find(edge.getId()) == edge_alignments.end())
            return std::string("");
        std::stringstream ss;
        const std::vector<ag::AlignmentChain<Contig, ag::Edge>> &als = edge_alignments.find(edge.getId())->second;
        if (als.empty()) {
            return std::string("");
        }
        size_t num = std::min<size_t>(10, als.size());
        ss << als[0].seg_from << "->" << als[0].seg_to.coordinaresStr();
        for (size_t i = 1; i < num; i++) {
            const ag::AlignmentChain<Contig, ag::Edge> &al = als[i];
            ss << "\\n" << al.seg_from << "->" << al.seg_to.coordinaresStr();
        }
        return ss.str();
    };
    return res;
}

std::function<std::string(const ag::Edge &edge)> ag::AlignedContigStorage::colorer(const std::string &color) const {
    std::function<std::string(const ag::Edge &edge)> res = [this, color](const ag::Edge &edge) {
        return edge_alignments.find(edge.getId()) == edge_alignments.end() || edge_alignments.find(edge.getId())->second.empty() ?
                   "black" : color;
    };
    return res;
}

void AlignedContigStorage::Fill(logging::Logger &logger, size_t threads, dbg::KmerIndex &index) {
    omp_set_num_threads(threads);
    ParallelRecordCollector<ag::AlignmentChain<Contig, Edge>> edge_records(threads);
    ParallelRecordCollector<ag::AlignmentChain<Contig, Vertex>> vertex_records(threads);
#pragma omp parallel for default(none) shared(contigs, edge_records, vertex_records, index)
    for(Contig &contig: contigs) {
        std::vector<ag::AlignmentChain<Contig, Edge>> path = index.carefulAlign(contig);
        for(ag::AlignmentChain<Contig, Edge> &al : path) {
            if (al.seg_to.left == 0) {
                Vertex &v = al.seg_to.contig().getStart();
                vertex_records.emplace_back(al.seg_from.contig(), v, al.seg_from.left, 0, v.size());
            }
            if (al.seg_to.right == al.seg_to.contig().truncSize()) {
                Vertex &v = al.seg_to.contig().getFinish();
                vertex_records.emplace_back(al.seg_from.contig(), v, al.seg_from.right, 0, v.size());
            }
            edge_records.emplace_back(al);
        }
    }
    std::vector<ag::AlignmentChain<Contig, ag::Edge>> erec_list = edge_records.collect();
    std::vector<ag::AlignmentChain<Contig, Vertex>> vrec_list = vertex_records.collect();
    __gnu_parallel::sort(erec_list.begin(), erec_list.end());
    __gnu_parallel::sort(vrec_list.begin(), vrec_list.end());
    vrec_list.erase(std::unique(vrec_list.begin(), vrec_list.end()), vrec_list.end());
    std::vector<std::pair<ag::ConstEdgeId , std::vector<ag::AlignmentChain<Contig, ag::Edge>>>> res;
    std::vector<std::pair<ag::ConstVertexId , std::vector<ag::AlignmentChain<Contig, ag::Vertex>>>> res_v;
    edge_alignments = GroupByContig(erec_list);
    vertex_alignments = GroupByContig<Vertex>(vrec_list);
}

void AlignedContigStorage::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    lock();
    std::vector<AlignmentChain<Contig, Vertex>> &supre_alignments = vertex_alignments.at(v.getId());
    for (AlignmentChain<Contig, Edge> &rec : edge_alignments.at(e.getId())) {
        supre_alignments.emplace_back(rec.seg_from.contig(), v, rec.seg_from.left,
                                      rec.seg_to.left, rec.seg_from.size() + e.getStart().size());
    }
    unlock();
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
    ag::AlignedContigStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.addContig(std::move(contig));
    }
    for(StringContig sc : io::SeqReader(references_lib)) {
        Contig tmp = sc.makeContig();
        storage.addContig(std::move(tmp));
    }
    if(paths.empty() && references_lib.empty()) {
        printer.printDot(dir / (stage_name + ".dot"), dbg);
        return;
    }
    dbg::KmerIndex index(dbg);
    index.fillAnchors(logger, threads, dbg, 500);
    storage.Fill(logger, threads, index);
    edge_printing_style = storage.edgeInfo() + edge_printing_style;
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
