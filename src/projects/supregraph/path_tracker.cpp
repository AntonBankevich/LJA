#include "path_tracker.hpp"

using namespace spg;

void spg::PrepareDLLPathTracker(logging::Logger &logger, size_t threads, dbg::SparseDBG &spg, size_t w,
            const io::Library &paths, ag::DLLAlignmentStorage &path_tracker) {
    std::vector<Contig> contigs = io::SeqReader(paths).readAllAsContigs();
    dbg::KmerIndex index(spg);
    index.fillAnchors(logger, threads, spg, w);
    for (Contig &contig : contigs) {
        std::vector<ag::AlignmentChain<Contig, ag::Edge>> al = index.carefulAlign(contig);
        path_tracker.addContig(contig, al);
    }
}

PathTracker::PathFigures &PathTracker::figuresFor(const std::string &contig_name) {
    auto it = paths.find(contig_name);
    if (it != paths.end())
        return it->second;
    return paths.emplace(contig_name, PathFigures(dir / contig_name)).first->second;
}

void PathTracker::drawFiguresForVertex(Vertex &vertex, const std::string &event_tag) {
    for (const std::string &contig_name : storage->passingForwardContigs(vertex)) {
        VERIFY(!startsWith(contig_name, "-"));
        std::experimental::filesystem::path fname = figuresFor(contig_name).nextFile(event_tag);
        ag::Printer p = printer + ag::VertexInfo::Colorer(ag::ConstMapping(vertex, "orange"));
        p.printDot(fname, ag::Component::neighbourhood(getFire<ag::AssemblyGraph>(),
                                                        std::vector<VertexId>{vertex.getId()}, radius, max_size));
    }
}

PathTracker::PathTracker(ag::ResolutionFire &fire, const ag::DLLAlignmentStorage &storage,
                          const ag::Printer &printer, std::experimental::filesystem::path dir,
                          size_t radius, size_t max_size) :
        ag::ResolutionListener(fire, "PathTracker"), storage(&storage), printer(printer),
        dir(std::move(dir)), radius(radius), max_size(max_size) {
    recreate_dir(this->dir);
}

void PathTracker::fireEdgeToSupreVertex(Vertex &v, Edge &e) {
    drawFiguresForVertex(v, "EdgeToSupreVertex_" + std::to_string(v.getInnerId()));
}

void PathTracker::fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) {
    drawFiguresForVertex(new_vertex, "MergePath_" + std::to_string(new_vertex.getInnerId()));
}

void PathTracker::fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) {
    drawFiguresForVertex(new_edge.getStart(), "MergePathToEdge_" + std::to_string(new_edge.getStart().getInnerId()));
    drawFiguresForVertex(new_edge.getFinish(), "MergePathToEdge_" + std::to_string(new_edge.getFinish().getInnerId()));
}

void PathTracker::fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                       const AlignmentForm &left_al, const AlignmentForm &right_al) {
    drawFiguresForVertex(new_edge.getStart(), "MergeTipsToEdge_" + std::to_string(new_edge.getStart().getInnerId()));
    drawFiguresForVertex(new_edge.getFinish(), "MergeTipsToEdge_" + std::to_string(new_edge.getFinish().getInnerId()));
}

void PathTracker::fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) {
    drawFiguresForVertex(split.getStart(), "SplitEdge_" + std::to_string(split.getStart().getInnerId()));
    drawFiguresForVertex(split.getFinish(), "SplitEdge_" + std::to_string(split.getFinish().getInnerId()));
}

void PathTracker::fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) {
    for (Vertex &new_vertex : resolution.newVertices())
        drawFiguresForVertex(new_vertex, "ResolveVertex_" + std::to_string(core.getInnerId()) + "_" +
                                 std::to_string(new_vertex.getInnerId()));
}
