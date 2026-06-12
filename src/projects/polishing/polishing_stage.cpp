#include "polishing_stage.hpp"
#include <common/cl_parser.hpp>
#include <assembly_graph/visualization.hpp>

#include "vertex_reduction.hpp"
#include "dbg/aln_reads_reader.hpp"

using namespace ag;

size_t Nx(const std::vector<size_t> &lens, size_t perc) {
    VERIFY(!lens.empty());
    size_t total = std::accumulate(lens.begin(), lens.end(), size_t(0));
    size_t pref_sum = 0;
    for(size_t len : lens) {
        pref_sum += len;
        if(pref_sum * 100 >= total * perc)
            return len;
    }
    return lens.back();
}

void PrintAssemblyStatistics(logging::Logger &logger, const std::vector<Contig> &contigs) {
    std::vector<size_t> lens;
    for(const Contig &contig : contigs) lens.emplace_back(contig.fullSize());
    std::sort(lens.begin(), lens.end(), std::greater<>());
    logger.info() << "Total contig length: " << std::accumulate(lens.begin(), lens.end(), size_t(0)) << std::endl;
    logger.info() << "Number of contigs: " << lens.size() << std::endl;
    if(lens.empty())
        return;
    logger.info() << "N50: " << Nx(lens, 50) << " N90: " << Nx(lens, 90) << std::endl;
}

bool CheckMergeRight(Vertex &vertex, const std::unordered_map<VertexId, Segment<Vertex>> &segs) {
    if(vertex.outDeg() == 0 || !vertex.isCore() || !vertex.isCanonical() || vertex.size() > 40000)
        return false;
    if(segs.at(vertex.getId()) != Segment(vertex))
        return false;
    for(Edge &edge : vertex) {
        Vertex &next = edge.getFinish();
        if(segs.at(next.getId()).left != vertex.size())
            return false;
    }
    return true;
}

bool ExtendLeft(Vertex &vertex, const std::unordered_map<VertexId, Segment<Vertex>> &segs) {
    if(vertex.inDeg() != 1 || !vertex.rc().front().isSuffix())
        return false;
    Vertex &prev = vertex.rc().front().getFinish().rc();
    return CheckMergeRight(prev, segs);
}

std::unordered_map <std::string, std::experimental::filesystem::path>
RunPolishing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
             const std::experimental::filesystem::path &gfa_file,
             const io::Library &corrected_reads, const io::Library &reads, size_t min_alignment, bool debug) {
    logger.info() << "Performing polishing and homopolymer uncompression" << std::endl;
    ag::AssemblyGraph graph = LoadSupregraphFromGFA(logger, threads, gfa_file);
    // std::vector<Segment<Vertex>> extra;
    // for(Vertex &vertex : graph.verticesUnique()) {
    //     if(!CheckMergeRight(vertex, segs)) {
    //         Segment<Vertex> seg = segs.at(vertex.getId());
    //         if(ExtendLeft(vertex, segs)) seg.left = 0;
    //         if(ExtendLeft(vertex.rc(), segs)) seg.right = vertex.size();
    //         extra.emplace_back(seg);
    //     }
    // }
    // for(Segment<Vertex> seg : extra) {
    //     segs[seg.contig().getId()] = seg;
    //     segs[seg.contig().rc().getId()] = seg.RC();
    // }
    DecompressingManager manager(graph, dir, min_alignment, 40000, debug);
    manager.ReduceAndUncompress(logger, threads, corrected_reads, reads);
    manager.calculateOverlaps(logger, threads);
    manager.printUncompressedGraph(logger, threads, dir/"final_graph.gfa");
    std::vector<Contig> assembly = manager.printAssembly(logger, threads);
    PrintAssemblyStatistics(logger, assembly);
    std::ofstream os_cut(dir / "assembly.fasta");
    for(Contig &contig : assembly) {
        os_cut << ">" << contig.getInnerId() << "\n" << contig.getSeq() << "\n";
    }
    os_cut.close();
    logger.info() << "Polished assembly results can be found in: " << (dir / "assembly.fasta") << std::endl;

    return {{"assembly", dir / "assembly.fasta"}, {"graph", dir / "final_graph.gfa"}};
}

std::unordered_map<std::string, std::experimental::filesystem::path>
PolishingPhase::innerRun(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                         bool debug, const AlgorithmParameterValues &parameterValues,
                         const std::unordered_map<std::string, io::Library> &input) {
    logger.info() << "Started homopolymer uncompression and polishing phase\n";
    size_t min_alignment = std::stoull(parameterValues.getValue("min-alignment"));
    return RunPolishing(logger, threads, dir, input.find("graph")->second.front(),
                        input.find("corrected_reads")->second, input.find("reads")->second, min_alignment, debug);
}
