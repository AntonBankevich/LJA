#include "polishing_stage.hpp"
#include <common/cl_parser.hpp>
#include <assembly_graph/visualization.hpp>
#include "dbg/aln_reads_reader.hpp"

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

using namespace ag;

Segment<Vertex> ExtendToSize(Segment<Vertex> seg, size_t min_size) {
    if(seg.size() < min_size) {
        size_t left_ext = seg.cutLeft();
        size_t right_ext = seg.cutRight();
        if(left_ext *2 < min_size + 1 - seg.size()) {
            right_ext = min_size - left_ext;
        } else if(left_ext *2 < min_size + 1 - seg.size()) {
            left_ext = min_size - right_ext;
        } else {
            left_ext = (min_size + 1 -seg.size()) / 2;
            right_ext = (min_size + 1 - seg.size()) / 2;
        }
        VERIFY(seg.size() + left_ext + right_ext >= min_size);
        return {seg.contig(), seg.left - std::min(seg.left, left_ext), std::min(seg.right + right_ext, seg.contig().size())};
    } else {
        return seg;
    }
}
std::unordered_map<VertexId, Segment<Vertex>> ConstructReduction(AssemblyGraph &graph, size_t min_size, size_t max_overlap) {
    std::unordered_map<VertexId, Segment<Vertex>> reduction;
    std::vector<VertexId> list = oneline::map(graph.verticesUnique().begin(), graph.verticesUnique().end(),
                                              IdTransformer<Vertex>());
    std::sort(list.begin(), list.end(), [](const VertexId &vid1, const VertexId &vid2) {return vid1->size() > vid2->size() || (vid1->size() == vid2->size() && vid1 > vid2);});
    VERIFY(list.empty() || list.front()->size() >= list.back()->size());
    for(Vertex &vertex : graph.vertices()) {reduction[vertex.getId()] = Segment(vertex);}
    for(VertexId vid : list) {
        if(vid->isCore()) continue;
        size_t left_cut = vid->size();
        size_t right_cut = vid->size();
        for(Edge &edge : *vid) {
            if(edge.isPrefix()) {
                Segment<Vertex> seg = reduction.at(edge.getFinish().getId());
                left_cut = std::min(left_cut, seg.left);
            } else {
                VERIFY(edge.isSuffix());
                right_cut = std::min(right_cut, edge.getFinish().size());
            }
        }
        for(Edge &edge : vid->rc()) {
            if(edge.isPrefix()) {
                Segment<Vertex> seg = reduction.at(edge.getFinish().getId());
                right_cut = std::min(right_cut, seg.left);
            } else {
                VERIFY(edge.isSuffix());
                left_cut = std::min(left_cut, edge.getFinish().size());
            }
        }
        if(left_cut == vid->size()) left_cut = 0;
        if(right_cut == vid->size()) right_cut = 0;
        size_t left = left_cut;
        size_t right = vid->size() - right_cut;
        VERIFY(left < right || vid->isOuter());
        Segment<Vertex> res(*vid, std::min(left, right), std::max(left, right));
        res = ExtendToSize(res, min_size);
        reduction[vid] = res;
        reduction[vid->rc().getId()] = res.RC();
    }
//    for(Vertex &v : graph.verticesUnique()) {
//        if(v.isCore() && v.size() < max_overlap) {
//            size_t right = 0;
//            for(Edge &edge : v.incoming()) {
//                right = std::max(reduction.at(edge.getStart().getId()).right, right);
//            }
//            if(right > 0)
//                right = std::min(v.size(), std::max(right, min_size));
//            reduction[v.getId()] = {v, 0, right};
//            reduction[v.rc().getId()] = reduction[v.getId()].RC();
//            for(Edge &edge : v) {
//                reduction[edge.getFinish().getId()].left = 0;
//                reduction[edge.getFinish().rc().getId()].right = edge.getFinish().size();
//            }
//        }
//    }
    return std::move(reduction);
}

std::vector<Segment<Vertex>> ExtendSegments(std::vector<Segment<Vertex>> &segs, size_t min_alignment) {
    std::vector<Segment<Vertex>> extended_segs;
    for(const Segment<Vertex> &seg : segs) {
        extended_segs.emplace_back(ExtendToSize(seg, min_alignment));
    }
    return extended_segs;
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

bool IsRedundant(Vertex &vertex, const std::unordered_map<VertexId, Segment<Vertex>> &segs) {
    size_t left = segs.at(vertex.getId()).left;
    for(Edge &edge : vertex) {
        if(edge.isPrefix() && segs.at(edge.getFinish().getId()).left <= left)
            return true;
    }
    return false;
}

std::unordered_map <std::string, std::experimental::filesystem::path>
RunPolishing(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
             const std::experimental::filesystem::path &gfa_file,
             const io::Library &corrected_reads, const io::Library &reads, size_t min_alignment, bool debug) {
    logger.info() << "Performing polishing and homopolymer uncompression" << std::endl;
    size_t dicompress = StringContig::max_dimer_size / 2;
    dbg::SeqReader reader(corrected_reads, logger, threads);
    ag::AssemblyGraph graph = LoadSupregraphFromGFA(logger, threads, gfa_file);
    std::unordered_map<VertexId, Segment<Vertex>> segs = ConstructReduction(graph, min_alignment, 20000);
    std::vector<Segment<Vertex>> extra;
    for(Vertex &vertex : graph.verticesUnique()) {
        if(!CheckMergeRight(vertex, segs)) {
            Segment<Vertex> seg = segs.at(vertex.getId());
            if(ExtendLeft(vertex, segs)) seg.left = 0;
            if(ExtendLeft(vertex.rc(), segs)) seg.right = vertex.size();
            extra.emplace_back(seg);
        }
    }
    for(Segment<Vertex> seg : extra) {
        segs[seg.contig().getId()] = seg;
        segs[seg.contig().rc().getId()] = seg.RC();
    }
    std::function<Contig(Vertex &)> seg_to_contig = [&segs](Vertex &v) {
        return Contig(segs.at(v.getId()).fullSeq(), std::to_string(segs.at(v.getId()).contig().getInnerId()));
    };
    std::vector<Contig> compressed = oneline::map(graph.verticesUnique().begin(), graph.verticesUnique().end(), seg_to_contig);
    auto res = PrintAlignments(logger, threads, compressed, reader.begin(), reader.end(), min_alignment, dir);
    std::vector<Contig> uncompressed = Polish(logger, threads, compressed, res.first, reads, dicompress);
    std::unordered_map<VertexId , Sequence> uncompression_results;
    IdIndex<Vertex> index(graph.vertices().begin(), graph.vertices().end());
    for(const Contig &contig : uncompressed) {
        VertexId vid = index.getById(Parse<Vertex::id_type>(contig.getInnerId())).getId();
        uncompression_results[vid] = contig.getSeq();
        uncompression_results[vid->rc().getId()] = !contig.getSeq();
    }
    printUncompressedResults(logger, threads, graph, segs, uncompression_results, dir, debug);
    std::vector<Contig> assembly;
    for(Vertex &vertex : graph.verticesUnique()) {
        if(!IsRedundant(vertex, segs))
            assembly.emplace_back(uncompression_results.at(vertex.getId()), "LJA_Contig_" + itos(vertex.getInnerId()));
    }
    std::sort(assembly.begin(), assembly.end(), [](const Contig &c1, const Contig &c2)->bool{return c1.fullSize() > c2.fullSize();});
    PrintAssemblyStatistics(logger, assembly);
    std::ofstream os_cut(dir / "assembly.fasta");
    for(Contig &contig : assembly) {
        os_cut << ">" << contig.getInnerId() << "\n" << contig.getSeq() << "\n";
    }
    os_cut.close();
    logger.info() << "Polished assembly results can be found in: " << (dir / "assembly.fasta") << std::endl;
    return {{"assembly", dir / "assembly.fasta"}, {"graph", dir / "mdbg.gfa"}};
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
