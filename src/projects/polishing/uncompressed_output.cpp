#include "uncompressed_output.hpp"

#include <sequences/edit_distance.hpp>
#include <common/logging.hpp>
#include <common/omp_utils.hpp>
#include <alignment/ksw_wrapper.hpp>
#include <utility>
#include <assembly_graph/visualization.hpp>

#include "homopolish.hpp"
#include "perfect_alignment.hpp"
#include "vertex_reduction.hpp"
#include "dbg/multi_graph.hpp"
#include "alignment/ksw_aligner.hpp"
#include "dbg/aln_reads_reader.hpp"

using namespace ag;

struct OverlapRecord {
    OverlapRecord(Edge &edge, Sequence seq_left, Sequence seq_right, AlignmentForm _cigar) :
                            edge(edge.getId()), seq_left(std::move(seq_left)), seq_right(std::move(seq_right)),
                                                                                                    cigar(std::move(_cigar)) {
        if(!cigar.empty() && cigar.back().type == 'I')
            cigar.pop_back();
        if(!cigar.empty() && cigar.front().type == 'D')
            cigar.pop_front();
    }

    EdgeId edge;
    Sequence seq_left, seq_right;
    AlignmentForm cigar;
    std::string cigarString() const {
        std::stringstream ss;
        for(auto &pair : cigar) {
            ss << pair.length;
            ss << pair.type;
        }
        if(cigar.empty())
            ss<< "0M";
        return ss.str();
    }

    size_t endSize() const {
        size_t res = 0;
        for(auto &pair : cigar) {
            if(pair.type != 'I')
                res += pair.length;
        }
        return res;
    }
    size_t startSize() const {
        size_t res = 0;
        for(auto &pair : cigar) {
            if(pair.type != 'D')
                res += pair.length;
        }
        return res;
    }
    std::pair<std::string, std::string> str() const {
        size_t pos_ref = seq_left.size() - endSize();
        size_t pos_q = 0;
        std::vector<char> s_ref;
        std::vector<char> s_q;
        for(const CigarPair &pair : cigar) {
            for(size_t i = 0; i < pair.length; i++) {
                if(pair.type == 'I')
                    s_ref.emplace_back('-');
                else
                    s_ref.emplace_back("ACGT"[seq_left[pos_ref + i]]);
                if(pair.type == 'D')
                    s_q.emplace_back('-');
                else
                    s_q.emplace_back("ACGT"[seq_right[pos_q + i]]);
            }
            if(pair.type != 'I')
                pos_ref += pair.length;
            if(pair.type != 'D')
                pos_q += pair.length;
        }
        return {std::string(s_q.begin(), s_q.end()), std::string(s_ref.begin(), s_ref.end())};
    }
};

size_t compressedPrefixSize(const Sequence &hpcPrefix, const Sequence &seq) {
    StringContig sc(seq.str(), "tmp");
    Sequence hpcSeq = sc.makeSequence();
    auto projection = bestPrefix(hpcPrefix, hpcSeq);
    size_t hpc_size = projection.first;
    size_t pos = 0;
    size_t hpc_pos = 0;
    while(hpc_pos < hpc_size) {
        while(pos + 1 < seq.size() && hpcSeq[hpc_pos] == seq[pos]) {
            pos++;
        }
        hpc_pos++;
        if(hpc_pos == hpc_size)
            break;
        if(pos < seq.size() && hpc_pos < hpc_size && hpcSeq[hpc_pos] != seq[pos]) {
            VERIFY(hpc_pos >= 2 && seq[pos] == hpcSeq[hpc_pos - 2]);
            hpc_pos -= 2;

        } else {
            VERIFY(hpcSeq[hpc_pos] == seq[pos]);
        }
    }
    return pos;
}

size_t homoSize(const Sequence &s, size_t pos) {
    while(pos > 0 && s[pos] == s[pos - 1])
        pos--;
    size_t res = 0;
    while(pos + 1 < s.size() && s[pos] == s[pos + 1]) {
        res += 1;
        pos += 1;
    }
    return res;
}
size_t leftHomoSize(const Sequence &s) {
    if(s.empty())
        return 0;
    return homoSize(s, 0);
}

size_t rightHomoSize(const Sequence &s) {
    if(s.empty())
        return 0;
    return homoSize(s, s.size() - 1);
}

AlignmentForm UncompressOverlap(const Sequence &hpcOverlap, const Sequence &left, const Sequence & right) {
    StringContig sc(left.str(), "left");
    size_t left_len = compressedPrefixSize(!hpcOverlap, !left);
    size_t right_len = compressedPrefixSize(hpcOverlap, right);
    Sequence left_seq = left.Subseq(left.size() - left_len);
    Sequence right_seq = right.Subseq(0, right_len);
    if(leftHomoSize(left_seq) > leftHomoSize(right_seq)) {
        left_seq = left_seq.Subseq(leftHomoSize(left_seq) - leftHomoSize(right_seq));
    }
    if(rightHomoSize(left_seq) < rightHomoSize(right_seq)) {
        right_seq = right_seq.Subseq(0, right_seq.size() - (rightHomoSize(right_seq) - rightHomoSize(left_seq)));
    }
    KSWAligner kswAligner(1, 5, 10, 2);
    AlignmentForm res;
    return kswAligner.globalAlignment(left_seq.str(), right_seq.str());
}

void DecompressingManager::ReduceAndUncompress(logging::Logger &logger, size_t threads,
        const io::Library &corrected_reads, const io::Library &reads) {
    size_t dicompress = StringContig::max_dimer_size / 2;
    dbg::SeqReader reader(corrected_reads, logger, threads);
    segs = ConstructReduction(graph(), min_overlap, max_repeat);
    std::function<Contig(Vertex &)> seg_to_contig = [this](Vertex &v) {
        return Contig(segs.at(v.getId()).fullSeq(), std::to_string(segs.at(v.getId()).contig().getInnerId()));
    };
    std::vector<Contig> compressed = oneline::map(graph().verticesUnique().begin(), graph().verticesUnique().end(), seg_to_contig);
    auto res = PrintAlignments(logger, threads, compressed, reader.begin(), reader.end(), min_overlap, dir);
    std::vector<Contig> uncompressed_contigs = Polish(logger, threads, compressed, res.first, reads, dicompress);
    IdIndex<Vertex> index(graph().vertices().begin(), graph().vertices().end());
    for(const Contig &contig : uncompressed_contigs) {
        VertexId vid = index.getById(Parse<Vertex::id_type>(contig.getInnerId())).getId();
        uncompressed[vid] = contig.getSeq();
        uncompressed[vid->rc().getId()] = !contig.getSeq();
    }
}

void DecompressingManager::printReduction(std::experimental::filesystem::path path) const {
    Printer printer(VertexPrintStyles::defaultDotInfo(), EdgePrintStyles::defaultDotInfo());
    const Vertex &v = *(graph().vertices().begin());
    std::function<std::string(const Vertex &)> vertex_labels = [this](const Vertex &v) {return std::to_string(segs.at(v.getId()).left) +" " + std::to_string(segs.at(v.getId()).right);};
    printer += VertexInfo::Labeler(vertex_labels);
    printer.printDot(path, graph());
}


void DecompressingManager::calculateOverlaps(logging::Logger &logger, size_t threads) {
    logger.info() << "Calculating overlaps between adjacent uncompressed edges" << std::endl;
    omp_set_num_threads(threads);
    std::vector<multigraph::EdgeId> e_ids;
    for (Edge &e: graph().edgesUnique())
        e_ids.push_back(e.getId());
#pragma omp parallel for default(none) shared(e_ids, segs)
    for(size_t i = 0; i < e_ids.size(); i++) {
        Edge &edge = *e_ids[i];
        Segment<Vertex> left_seg = segs.at(edge.getStart().getId());
        Segment<Vertex> right_seg = segs.at(edge.getFinish().getId());
        size_t shift = edge.rc().truncSize();
        Sequence overlap = edge.getStart().getSeq().Subseq(std::max(left_seg.left, right_seg.left + shift),
                                                           std::min(left_seg.right, right_seg.right + shift));
        AlignmentForm cigar = UncompressOverlap(overlap, uncompressed.at(edge.getStart().getId()),
                                                         uncompressed.at(edge.getFinish().getId()));
        OverlapRecord overlapRecord(edge, uncompressed.at(edge.getStart().getId()),
                                    uncompressed.at(edge.getFinish().getId()), cigar);
#pragma omp critical
        overlap_alignment[edge.getId()] = cigar;
#pragma omp critical
        overlap_alignment[edge.rc().getId()] = cigar.Reverse().RC();
    }
}

void DecompressingManager::printUncompressedGraph(logging::Logger &logger, size_t threads, std::experimental::filesystem::path path) {
    logger.info() << "Printing polished gfa file to " << (path) << std::endl;
    std::ofstream os;
    os.open(path);
    os << "H\tVN:Z:1.0" << std::endl;
    std::unordered_map<multigraph::Edge *, std::string> eids;
    for(Vertex &vertex : graph().verticesUnique()){
        os << "S\t" << vertex.getId() << "\t" << uncompressed.at(vertex.getId()) << "\n";
    }
    for(Edge &edge : graph().edgesUnique()) {
        bool inc_sign = edge.getStart().isCanonical();
        VertexId incId = inc_sign ? edge.getStart().getId() : edge.getStart().rc().getId();
        bool out_sign = edge.getFinish().isCanonical();
        VertexId outId = out_sign ? edge.getFinish().getId() : edge.getFinish().rc().getId();
        os << "L\t" << incId << "\t" << (inc_sign ? "+" : "-") << "\t" << outId << "\t"
           << (out_sign ? "+" : "-") << "\t" << overlap_alignment[edge.getId()].toCigarString() << "\n";
    }
    os.close();
}

std::vector<Contig> DecompressingManager::printAssembly(logging::Logger &logger, size_t threads) {
    std::vector<Contig> assembly;
    std::unordered_map<VertexId, size_t> cut;
    for (Vertex &vertex : graph().vertices()) {cut[vertex.getId()] = 0;}
    for(Vertex &vertex : graph().vertices()) {
        if (vertex.outDeg()> 1 && segs[vertex.getId()].right == vertex.size())
            for (Edge &edge : vertex) {
                cut[edge.getFinish().getId()] = overlap_alignment[edge.getId()].targetLength();
            }
    }
    for(Vertex &vertex : graph().verticesUnique()) {
        Sequence seq = uncompressed.at(vertex.getId());
        size_t left = cut.at(vertex.getId());
        size_t right = cut.at(vertex.rc().getId());
        if (seq.size() > left + right + 5000 && !(vertex.isCore() && segs[vertex.getId()].size() == min_overlap)) {
            assembly.emplace_back(seq.Subseq(left, seq.size() - right), itos(vertex.getInnerId()));

        }
    }
    std::sort(assembly.begin(), assembly.end(), [](const Contig &c1, const Contig &c2)->bool{return c1.fullSize() > c2.fullSize();});
    return std::move(assembly);
}

// void printUncompressedResults(logging::Logger &logger, size_t threads, multigraph::MultiGraph &graph,
//                               const std::unordered_map<VertexId, Segment<Vertex>> &segs,
//                               const std::unordered_map<VertexId , Sequence> &uncompression_results,
//                               const std::experimental::filesystem::path &out_dir, bool debug) {
//     logger.info() << "Printing polished gfa file to " << (out_dir / "mdbg.gfa") << std::endl;
//     std::ofstream os;
//     os.open(out_dir / "mdbg.gfa");
//     os << "H\tVN:Z:1.0" << std::endl;
//     std::unordered_map<multigraph::Edge *, std::string> eids;
//     for(Vertex &vertex : graph.verticesUnique()){
//         os << "S\t" << vertex.getId() << "\t" << uncompression_results.at(vertex.getId()) << "\n";
//     }
//     for(OverlapRecord &rec : cigars_collection) {
//         bool inc_sign = rec.edge->getStart().isCanonical();
//         VertexId incId = inc_sign ? rec.edge->getStart().getId() : rec.edge->getStart().rc().getId();
//         bool out_sign = rec.edge->getFinish().isCanonical();
//         VertexId outId = out_sign ? rec.edge->getFinish().getId() : rec.edge->getFinish().rc().getId();
//         os << "L\t" << incId << "\t" << (inc_sign ? "+" : "-") << "\t" << outId << "\t"
//            << (out_sign ? "+" : "-") << "\t" << rec.cigarString() << "\n";
//     }
//     os.close();
// }
