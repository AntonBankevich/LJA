#include "uncompressed_output.hpp"

#include <sequences/edit_distance.hpp>
#include <common/logging.hpp>
#include <common/omp_utils.hpp>
#include <ksw2/ksw_wrapper.hpp>
#include "multi_graph.hpp"

struct OverlapRecord {
    OverlapRecord(multigraph::Edge *left, multigraph::Edge *right, const Sequence &seq_left, const Sequence &seq_right, std::vector<cigar_pair> _cigar) : left(left), right(right), seq_left(seq_left), seq_right(seq_right),
                                                                                                    cigar(std::move(_cigar)) {
        if(!cigar.empty() && cigar.back().type == 'I')
            cigar.pop_back();
        if(!cigar.empty() && cigar.front().type == 'D')
            cigar = {cigar.begin() + 1, cigar.end()};
    }

    multigraph::Edge *left;
    multigraph::Edge *right;
    Sequence seq_left, seq_right;
    std::vector<cigar_pair> cigar;
    std::string cigarString() const {
        std::stringstream ss;
        for(auto &pair : cigar) {
            ss << pair.length;
            ss << pair.type;
        }
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
        for(const cigar_pair &pair : cigar) {
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

std::vector<cigar_pair> UncompressOverlap(const Sequence &hpcOverlap, const Sequence &left, const Sequence & right) {
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
    return kswAligner.iterativeBandAlign(left_seq.str(), right_seq.str(), 5, 100, 0.01);
}

std::vector<Contig> printUncompressedResults(logging::Logger &logger, size_t threads, multigraph::MultiGraph &graph,
                              const std::vector<Contig> &uncompressed, const std::experimental::filesystem::path &out_dir, bool debug) {
    logger.info() << "Calculating overlaps between adjacent uncompressed edges" << std::endl;
    std::unordered_map<int, Sequence> uncompression_results;
    for(const Contig &contig : uncompressed) {
        uncompression_results[std::stoi(contig.id)] = contig.seq;
        uncompression_results[-std::stoi(contig.id)] = !contig.seq;
    }
    ParallelRecordCollector<OverlapRecord> cigars_collection(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(graph, cigars_collection, uncompression_results, logger, debug, std::cout)
    for(size_t i = 0; i < graph.vertices.size(); i++) {
        multigraph::Vertex &vertex = *graph.vertices[i];
        if(!vertex.isCanonical())
            continue;
        for (multigraph::Edge *out_edge : vertex.outgoing) {
            for (multigraph::Edge *inc_edge : vertex.rc->outgoing) {
                VERIFY_OMP(out_edge->getSeq().startsWith(vertex.seq));
                VERIFY_OMP(inc_edge->getSeq().startsWith(!vertex.seq));
                std::vector<cigar_pair> cigar = UncompressOverlap(graph.vertices[i]->seq, uncompression_results[inc_edge->rc->getId()],
                                                                  uncompression_results[out_edge->getId()]);
                OverlapRecord overlapRecord(inc_edge->rc, out_edge, uncompression_results[inc_edge->rc->getId()],
                                                                  uncompression_results[out_edge->getId()], cigar);
                if(debug) {
#pragma omp critical
                    {
                        logger.debug() << inc_edge->getId() << " " << std::endl << out_edge->getId() << " " <<overlapRecord.cigarString() << " " << overlapRecord.startSize() << " " << overlapRecord.endSize() << std::endl;
                        std::pair<std::string, std::string> al = overlapRecord.str();
                        logger.debug() << al.first << "\n" << al.second << std::endl;
                    }
                }
                cigars_collection.emplace_back(std::move(overlapRecord));
            }
        }
    }
    logger.info() << "Printing final gfa file to " << (out_dir / "mdbg.gfa") << std::endl;
    std::ofstream os;
    os.open(out_dir / "mdbg.gfa");
    os << "H\tVN:Z:1.0" << std::endl;
    std::unordered_map<multigraph::Edge *, std::string> eids;
    for(auto p : graph.edges){
        multigraph::Edge *edge = p.second;
        if (edge->isCanonical()) {
            os << "S\t" << itos(edge->getId()) << "\t" << uncompression_results[edge->getId()] << "\n";
        }
    }
    for(OverlapRecord &rec : cigars_collection) {
        bool inc_sign = rec.left->isCanonical();
        int incId = inc_sign ? rec.left->getId() : rec.left->rc->getId();
        bool out_sign = rec.right->isCanonical();
        int outId = out_sign ? rec.right->getId() : rec.right->rc->getId();
        os << "L\t" << incId << "\t" << (inc_sign ? "+" : "-") << "\t" << outId << "\t"
           << (out_sign ? "+" : "-") << "\t" << rec.cigarString() << "\n";
    }
    os.close();
    std::ofstream os_cut;
    std::unordered_map<multigraph::Vertex *, size_t> cut; //Choice of vertex side for cutting
    for(auto p : graph.vertices) {
        multigraph::Vertex *v = p.second;
        if(v->seq <= !v->seq) {
            if(v->outDeg() == 1) {
                cut[v] = 0;
            } else {
                cut[v] = 1;
            }
            cut[v->rc] = 1 - cut[v];
        }
    }
    std::unordered_map<multigraph::Edge*, size_t> cuts; //Sizes of cuts from the edge start
    for(auto p : graph.edges) {
        multigraph::Edge *e = p.second;
        cuts[e] = 0;
    }
    for(OverlapRecord &rec : cigars_collection) {
        cuts[rec.left->rc] = cut[rec.left->rc->start] * rec.endSize();
        cuts[rec.right] = cut[rec.right->start] * rec.startSize();
    }
    std::vector<Contig> res;
    for(auto p : graph.edges) {
        multigraph::Edge *edge = p.second;
        if(edge->isCanonical()) {
            //TODO make canonical be the same as positive id
            size_t cut_left = cuts[edge];
            size_t cut_right = cuts[edge->rc];
            Sequence seq = uncompression_results[edge->getId()];
            if(cut_left + cut_right >= seq.size()) {
                continue;
            }
            res.emplace_back(seq.Subseq(cut_left, seq.size() - cut_right), itos(edge->getId()));
        }
    }
    return std::move(res);
}
