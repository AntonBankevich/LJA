#include "uncompressed_output.hpp"

#include <sequences/edit_distance.hpp>
#include <common/logging.hpp>
#include <common/omp_utils.hpp>
#include <ksw2/ksw_wrapper.hpp>
#include <utility>
#include "dbg/multi_graph.hpp"
using namespace multigraph;

struct OverlapRecord {
    OverlapRecord(multigraph::Edge &left, multigraph::Edge &right, Sequence seq_left, Sequence seq_right, std::vector<cigar_pair> _cigar) :
                            left(left.getId()), right(right.getId()), seq_left(std::move(seq_left)), seq_right(std::move(seq_right)),
                                                                                                    cigar(std::move(_cigar)) {
        if(!cigar.empty() && cigar.back().type == 'I')
            cigar.pop_back();
        if(!cigar.empty() && cigar.front().type == 'D')
            cigar = {cigar.begin() + 1, cigar.end()};
    }

    EdgeId left;
    EdgeId right;
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
    return kswAligner.iterativeBandAlign(left_seq.str(), right_seq.str(), 5, 100, 0.01, 0);
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
    std::vector<multigraph::VertexId> v_ids;
    for (Vertex &v: graph.vertices())
        if(v.isCanonical())
            v_ids.push_back(v.getId());
#pragma omp parallel for default(none) shared(graph, v_ids, cigars_collection, uncompression_results, logger, debug, std::cout)
    for(size_t i = 0; i < v_ids.size(); i++) {
        Vertex &vertex = *v_ids[i];
        for (Edge &out_edge : vertex) {
            for (Edge &inc_edge : vertex.rc()) {
                VERIFY_OMP(out_edge.getSeq().startsWith(vertex.getSeq()));
                VERIFY_OMP(inc_edge.getSeq().startsWith(!vertex.getSeq()));
                std::vector<cigar_pair> cigar = UncompressOverlap(vertex.getSeq(), uncompression_results[inc_edge.rc().getId().innerId()],
                                                                  uncompression_results[out_edge.getId().innerId()]);
                OverlapRecord overlapRecord(inc_edge.rc(), out_edge, uncompression_results[inc_edge.rc().getId().innerId()],
                                                                  uncompression_results[out_edge.getId().innerId()], cigar);
                if(debug) {
#pragma omp critical
                    {
                        logger.debug() << inc_edge.getId() << " " << std::endl << out_edge.getId() << " "
                                << overlapRecord.cigarString() << " " << overlapRecord.startSize() << " "
                                << overlapRecord.endSize() << std::endl;
                        std::pair<std::string, std::string> al = overlapRecord.str();
                        logger.debug() << al.first << "\n" << al.second << std::endl;
                    }
                }
                cigars_collection.emplace_back(std::move(overlapRecord));
            }
        }
    }
    logger.info() << "Printing polished gfa file to " << (out_dir / "mdbg.gfa") << std::endl;
    std::ofstream os;
    os.open(out_dir / "mdbg.gfa");
    os << "H\tVN:Z:1.0" << std::endl;
    std::unordered_map<multigraph::Edge *, std::string> eids;
    for(Edge &edge : graph.edges()){
        if (edge.isCanonical()) {
            os << "S\t" << edge.getId() << "\t" << uncompression_results[edge.getId().innerId()] << "\n";
        }
    }
    for(OverlapRecord &rec : cigars_collection) {
        bool inc_sign = rec.left->isCanonical();
        EdgeId incId = inc_sign ? rec.left->getId() : rec.left->rc().getId();
        bool out_sign = rec.right->isCanonical();
        EdgeId outId = out_sign ? rec.right->getId() : rec.right->rc().getId();
        os << "L\t" << incId << "\t" << (inc_sign ? "+" : "-") << "\t" << outId << "\t"
           << (out_sign ? "+" : "-") << "\t" << rec.cigarString() << "\n";
    }
    os.close();
    std::ofstream os_cut;
    std::unordered_map<VertexId, size_t> cut; //Choice of vertex side for cutting
    for(Vertex &v : graph.vertices()) {
        if(v.isCanonical()) {
            if(v.outDeg() == 1) {
                cut[v.getId()] = 0;
            } else {
                cut[v.getId()] = 1;
            }
            cut[v.rc().getId()] = 1 - cut[v.getId()];
        }
    }
    std::unordered_map<EdgeId, size_t> cuts; //Sizes of cuts from the edge start
    for(Edge &e : graph.edges()) {
        cuts[e.getId()] = 0;
    }
    for(OverlapRecord &rec : cigars_collection) {
        cuts[rec.left->rc().getId()] = cut[rec.left->rc().start().getId()] * rec.endSize();
        cuts[rec.right] = cut[rec.right->start().getId()] * rec.startSize();
    }
    std::vector<Contig> res;
    for(Edge &edge : graph.edges()) {
        if(edge.isCanonical()) {
            //TODO make canonical be the same as positive id
            size_t cut_left = cuts[edge.getId()];
            size_t cut_right = cuts[edge.rc().getId()];
            Sequence seq = uncompression_results[edge.getId().innerId()];
            if(cut_left + cut_right >= seq.size()) {
                continue;
            }
            res.emplace_back(seq.Subseq(cut_left, seq.size() - cut_right), itos(edge.getId().innerId()));
        }
    }
    return std::move(res);
}
