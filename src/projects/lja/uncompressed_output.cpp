#include "uncompressed_output.hpp"

#include <sequences/edit_distance.hpp>
#include <common/logging.hpp>
#include <common/omp_utils.hpp>
#include <ksw2/ksw_wrapper.hpp>
#include "multi_graph.hpp"

struct OverlapRecord {
    OverlapRecord(multigraph::Edge *left, multigraph::Edge *right, std::vector<cigar_pair> cigar) : left(left), right(right),
                                                                                                    cigar(std::move(cigar)) {
        if(!cigar.empty() && cigar.back().type == 'D')
            cigar.pop_back();
        if(!cigar.empty() && cigar.front().type == 'I')
            cigar = {cigar.begin() + 1, cigar.end()};
    }

    multigraph::Edge *left;
    multigraph::Edge *right;
    std::vector<cigar_pair> cigar;
    std::string cigarString() const {
        std::stringstream ss;
        for(auto &pair : cigar) {
            if(pair.length > 1) {
                ss << pair.length;
            }
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
};

size_t compressedPrefixSize(const Sequence &hpcPrefix, const Sequence &seq) {
    StringContig sc(seq.str(), "tmp");
    Sequence hpcSeq = sc.makeSequence();
    auto projection = bestPrefix(hpcPrefix, hpcSeq);
    size_t hpc_size = projection.first;
    size_t pos = 0;
    size_t hpc_pos = 0;
    while(hpc_pos < hpc_size) {
        while(hpcSeq[hpc_pos] == seq[pos]) {
            pos++;
        }
        hpc_pos++;
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
    KSWAligner kswAligner(1, 10, 10, 1);
    return kswAligner.iterativeBandAlign(left_seq.str(), right_seq.str(), 5, 100, 0.01);
}

void printUncompressedResults(logging::Logger &logger, size_t threads, multigraph::MultiGraph &graph,
                              const std::vector<Contig> &uncompressed, const std::experimental::filesystem::path &out_dir) {
    std::unordered_map<int, Sequence> uncompression_results;
    for(const Contig &contig : uncompressed) {
        uncompression_results[std::stoi(contig.id)] = contig.seq;
        uncompression_results[-std::stoi(contig.id)] = !contig.seq;
    }
    ParallelRecordCollector<OverlapRecord> cigars_collection(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(graph, cigars_collection, uncompression_results)
    for(size_t i = 0; i < graph.vertices.size(); i++) {
        multigraph::Vertex &vertex = *graph.vertices[i];
        if(!vertex.isCanonical())
            continue;
        for (multigraph::Edge *out_edge : vertex.outgoing) {
            for (multigraph::Edge *inc_edge : vertex.rc->outgoing) {
                VERIFY(out_edge->getSeq().startsWith(vertex.seq));
                VERIFY(inc_edge->getSeq().startsWith(!vertex.seq));
                std::vector<cigar_pair> cigar = UncompressOverlap(graph.vertices[i]->seq, uncompression_results[inc_edge->rc->getId()],
                                                                  uncompression_results[out_edge->getId()]);
                cigars_collection.emplace_back(inc_edge, out_edge->rc, cigar);
            }
        }
    }
    std::ofstream os;
    os.open(out_dir / "mdbg.gfa");
    os << "H\tVN:Z:1.0" << std::endl;
    std::unordered_map<multigraph::Edge *, std::string> eids;
    for(multigraph::Edge *edge : graph.edges){
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
    os_cut.open(out_dir / "assembly.fasta");
    std::unordered_map<multigraph::Vertex *, size_t> cut;
    for(multigraph::Vertex *v : graph.vertices) {
        if(v->seq <= !v->seq) {
            if(v->outDeg() == 1) {
                cut[v] = 0;
            } else {
                cut[v] = 1;
            }
            cut[v->rc] = 1 - cut[v];
        }
    }
    std::unordered_map<multigraph::Edge*, size_t> cuts;
    for(OverlapRecord &rec : cigars_collection) {
        cuts[rec.left->rc] *= rec.endSize();
        cuts[rec.right] = rec.startSize();
    }

    std::vector<Contig> res;
    size_t cnt = 1;
    for(multigraph::Edge *edge : graph.edges) {
        if(edge->isCanonical()) {
            size_t cut_left = 0;
            if(edge->start->inDeg() != 0) {
                cut_left = cuts[edge];
            }
            cut_left *= cut[edge->start];
            size_t cut_right = 0;
            if(edge->end->outDeg() != 0) {
                cut_right = cuts[edge->rc];
            }
            cut_right *= (1 - cut[edge->end]);
            if(cut_left + cut_right >= edge->size()) {
                continue;
            }
            res.emplace_back(edge->getSeq().Subseq(cut_left, edge->size() - cut_right), itos(cnt));
            cnt++;
        }
    }
}