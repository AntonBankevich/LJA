#include "ReadsAligner.h"
#include "graph_alignment.hpp"
#include "alignment_correction.hpp"
#include <assembly_graph/paths.hpp>
#include <dbg/multi_graph.hpp>
#include <common/id_index.hpp>
#include <common/logging.hpp>
#include <common/cl_parser.hpp>
#include <sequences/seqio.hpp>
#include <alignment/ksw_wrapper.hpp>
#include <utility>
#include <queue>
#include <alignment/ksw_aligner.hpp>

using namespace multigraph;

std::unordered_map<std::string, std::vector<nano::GraphContig>>
AlignOnt(logging::Logger &logger, const size_t threads, const std::experimental::filesystem::path &dir,
         const MultiGraph &mg, const io::Library &ont_reads, bool reuse_alignment) {
    const unsigned int BATCH_SIZE = 10000;
    io::SeqReader reader(ont_reads);
    const std::experimental::filesystem::path &input_gfa = dir / "input.gfa";
    MultiGraphHelper::printEdgeGFA(mg, input_gfa);
    logger.info() << "Data loaded" << std::endl;

    std::unordered_map<std::string, Contig> batch;
    nano::ReadsAlignerGA reads_aligner(mg);

    int batch_num = 0;
    StringContig::homopolymer_compressing = true;
    StringContig::min_dimer_to_compress = 32;
    StringContig::max_dimer_size = 32;
    for(StringContig scontig : reader) {
        Contig contig = scontig.makeContig();
        batch[contig.getInnerId()] = contig;
        if (batch.size() > BATCH_SIZE) {
            if(!reuse_alignment)
                reads_aligner.Align(batch, input_gfa, dir, threads, batch_num);
            batch_num++;
            batch.clear();
        }
    }
    if(!reuse_alignment)
        reads_aligner.Align(batch, input_gfa, dir, threads, batch_num);
    batch_num++;
    batch.clear();
    std::unordered_map<std::string, std::vector<nano::GraphContig>> result;
    for (int i = 0; i < batch_num; i++) {
        for(auto &it : reads_aligner.ExtractPaths(dir, i)) {
            result[it.first].insert(result[it.first].end(), it.second.begin(), it.second.end());
        }
    }
    return std::move(result);
}

size_t Score(const AlignmentForm &cigar, const Sequence &from_seq, const Sequence &to_seq,
             const std::vector<std::pair<size_t, size_t>> &to_ignore) {
    size_t from_pos = 0;
    size_t to_pos = 0;
    size_t diff = 0;
    size_t cur = 0;
    for(const CigarPair &cp: cigar) {
        if(cp.type == 'M') {
            for(size_t i = 0; i < cp.length; i++) {
                if(cur < to_ignore.size() && from_pos + i >= to_ignore[cur].second)
                    cur++;
                if((cur == to_ignore.size() || from_pos + i < to_ignore[cur].first) && to_seq[to_pos + i] != from_seq[from_pos + i])
                    diff++;
            }
            from_pos += cp.length;
            to_pos += cp.length;
        }
        if(cp.type == 'D') {
            if((cur == to_ignore.size() || from_pos < to_ignore[cur].first)) {
                diff += cp.length;
            }
            to_pos += cp.length;
        }
        if(cp.type == 'I') {
            for(size_t i = 0; i < cp.length; i++) {
                if(cur < to_ignore.size() && from_pos + i >= to_ignore[cur].second)
                    cur++;
                if((cur == to_ignore.size() || from_pos + i < to_ignore[cur].first))
                    diff++;
            }
            from_pos += cp.length;
        }
    }
    return diff;
}

AlignmentForm defaultAlignExtension(const Sequence &from_seq, const Sequence &to_seq) {
    KSWAligner kswAligner(1, 5, 5, 3);
    return kswAligner.extendAlignment(to_seq.str(), from_seq.str());
}

AlignmentForm defaultAlign(const Sequence &from_seq, const Sequence &to_seq, bool extension = false) {
    if(extension) {
        return defaultAlignExtension(from_seq, to_seq);
    } else {
        KSWAligner kswAligner(1, 5, 5, 3);
        return kswAligner.globalAlignment(to_seq.str(), from_seq.str());
    }
}

//size_t Score(const Sequence &from_seq, const Sequence &to_seq, const std::vector<std::pair<size_t, size_t>> &to_ignore, bool extension = false) {
//    std::vector<cigar_pair> cigar = defaultAlign(from_seq, to_seq, extension);
//    return Score(cigar, from_seq, to_seq, to_ignore);
//}


std::vector<std::pair<size_t, size_t>> Mark(const Sequence &seq) {
    size_t w = 20;
    size_t threshold = 16;
    size_t merge_dist = 15;
    if(seq.size() < w) {
        return {};
    }
    std::vector<std::pair<size_t, size_t>> res;
    std::vector<size_t> cnt(4, 0);

    for(size_t i = 0; i < w; i++) {
        cnt[seq[i]]++;
    }
    for(size_t i = w; i < seq.size(); i++) {
        if(cnt[0] + cnt[1] > threshold || cnt[0] + cnt[2] > threshold || cnt[0] + cnt[3] > threshold
             || cnt[2] + cnt[1] > threshold || cnt[3] + cnt[1] > threshold || cnt[2] + cnt[3] > threshold) {
            if(res.empty() || i - w > res.back().second + merge_dist) {
                res.emplace_back(i - w, i);
            } else {
                res.back() = {res.back().first, i};
            }
        }
        cnt[seq[i]]++;
        cnt[seq[i - w]]--;
    }
    return std::move(res);
}

std::pair<size_t, size_t> Score2(const Sequence &read, const Sequence &p1, const Sequence &p2) {
    std::vector<std::pair<size_t, size_t>> to_ignore = Mark(read);
    AlignmentForm cigar1 = defaultAlign(read, p1);
    AlignmentForm cigar2 = defaultAlign(read, p2);
    if(cigar1.queryLength() == 0 || cigar2.queryLength() == 0)
        return {0, 0};
    return {Score(cigar1, read, p1, to_ignore), Score(cigar2, read, p2, to_ignore)};
}

multigraph::GraphPath ContigToPath(const nano::GraphContig &al, MultiGraph &graph) {
    size_t skip_left = al.gStart;
    size_t skip_right = 0;
    size_t len = al.gEnd - al.gStart;
    SequenceBuilder sb;
    multigraph::GraphPath path;
    IdIndex<Edge> index(graph.edges().begin(), graph.edges().end());
    for(std::string code : al.path) {
        bool rc = (code.back() == '-');
        Edge::id_type eid = Parse<Edge::id_type>(code.substr(0, code.size() - 1));
        Edge &edge = rc ? index.getById(eid).rc() : index.getById(eid);
        if(path.empty()) {
            if(skip_left >= edge.fullSize() - edge.getFinish().size()) {
                skip_left -= edge.fullSize() - edge.getFinish().size();
            } else {
                path += edge;
                size_t to_add = edge.fullSize() - skip_left;
                if(len <= to_add) {
                    skip_right = to_add - len;
                    len = 0;
                    break;
                } else {
                    len -= to_add;
                }
            }
        } else {
            path += edge;
            size_t to_add = edge.fullSize() - edge.getStart().size();
            if (len <= to_add) {
                skip_right = to_add - len;
                len = 0;
                break;
            } else {
                len -= to_add;
            }
        }
    }
    VERIFY(path.empty() || len == 0);
    path.cutFront(skip_left);
    path.cutBack(skip_right);
    return std::move(path);
}

std::vector<std::pair<size_t, size_t>> Nails(const multigraph::MultiGraph &graph, const Sequence &from_seq,
                                             const multigraph::GraphPath &mpath, const AlignmentForm &cigar) {
    Sequence to_seq = mpath.Seq();
    std::vector<std::pair<size_t, size_t>> vertices;
    size_t skip_left = mpath.leftSkip();
    size_t cur_pos = 0;
    for(multigraph::MGEdge & edge: mpath.edges()) {
        VERIFY(skip_left < edge.fullSize());
        size_t to_add = edge.fullSize() - skip_left;
        cur_pos += to_add;
        vertices.emplace_back(cur_pos - edge.getFinish().size(), cur_pos);
        skip_left = edge.getFinish().size();
    }
    vertices.pop_back();
    std::vector<std::pair<size_t, size_t>> res(vertices.size(), std::make_pair<size_t, size_t>(0, 0));
    std::vector<size_t> best_len(vertices.size(), 0);
    size_t from_pos = 0;
    size_t to_pos = 0;
    size_t cur = 0;
    for(const CigarPair &cp: cigar) {
        if(cp.type == CigarEvent::M) {
            size_t match = 0;
            for(size_t i = 0; i < cp.length; i++) {
                if(to_seq[to_pos + i] == from_seq[from_pos + i])
                    match++;
            }
            size_t cur1 = cur;
            while(cur1 < vertices.size() && vertices[cur1].first < to_pos + cp.length) {
                size_t left = std::max(vertices[cur1].first, to_pos);
                size_t right = std::min(vertices[cur1].second, to_pos + cp.length);
                VERIFY(left < right);
                size_t p = (right - left) / 2;
                while (p < right - left && from_seq[left - to_pos + from_pos + p] != to_seq[left + p]) {
                    p++;
                }
                if(p == right - left) {
                    p = 0;
                    while (p < right - left && from_seq[left - to_pos + from_pos + p] != to_seq[left + p]) {
                        p++;
                    }
                }
                if(p < right - left && cp.length > best_len[cur1]) {
                    res[cur1] = {left - to_pos + from_pos + p, left + p - vertices[cur1].first};
                    best_len[cur1] = cp.length;
                }
                cur1++;
            }
            from_pos += cp.length;
            to_pos += cp.length;
        }
        if(cp.type == CigarEvent::D) {
            to_pos += cp.length;
        }
        if(cp.type == CigarEvent::I) {
            from_pos += cp.length;
        }
        while(cur < vertices.size() && vertices[cur].second <= to_pos) {
            cur++;
        }
    }
    VERIFY(res.size() == vertices.size());
    for(size_t i = 0; i + 1< res.size(); i++) {
        if(res[i+1].first < res[i].first) {
            res[i+1].second += res[i].first - res[i+1].first;
            res[i+1].first = res[i].first;
        }
    }
    return std::move(res);
}


//void AnalyseAndPrint(const Sequence &from_seq, const Sequence &to_seq, int end_bonus = 0) {
//    std::vector<cigar_pair> cigar = defaultAlign(from_seq, to_seq, end_bonus);
//    std::vector<std::pair<size_t, size_t>> to_ignore = Mark(from_seq);
//    std::cout << CigarToString(cigar) << std::endl;
//    std::cout << Score(cigar, from_seq, to_seq, to_ignore) << " " << Score(cigar, from_seq, to_seq, {}) << std::endl;
//    std::vector<string> res = cigarToAlignmentString(from_seq, to_seq, cigar, to_ignore);
//    std::cout << res[0] << "\n" << res[1] << "\n" << res[2] << std::endl;
//}

std::vector<std::string> mixAndShorten(const std::vector<std::string> &s1, const std::vector<std::string> &s2) {
    std::vector<char> ref1;
    std::vector<char> m1;
    std::vector<char> read;
    std::vector<char> m2;
    std::vector<char> ref2;
    size_t cur1 = 0;
    size_t cur2 = 0;
    size_t cnt = 0;
    std::string insert = "======";
    while(cur1 < s1[0].size() || cur2 < s2[0].size()) {
        if(cur1 < s1[0].size() && cur2 < s2[0].size() && s1[0][cur1] == s2[0][cur2]) {
            ref1.emplace_back(s1[2][cur1]);
            m1.emplace_back(s1[1][cur1]);
            read.emplace_back(s1[0][cur1]);
            m2.emplace_back(s2[1][cur2]);
            ref2.emplace_back(s2[2][cur2]);
            cur1++;
            cur2++;
            if(ref1.back() == ref2.back()) {
                cnt++;
            } else {
                cnt = 0;
            }
        } else if(cur1 < s1[0].size() && s1[0][cur1] == '-') {
            ref1.emplace_back(s1[2][cur1]);
            m1.emplace_back(s1[1][cur1]);
            read.emplace_back(s1[0][cur1]);
            if(s1[1][cur1] == '*')
                m2.emplace_back('*');
            else
                m2.emplace_back('|');
            ref2.emplace_back('-');
            cur1++;
            cnt = 0;
        } else if(cur2 < s2[0].size() && s2[0][cur2] == '-') {
            ref1.emplace_back('-');
            if(s2[1][cur2] == '*')
                m1.emplace_back('*');
            else
                m1.emplace_back('|');
            read.emplace_back(s2[0][cur2]);
            m2.emplace_back(s2[1][cur2]);
            ref2.emplace_back(s2[2][cur2]);
            cur2++;
            cnt = 0;
        } else if(cur1 == s1[0].size()){
            read.emplace_back(s2[0][cur2]);
            m2.emplace_back(s2[1][cur2]);
            ref2.emplace_back(s2[2][cur2]);
            cur2++;
            cnt = 0;
        } else if(cur2 == s2[0].size()) {
            ref1.emplace_back(s1[2][cur1]);
            m1.emplace_back(s1[1][cur1]);
            read.emplace_back(s1[0][cur1]);
            cur1++;
            cnt = 0;
        } else {
            VERIFY(false);
        }
        if(cnt == 20) {
            ref1.insert(ref1.end() - 10, insert.begin(), insert.end());
            m1.insert(m1.end() - 10, insert.begin(), insert.end());
            read.insert(read.end() - 10, insert.begin(), insert.end());
            m2.insert(m2.end() - 10, insert.begin(), insert.end());
            ref2.insert(ref2.end() - 10, insert.begin(), insert.end());
        } else if(cnt > 20) {
            ref1.erase(ref1.end() - 10);
            m1.erase(m1.end() - 10);
            read.erase(read.end() - 10);
            m2.erase(m2.end() - 10);
            ref2.erase(ref2.end() - 10);
        }
    }
    return {std::string(ref1.begin(), ref1.end()), std::string(m1.begin(), m1.end()), std::string(read.begin(), read.end()),
            std::string(m2.begin(), m2.end()), std::string(ref2.begin(), ref2.end())};
}

void AnalyseAndPrint(const Sequence &from_seq, const Sequence &to_seq1, const Sequence &to_seq2, bool extension = false) {
    AlignmentForm cigar1 = defaultAlign(from_seq, to_seq1, extension);
    AlignmentForm cigar2 = defaultAlign(from_seq, to_seq2, extension);
    std::vector<std::pair<size_t, size_t>> to_ignore = Mark(from_seq);
    std::vector<string> res1 = cigar1.toString(from_seq, to_seq1, to_ignore);
    std::vector<string> res2 = cigar2.toString(from_seq, to_seq2, to_ignore);
    std::vector<std::string> res = mixAndShorten(res1, res2);
    std::cout << join("\n", res) << std::endl;
}

bool CheckAndReroute(const Sequence &read_seq, const std::vector<std::pair<size_t, size_t>> &nails, multigraph::GraphPath &path,
                     AlignmentForm &al, const Detour &detour){
    VERIFY(path.getVertex(detour.start) == detour.path.start());
    VERIFY(path.getVertex(detour.end) == detour.path.finish());
    multigraph::GraphPath correction = detour.path;
    multigraph::GraphPath initial = path.subPath(detour.start, detour.end);
    Sequence correctionSeq = correction.Seq();
    if(correctionSeq.size() < nails[detour.start - 1].second + (detour.path.finish().size() - nails[detour.end - 1].second - 1))
           return false;
    correctionSeq = correctionSeq.Subseq(nails[detour.start - 1].second, correctionSeq.size() -
                                                                         (detour.path.finish().size() - nails[detour.end - 1].second - 1));
    Sequence initialSeq = initial.Seq();
    initialSeq = initialSeq.Subseq(nails[detour.start - 1].second, initialSeq.size() -
                                                                   (initial.finish().size() - nails[detour.end - 1].second - 1));
    Sequence subreadSeq = read_seq.Subseq(nails[detour.start - 1].first, nails[detour.end - 1].first + 1);
    std::vector<std::pair<size_t, size_t>> to_ignore = Mark(subreadSeq);
    AlignmentForm initialCigar = defaultAlign(subreadSeq, initialSeq);
    AlignmentForm correctionCigar = defaultAlign(subreadSeq, correctionSeq);
    if(initialCigar.queryLength() == 0 || correctionCigar.queryLength() == 0)
        return false;
    size_t initialScore = Score(initialCigar, subreadSeq, initialSeq, to_ignore);
    size_t correctionScore = Score(correctionCigar, subreadSeq, correctionSeq, to_ignore);
    if(correctionScore < initialScore) {
//        std::cout << "Changed path " << initialScore << " " << correctionScore << " " << std::endl;
//        std::cout << initial.lenStr() << std::endl;
//        std::cout << correction.lenStr() << std::endl;
//        std::cout << join("\n", mixAndShorten(initialCigar.toString(subreadSeq, initialSeq, to_ignore),
//                                              correctionCigar.toString(subreadSeq, correctionSeq, to_ignore))) << std::endl;
        auto lit = al.columnByQpos(nails[detour.start - 1].first);
        auto rit = al.columnByQpos(nails[detour.end - 1].first);
        ++rit;
        AlignmentForm new_al;
        new_al += al.Prefix(lit);
        new_al += correctionCigar;
        new_al += al.Suffix(rit);
        al = std::move(new_al);
        path = path.reroute(detour.start, detour.end, detour.path);
        VERIFY(al.queryLength() == read_seq.size());
        VERIFY(al.targetLength() == path.len());
        std::cout << std::endl;
//        AnalyseAndPrint(alignedSubread, alignedInitial.Seq());
//        AnalyseAndPrint(alignedSubread, alignedDetour.Seq());
//        AnalyseAndPrint(alignedInitial.Seq(), alignedDetour.Seq());
//        AnalyseAndPrint(alignedSubread, alignedInitialSeq, alignedDetourSeq);
        return true;
    }
    return false;
}

int BulgeFinder::INF = std::numeric_limits<int>::max() / 2;

size_t BulgeFinder::getMinDist(const MGVertex &v1, const MGVertex &v2) {
    if(min_dist.find(v1.getId()) == min_dist.end()) {
        typedef std::pair<size_t, multigraph::ConstVertexId> StoredValue;
        std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
        std::unordered_map<multigraph::ConstVertexId, size_t> res;
        queue.emplace(0, v1.getId());
        while(!queue.empty()) {
            StoredValue next = queue.top();
            queue.pop();
            multigraph::ConstVertexId nextVertex = next.second;
            if(res.find(nextVertex) == res.end()) {
                res[nextVertex] = next.first;
                for(const multigraph::MGEdge &edge : *nextVertex) {
                    size_t new_len = next.first + edge.truncSize();
                    if(new_len <= max_size)
                        queue.emplace(new_len, edge.getFinish().getId());
                }
            }
        }
        min_dist[v1.getId()] = std::move(res);
    }
    const std::unordered_map<multigraph::ConstVertexId, size_t> &m2 = min_dist.at(v1.getId());
    if(m2.find(v2.getId()) == m2.end())
        return INF;
    else
        return m2.at(v2.getId());
}

std::vector<Detour> BulgeFinder::findSimpleBulges(const multigraph::GraphPath &path) {
    std::vector<Detour> res;
    for(size_t i = 1; i + 1 < path.size(); i++) {
        multigraph::MGVertex & start = path.getVertex(i);
        multigraph::MGVertex & end = path.getVertex(i + 1);
        for(multigraph::MGEdge &edge : start) {
            if(edge != path.getEdge(i) && edge.getFinish() == end) {
                res.emplace_back(Detour(i, i + 1, {edge}));
            }
        }
    }
    return std::move(res);
}

bool BulgeFinder::recursiveFindBulges(std::vector<multigraph::GraphPath> &bulges, multigraph::GraphPath &bulge, const MGEdge &last_edge,
                                      size_t clen, size_t tlen) {
    if(bulge.finish() == last_edge.getFinish() && bulge.backEdge() != last_edge && tlen <= clen + max_diff <= tlen + 2 * max_diff) {
        bulges.push_back(bulge);
        if(bulges.size() >= 20)
            return false;
    }
    for(multigraph::MGEdge &edge : bulge.finish()) {
        size_t new_len = clen + edge.truncSize();
        size_t min_bulge_len = new_len + getMinDist(edge.getFinish(), last_edge.getFinish());
        if(min_bulge_len <= tlen + max_diff) {
            bulge += edge;
            bool res = recursiveFindBulges(bulges, bulge, last_edge, new_len, tlen);
            bulge.pop_back();
            if(!res)
                return false;
        }
    }
    return true;
}

bool BulgeFinder::recursiveFindBulge(std::vector<Detour> &res, const multigraph::GraphPath &path, size_t start, size_t end,
                                     size_t path_len) {
    if(path_len == INF) {
        path_len = 0;
        for (size_t i = start; i < end; i++) {
            path_len += path.getEdge(i).truncSize();
        }
    }
    if (path_len > max_size) {
        return true;
    }
    std::vector<multigraph::GraphPath> bulges;
    bool found_all = true;
    for(multigraph::MGEdge &edge : path.getVertex(start)) {
        if(edge == path.getEdge(start))
            continue;
        multigraph::GraphPath bulge(edge);
        found_all &=!recursiveFindBulges(bulges, bulge, path.getEdge(end - 1), edge.truncSize(), path_len);
    }
    for(multigraph::GraphPath &bulge : bulges) {
        res.emplace_back(start, end, std::move(bulge));
    }
    return found_all;
}

std::vector<Detour> BulgeFinder::findBulges(const multigraph::GraphPath &path) {
    std::vector<Detour> res;
    bool found_all = true;
    for(size_t i = 1; i + 1 < path.size(); i++) {
        size_t len = 0;
        for(size_t j = i + 1; j + 1 < path.size(); j++) {
            len += path.getEdge(j - 1).truncSize();
            if(len > max_size + path.getVertex(j).size()) {
                break;
            }
            found_all &= recursiveFindBulge(res, path, i, j, len);
        }
    }
    if(res.size() > 0)
        std::cout << path.size() << " " << res.size() << " " << found_all << std::endl;
    return std::move(res);
}

bool changeEnd(multigraph::GraphPath &mpath, const Sequence &from_seq, const std::vector<std::pair<size_t, size_t>> &nails, AlignmentForm &al) {
    MGVertex &last_junction = mpath.getVertex(mpath.size() - 1);
    bool changed = false;
    for(MGEdge &edge: last_junction) {
        if(edge == mpath.backEdge())
            continue;
//        std::cout << "Change end attempt" << std::endl;
        if(edge.fullSize() < mpath.getEdge(mpath.size() - 1).fullSize() - mpath.rightSkip())
            continue;
        Sequence subreadSeq = from_seq.Subseq(nails.back().first, from_seq.size());
        Sequence initialSeq = mpath.backEdge().getSeq().Subseq(nails.back().second,
                                                                      mpath.getEdge(mpath.size() - 1).fullSize() -
                                                                      mpath.rightSkip());
        Sequence correctionSeq = edge.getSeq().Subseq(nails.back().second, edge.fullSize());
        correctionSeq = correctionSeq.Subseq(0, std::min(correctionSeq.size(), subreadSeq.size() + 1000));
        std::vector<std::pair<size_t, size_t>> to_ignore = Mark(subreadSeq);
        AlignmentForm initialCigar = defaultAlignExtension(subreadSeq, initialSeq);
        AlignmentForm correctionCigar = defaultAlignExtension(subreadSeq, correctionSeq);
        if(initialCigar.queryLength() == 0 || correctionCigar.queryLength() == 0)
            return false;
        size_t score1 = Score(initialCigar, subreadSeq, initialSeq, to_ignore);
        size_t score2 = Score(correctionCigar, subreadSeq, correctionSeq, to_ignore);
//        std::cout << join("\n", mixAndShorten(initialCigar.toString(subreadSeq, initialSeq, to_ignore),
//                                              correctionCigar.toString(subreadSeq, correctionSeq, to_ignore))) << std::endl;
        if(initialCigar.queryLength() >= subreadSeq.size() - 5 && correctionCigar.queryLength() >= subreadSeq.size() - 5 && score1 > score2) {
            VERIFY(correctionCigar.queryLength() == subreadSeq.size());
//            std::cout << "Changed path end " << score1 << " " << score2 << " " << Score(initialCigar, subreadSeq, initialSeq, {}) << " " << Score(correctionCigar, subreadSeq, correctionSeq, {}) << std::endl;
//            std::cout << join("\n", mixAndShorten(initialCigar.toString(subreadSeq, initialSeq, to_ignore),
//                                                  correctionCigar.toString(subreadSeq, correctionSeq, to_ignore))) << std::endl;
            mpath.pop_back();
            mpath += edge;
            mpath.cutBack(edge.fullSize() - correctionCigar.targetLength() - nails.back().second);
            auto start = al.columnByQpos(nails.back().first);
            al = al.Prefix(start) + correctionCigar;
            VERIFY(al.queryLength() == from_seq.size());
            VERIFY(al.targetLength() == mpath.len());
            changed = true;
        }
    }
    return changed;
}

multigraph::GraphPath FixPath(const nano::GraphContig &graphContig, BulgeFinder &bulgeFinder, MultiGraph &graph) {
    Sequence from_seq = graphContig.read_str.getSeq().Subseq(graphContig.qStart, graphContig.qEnd);
    multigraph::GraphPath mpath = ContigToPath(graphContig, graph);

    AlignmentForm cigar(graphContig.cigar);
    bool changed = true;
    size_t cnt = 0;
    while(changed && cnt < 10) {
        changed = false;
        std::vector<std::pair<size_t, size_t>> nails = Nails(graph, from_seq, mpath, cigar);
        multigraph::GraphPath correction = mpath;
        for(const Detour &detour : bulgeFinder.findBulges(correction)) {
//            std::cout << "Detour attempt " << detour.path.size() << " " << detour.end - detour.start << std::endl;
            if(CheckAndReroute(from_seq, nails, correction, cigar, detour)) {
                mpath = std::move(correction);
//                cigar = defaultAlign(from_seq, mpath.Seq());
                changed = true;
                cnt++;
                break;
            }
        }
    }
    if(mpath.size() > 1) {
        std::vector<std::pair<size_t, size_t>> nails = Nails(graph, from_seq, mpath, cigar);
        changed = changeEnd(mpath, from_seq, nails, cigar);
//        if(changed)
//            cigar = defaultAlign(from_seq, mpath.Seq());
        mpath = mpath.RC();
        cigar = cigar.RC();
        nails = Nails(graph, !from_seq, mpath, cigar);
        changed = changeEnd(mpath, !from_seq, nails, cigar);
        mpath = mpath.RC();
        cigar = cigar.RC();
    }
    return std::move(mpath);
}