#include <dbg/multi_graph.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>
#include <ksw2/ksw_wrapper.hpp>
#include <utility>
#include <queue>
#include <dbg/paths.hpp>
#include "ReadsAligner.h"

using namespace multigraph;
std::unordered_map<std::string, std::vector<nano::GraphContig>> AlignOnt(logging::Logger &logger, const size_t threads,
                 const std::experimental::filesystem::path &dir, const multigraph::MultiGraph &mg, const io::Library &ont_reads, bool reuse_alignment = false){
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

std::vector<int> CodeToPath(const std::vector<std::string> &code) {
    std::vector<int> res;
    for(const std::string &name : code) {
        int id = std::stoi(name.substr(0, name.size() - 1));
        if (name.back() == '-')
            id = -id;
        res.emplace_back(id);
    }
    return std::move(res);
}

Sequence GraphSeq(multigraph::MultiGraph &graph, const nano::GraphContig &al) {
    size_t skip_left = al.gStart;
    size_t len = al.gEnd - al.gStart;
    SequenceBuilder sb;
    std::vector<int> ids = CodeToPath(al.path);
    for(int id : ids) {
        multigraph::Edge &edge = *graph.getEdgeById(id);
        VERIFY(skip_left < edge.size());
        Sequence to_add = edge.getSeq().Subseq(skip_left);
        if(len <= to_add.size()) {
            sb.append(to_add.Subseq(0, len));
            len = 0;
            skip_left = 0;
            break;
        } else {
            sb.append(to_add);
            skip_left = edge.getFinish().size();
            len -= to_add.size();
        }
    }
    VERIFY(len == 0);
    return sb.BuildSequence();
}
size_t Score(const std::vector<cigar_pair> &cigar, const Sequence &from_seq, const Sequence &to_seq,
             const std::vector<std::pair<size_t, size_t>> &to_ignore) {
    size_t from_pos = 0;
    size_t to_pos = 0;
    size_t diff = 0;
    size_t cur = 0;
    for(const cigar_pair &cp: cigar) {
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

std::pair<size_t, size_t> FromToLen(const std::vector<cigar_pair> &cigar) {
    size_t from_len = 0;
    size_t to_len = 0;
    for(const cigar_pair &cp: cigar) {
        if(cp.type != 'I')
            to_len += cp.length;
        if(cp.type != 'D')
            from_len += cp.length;
    }
    return {from_len, to_len};
}

std::vector<cigar_pair> defaultAlign(const Sequence &from_seq, const Sequence &to_seq, int end_bonus = 0) {
    KSWAligner kswAligner(1, 5, 5, 3);
    return kswAligner.iterativeBandAlign(to_seq.str(), from_seq.str(), 5, 100, 0.01, end_bonus);
}

size_t Score(const Sequence &from_seq, const Sequence &to_seq, const std::vector<std::pair<size_t, size_t>> &to_ignore, int end_bonus = 0) {
    std::vector<cigar_pair> cigar = defaultAlign(from_seq, to_seq, end_bonus);
    return Score(cigar, from_seq, to_seq, to_ignore);
}

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

std::pair<size_t, size_t> Score2(const Sequence &read, const Sequence &p1, const Sequence &p2, int end_bonus = 0) {
    std::vector<std::pair<size_t, size_t>> to_ignore = Mark(read);
    return {Score(read, p1, to_ignore, end_bonus), Score(read, p2, to_ignore, end_bonus)};
}

std::string CigarToString(const std::vector<cigar_pair> &cigar) {
    std::stringstream ss;
    for(const cigar_pair &p : cigar) {
        if(p.length != 1)
            ss << p.length;
        ss << p.type;
    }
    return ss.str();
}


std::vector<cigar_pair> GAStringToCigar(const std::string &s) {
    size_t n = 0;
    std::vector<cigar_pair> res;
    size_t pos = s.find_last_of(':') + 1;
    for(size_t i = pos; i < s.size(); i++){
        char c = s[i];
        if(c >= '0' && c <= '9') {
            n = n * 10 + c - '0';
        } else {
            if (n == 0)
                n = 1;
            if (c == '=' || c == 'X') {
                c = 'M';
            }
            if(!res.empty() && c == res.back().type) {
                res.back() = {c, res.back().length + n};
            } else {
                res.emplace_back(c, n);
            }
            n = 0;
        }
    }
    return std::move(res);
}
MGGraphPath ContigToPath(const nano::GraphContig &al, multigraph::MultiGraph &graph) {
    size_t skip_left = al.gStart;
    size_t skip_right = 0;
    size_t len = al.gEnd - al.gStart;
    SequenceBuilder sb;
    std::vector<int> ids = CodeToPath(al.path);
    std::vector<Edge *> path;
    for(int id : ids) {
        multigraph::Edge &edge = *graph.getEdgeById(id);
        if(path.empty()) {
            if(skip_left >= edge.size() - edge.getFinish().size()) {
                skip_left -= edge.size() - edge.getFinish().size();
            } else {
                path.emplace_back(&edge);
                size_t to_add = edge.size() - skip_left;
                if(len <= to_add) {
                    skip_right = to_add - len;
                    len = 0;
                    break;
                } else {
                    len -= to_add;
                }
            }
        } else {
            path.emplace_back(&edge);
            size_t to_add = edge.size() - edge.getStart().size();
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
    return {path.front()->getStart(), {path.begin(), path.end()}, skip_left, skip_right};
}

std::vector<std::pair<size_t, size_t>> Nails(const multigraph::MultiGraph &graph, const Sequence &from_seq,
                                             const MGGraphPath &mpath, const std::vector<cigar_pair> &cigar) {
    Sequence to_seq = mpath.Seq();
    std::vector<std::pair<size_t, size_t>> vertices;
    size_t skip_left = mpath.leftSkip();
    size_t cur_pos = 0;
    for(multigraph::Edge & edge: mpath.edges()) {
        VERIFY(skip_left < edge.size());
        size_t to_add = edge.size() - skip_left;
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
    for(const cigar_pair &cp: cigar) {
        if(cp.type == 'M') {
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
                while (p < cp.length && from_seq[left - to_pos + from_pos + p] != to_seq[left + p]) {
                    p++;
                }
                if(p == cp.length) {
                    p = 0;
                    while (p < cp.length && from_seq[left - to_pos + from_pos + p] != to_seq[left + p]) {
                        p++;
                    }
                }
                if(p < cp.length && cp.length > best_len[cur1]) {
                    res[cur1] = {left - to_pos + from_pos + p, left + p - vertices[cur1].first};
                    best_len[cur1] = cp.length;
                }
                cur1++;
            }
            from_pos += cp.length;
            to_pos += cp.length;
        }
        if(cp.type == 'D') {
            to_pos += cp.length;
        }
        if(cp.type == 'I') {
            from_pos += cp.length;
        }
        while(cur < vertices.size() && vertices[cur].second <= to_pos) {
            cur++;
        }
    }
    VERIFY(res.size() == vertices.size());
    return std::move(res);
}

struct Detour {
    size_t start;
    size_t end;
    MGGraphPath path;

    Detour(size_t start, size_t anEnd, MGGraphPath path) : start(start), end(anEnd),
                                                                                      path(std::move(path)) {}
};

class BulgeFinder {
private:
    const multigraph::MultiGraph *mg;
    std::unordered_map<ConstVertexId, std::unordered_map<ConstVertexId, size_t>> min_dist;
    size_t max_size;
    size_t max_diff;
    static int INF;

    size_t getMinDist(const Vertex &v1, const Vertex &v2) {
        if(min_dist.find(v1.getId()) == min_dist.end()) {
            typedef std::pair<size_t, ConstVertexId> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            std::unordered_map<ConstVertexId, size_t> res;
            queue.emplace(0, v1.getId());
            while(!queue.empty()) {
                StoredValue next = queue.top();
                queue.pop();
                ConstVertexId nextVertex = next.second;
                if(res.find(nextVertex) == res.end()) {
                    res[nextVertex] = next.first;
                    for(const multigraph::Edge &edge : *nextVertex) {
                        size_t new_len = next.first + edge.truncSize();
                        if(new_len <= max_size)
                            queue.emplace(new_len, edge.getFinish().getId());
                    }
                }
            }
            min_dist[v1.getId()] = std::move(res);
        }
        const std::unordered_map<ConstVertexId, size_t> &m2 = min_dist.at(v1.getId());
        if(m2.find(v2.getId()) == m2.end())
            return INF;
        else
            return m2.at(v2.getId());
    }
public:
    BulgeFinder(multigraph::MultiGraph &mg, size_t max_size, size_t max_diff) : mg(&mg), max_size(max_size), max_diff(max_diff) {
    }

    std::vector<Detour> findSimpleBulges(const MGGraphPath &path) {
        std::vector<Detour> res;
        for(size_t i = 1; i + 1 < path.size(); i++) {
            multigraph::Vertex & start = path.getVertex(i);
            multigraph::Vertex & end = path.getVertex(i+1);
            for(multigraph::Edge &edge : start) {
                if(edge != path.getEdge(i) && edge.getFinish() == end) {
                    res.emplace_back(Detour(i, i + 1, {edge}));
                }
            }
        }
        return std::move(res);
    }

    bool recursiveFindBulges(std::vector<MGGraphPath> &bulges, MGGraphPath &bulge,
                             const multigraph::Edge &last_edge, size_t clen, size_t tlen) {
        if(bulge.finish() == last_edge.getFinish() && bulge.backEdge() != last_edge && tlen <= clen + max_diff <= tlen + 2 * max_diff) {
            bulges.push_back(bulge);
            if(bulges.size() >= 20)
                return false;
        }
        for(multigraph::Edge &edge : bulge.finish()) {
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

    bool recursiveFindBulge(std::vector<Detour> &res, const MGGraphPath &path, size_t start, size_t end, size_t path_len = INF) {
        if(path_len == INF) {
            path_len = 0;
            for (size_t i = start; i < end; i++) {
                path_len += path.getEdge(i).truncSize();
            }
        }
        if (path_len > max_size) {
            return true;
        }
        std::vector<MGGraphPath> bulges;
        bool found_all = true;
        for(multigraph::Edge &edge : path.getVertex(start)) {
            if(edge == path.getEdge(start))
                continue;
            MGGraphPath bulge(edge);
            found_all &=!recursiveFindBulges(bulges, bulge, path.getEdge(end), edge.truncSize(), path_len);
        }
        for(MGGraphPath &bulge : bulges) {
            res.emplace_back(start, end, std::move(bulge));
        }
        return found_all;
    }

    std::vector<Detour> findBulges(const MGGraphPath &path) {
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
};

void AnalyseAndPrint(const Sequence &from_seq, const Sequence &to_seq, int end_bonus = 0) {
    std::vector<cigar_pair> cigar = defaultAlign(from_seq, to_seq, end_bonus);
    std::vector<std::pair<size_t, size_t>> to_ignore = Mark(from_seq);
    std::cout << CigarToString(cigar) << std::endl;
    std::cout << Score(cigar, from_seq, to_seq, to_ignore) << " " << Score(cigar, from_seq, to_seq, {}) << std::endl;
    size_t from_pos = 0;
    size_t to_pos = 0;
    size_t cur = 0;
    std::vector<char> s1, s2, m;
    for(const cigar_pair &cp: cigar) {
        for(size_t t = 0; t < cp.length; t++) {
            if(cur < to_ignore.size() && to_ignore[cur].second == from_pos)
                cur++;
            if(cur < to_ignore.size() && from_pos < to_ignore[cur].second && from_pos >= to_ignore[cur].first) {
                m.emplace_back('*');
            } else if(cp.type == 'M' && from_seq[from_pos] == to_seq[to_pos])
                m.emplace_back('|');
            else
                m.emplace_back('.');
            if(cp.type == 'D') {
                s1.emplace_back('-');
            } else {
                s1.emplace_back("ACGT"[from_seq[from_pos]]);
                from_pos++;
            }
            if(cp.type == 'I') {
                s2.emplace_back('-');
            } else {
                s2.emplace_back("ACGT"[to_seq[to_pos]]);
                to_pos++;
            }
        }
    }
    std::cout << std::string(s1.begin(), s1.end()) << std::endl;
    std::cout << std::string(m.begin(), m.end()) << std::endl;
    std::cout << std::string(s2.begin(), s2.end()) << std::endl;
}

bool CheckAndReroute(const Sequence &read_seq, const std::vector<std::pair<size_t, size_t>> &nails, MGGraphPath &path,
                     const Detour &detour){
    MGGraphPath alignedDetour = detour.path;
    MGGraphPath alignedInitial = path.subPath(detour.start, detour.end);
    alignedDetour.cutFront(nails[detour.start - 1].second).cutBack(detour.path.finish().size() - nails[detour.end - 1].second - 1);
    alignedInitial.cutFront(nails[detour.start - 1].second).cutBack(detour.path.finish().size() - nails[detour.end - 1].second - 1);
    Sequence alignedSubread = read_seq.Subseq(nails[detour.start - 1].first, nails[detour.end - 1].first + 1);
    std::pair<size_t, size_t> scores = Score2(alignedSubread, alignedInitial.Seq(), alignedDetour.Seq());
    if(scores.second < scores.first) {
        std::cout << "Changed path " << scores.first << " " << scores.second << " " << std::endl;
        std::cout << path.str() << std::endl;
        std::cout << alignedInitial.str() << std::endl;
        std::cout << alignedSubread.str() << std::endl;
        path = path.reroute(detour.start, detour.end, detour.path);
        std::cout << std::endl;
        AnalyseAndPrint(alignedSubread, alignedInitial.Seq());
        AnalyseAndPrint(alignedSubread, alignedDetour.Seq());
        AnalyseAndPrint(alignedInitial.Seq(), alignedDetour.Seq());
        return true;
    }
    return false;
}

int BulgeFinder::INF = std::numeric_limits<int>::max() / 2;

bool changeEnd(MGGraphPath &mpath, const Sequence &from_seq, const std::vector<std::pair<size_t, size_t>> &nails) {
    Vertex &last_junction = mpath.getVertex(mpath.size() - 1);
    bool changed = false;
    for(Edge &edge: last_junction) {
        if(edge == mpath.backEdge())
            continue;
        if(edge.size() < mpath.getEdge(mpath.size() - 1).size() - mpath.rightSkip())
            continue;
        Sequence s = from_seq.Subseq(nails.back().first, from_seq.size());
        Sequence s1 = mpath.getEdge(mpath.size() - 1).getSeq().Subseq(nails.back().second,
                                                                   mpath.getEdge(mpath.size() - 1).size() -
                                                                                        mpath.rightSkip());
        Sequence s2 = edge.getSeq().Subseq(nails.back().second, edge.size());
        s2 = s2.Subseq(0, std::min(s2.size(), s.size() + 1000));
        std::vector<std::pair<size_t, size_t>> to_ignore = Mark(s);
        std::vector<cigar_pair> al1 = defaultAlign(s, s1, 1000000);
        std::vector<cigar_pair> al2 = defaultAlign(s, s2, 1000000);
        size_t score1 = Score(al1, s, s1, to_ignore);
        size_t score2 = Score(al2, s, s2, to_ignore);
        std::pair<size_t, size_t> fromto1 = FromToLen(al1);
        std::pair<size_t, size_t> fromto2 = FromToLen(al2);
//            if(fromto1.first < s.size()) {
//                AnalyseAndPrint(s, s1, 1000000);
//            }
        if(fromto1.first >= s.size() - 5 && fromto2.first >= s.size() - 5 &&  score1 > score2) {
            std::cout << "Changed path end " << score1 << " " << score2 << " " << Score(s, s1, {}) << " " << Score(s, s2, {}) << std::endl;
            mpath.pop_back();
            mpath += edge;
            mpath.cutBack(edge.size() - fromto2.second - nails.back().second);
            AnalyseAndPrint(s, s1, 1000000);
            AnalyseAndPrint(s, s2, 1000000);
            changed = true;
        }
    }
    return changed;
}

void FixPath(const nano::GraphContig &graphContig, BulgeFinder &bulgeFinder, multigraph::MultiGraph &graph) {
    std::vector<int> ids = CodeToPath(graphContig.path);
    Sequence from_seq = graphContig.read_str.getSeq().Subseq(graphContig.qStart, graphContig.qEnd);
    MGGraphPath mpath = ContigToPath(graphContig, graph);

    std::vector<cigar_pair> cigar = GAStringToCigar(graphContig.cigar);
    bool changed = true;
    while(changed) {
        changed = false;
        std::vector<std::pair<size_t, size_t>> nails = Nails(graph, from_seq, mpath, cigar);
        MGGraphPath correction = mpath;
        for(const Detour &detour : bulgeFinder.findBulges(correction)) {
            std::cout << detour.path.size() << " " << detour.end - detour.start << std::endl;
            if(CheckAndReroute(from_seq, nails, correction, detour)) {
                mpath = std::move(correction);
                cigar = defaultAlign(from_seq, mpath.Seq());
                changed = true;
                break;
            }
        }
    }
    if(mpath.size() > 1) {
        std::vector<std::pair<size_t, size_t>> nails = Nails(graph, from_seq, mpath, cigar);
        changed = changeEnd(mpath, from_seq, nails);
        if(changed)
            cigar = defaultAlign(from_seq, mpath.Seq());
        mpath = mpath.RC();
        cigar = RcCigar(cigar);
        nails = Nails(graph, !from_seq, mpath, cigar);
        changed = changeEnd(mpath, !from_seq, nails);
        mpath = mpath.RC();
    }
}




int main(int argc, char **argv) {
    AlgorithmParameters parameters({"threads=", "output-dir=", "reuse-alignment"}, {"nano", "graph"}, "");
    CLParser parser(parameters, {"o=output-dir", "t=threads"});
    AlgorithmParameterValues parameterValues = parser.parseCL(argc, argv);
    if (!parameterValues.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parameterValues.checkMissingValues() << "\n" << std::endl;
        std::cout << parameterValues.helpMessage() << std::endl;
        return 1;
    }
    io::Library reads = oneline::initialize<std::experimental::filesystem::path>(parameterValues.getListValue("nano"));
    std::experimental::filesystem::path graph = parameterValues.getListValue("graph")[0];
    std::experimental::filesystem::path dir = parameterValues.getValue("output-dir");
    size_t threads = std::stoull(parameterValues.getValue("threads"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "align");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), logging::trace);
    logger.trace() << "Command line:" << std::endl;
    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    logging::logGit(logger, dir / "version.txt");

    bool reuse = parameterValues.getCheck("reuse-alignment");
    multigraph::MultiGraph mg = MultiGraphHelper::LoadGFA(graph, false);
    mg = MultiGraphHelper::TransformToEdgeGraph(mg);

    logger.info() << "Performing alignment" << std::endl;
    std::unordered_map<std::string, std::vector<nano::GraphContig>> result = AlignOnt(logger, threads, dir, mg, reads, reuse);
    logger.info() << "Alignment finished. Correcting paths." << std::endl;
    std::cout << result.size() << std::endl;

    BulgeFinder bulgeFinder(mg, 5000, 1000);
    for(auto &it : result) {
        for (nano::GraphContig &al: it.second) {
            MGGraphPath path = ContigToPath(al, mg);
//            AnalyseAndPrint(al.read_str.seq.Subseq(al.qStart, al.qEnd), GraphSeq(mg, al));
            FixPath(al, bulgeFinder, mg);
        }
    }
    return 0;
}
