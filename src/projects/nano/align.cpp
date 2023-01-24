#include <lja/multi_graph.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>
#include <ksw2/ksw_wrapper.hpp>
#include <utility>
#include <queue>
#include "ReadsAligner.h"

std::unordered_map<std::string, std::vector<nano::GraphContig>> AlignOnt(logging::Logger &logger, const size_t threads,
                 const std::experimental::filesystem::path &dir, const multigraph::MultiGraph &mg, const io::Library &ont_reads){
    const unsigned int BATCH_SIZE = 10000;
    io::SeqReader reader(ont_reads);
    const std::experimental::filesystem::path &input_gfa = dir / "input.gfa";
    mg.printEdgeGFA(input_gfa);
    logger.info() << "Data loaded" << std::endl;

    std::unordered_map<std::string, Contig> batch;
    nano::ReadsAlignerGA reads_aligner(mg);

    int batch_num = 0;
    StringContig::homopolymer_compressing = true;
    StringContig::min_dimer_to_compress = 32;
    StringContig::max_dimer_size = 32;
    for(StringContig scontig : reader) {
        Contig contig = scontig.makeContig();
        batch[contig.id] = contig;
        if (batch.size() > BATCH_SIZE) {
            batch_num++;
            reads_aligner.Align(batch, input_gfa, dir, threads, batch_num);
            batch.clear();
        }
    }
    reads_aligner.Align(batch, input_gfa, dir, threads, ++batch_num);
    batch.clear();
    std::unordered_map<std::string, std::vector<nano::GraphContig>> result;
    for (int i = 1; i < batch_num + 1; ++ i) {
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

Sequence GraphSeq(const multigraph::MultiGraph &graph, const nano::GraphContig &al) {
    size_t skip_left = al.gStart;
    size_t len = al.gEnd - al.gStart;
    SequenceBuilder sb;
    std::vector<int> ids = CodeToPath(al.path);
    for(int id : ids) {
        const multigraph::Edge &edge = graph.edges.at(id);
        VERIFY(skip_left < edge.size());
        Sequence to_add = edge.getSeq().Subseq(skip_left);
        if(len <= to_add.size()) {
            sb.append(to_add.Subseq(0, len));
            len = 0;
            skip_left = 0;
            break;
        } else {
            sb.append(to_add);
            skip_left = edge.end->size();
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

std::vector<cigar_pair> defaultAlign(const Sequence &from_seq, const Sequence &to_seq) {
    KSWAligner kswAligner(1, 5, 10, 2);
    return kswAligner.iterativeBandAlign(to_seq.str(), from_seq.str(), 5, 100, 0.01);
}

size_t Score(const Sequence &from_seq, const Sequence &to_seq, const std::vector<std::pair<size_t, size_t>> &to_ignore) {
    std::vector<cigar_pair> cigar = defaultAlign(from_seq, to_seq);
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

std::pair<size_t, size_t> Score2(const Sequence &read, const Sequence &p1, const Sequence &p2) {
    std::vector<std::pair<size_t, size_t>> to_ignore = Mark(read);
    return {Score(read, p1, to_ignore), Score(read, p2, to_ignore)};
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

class MPath {
private:
    std::vector<multigraph::Edge *> path;
    size_t skip_left;
    size_t skip_right;
public:
    MPath(const nano::GraphContig &al, multigraph::MultiGraph &graph) {
        skip_left = al.gStart;
        size_t len = al.gEnd - al.gStart;
        SequenceBuilder sb;
        std::vector<int> ids = CodeToPath(al.path);
        for(int id : ids) {
            multigraph::Edge &edge = graph.edges.at(id);
            if(path.empty()) {
                if(skip_left >= edge.size() - edge.end->size()) {
                    skip_left -= edge.size() - edge.end->size();
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
                size_t to_add = edge.size() - edge.start->size();
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
    }

    MPath(std::vector<multigraph::Edge *> path, size_t skipLeft, size_t skipRight) : path(std::move(path)),
                                                                                            skip_left(skipLeft),
                                                                                            skip_right(skipRight) {}

    Sequence getSeq() const {
        SequenceBuilder sb;
        if(path.empty())
            return {};
        if(path.size() == 1) {
            Sequence res =  path[0]->getSeq().Subseq(skip_left, path[0]->size() - skip_right);
            return res;
        }
        sb.append(path[0]->getSeq().Subseq(skip_left));
        for(size_t i = 1; i + 1 < path.size(); i++) {
            sb.append(path[i]->getSeq().Subseq(path[i]->start->size()));
        }
        sb.append(path.back()->getSeq().Subseq(path.back()->start->size(), path.back()->size() - skip_right));
        Sequence seq = sb.BuildSequence();
        return seq;
    }

    std::vector<multigraph::Edge *>::const_iterator begin() const {return path.begin();}
    std::vector<multigraph::Edge *>::const_iterator end() const {return path.end();}
    size_t getSkipLeft() const {return skip_left;}
    size_t getSkipRight() const {return skip_right;}
};


std::vector<std::pair<size_t, size_t>> Nails(const multigraph::MultiGraph &graph, const Sequence &from_seq,
                                 const MPath &mpath, const std::vector<cigar_pair> &cigar) {
    Sequence to_seq = mpath.getSeq();
    std::vector<std::pair<size_t, size_t>> vertices;
    size_t skip_left = mpath.getSkipLeft();
    size_t cur_pos = 0;
    for(multigraph::Edge * edge: mpath) {
        VERIFY(skip_left < edge->size());
        size_t to_add = edge->size() - skip_left;
        cur_pos += to_add;
        vertices.emplace_back(cur_pos - edge->end->size(), cur_pos);
        skip_left = edge->end->size();
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

template<class I>
Sequence BuildSequence(I begin, const I &end, size_t trim_left, size_t trim_right) {
    size_t len = 0;
    for(I tmp = begin; tmp != end; tmp++) {
        multigraph::Edge &edge = **begin;
        len += edge.size() - edge.start->size();
    }
    len += (*begin)->start->size();
    len -= trim_left + trim_right;
    size_t skip_left = trim_left;
    SequenceBuilder sb;
    while(begin != end) {
        multigraph::Edge &edge = **begin;
        if(skip_left < edge.size()) {
            skip_left -= edge.size() - edge.end->size();
            continue;
        }
        Sequence to_add = edge.getSeq().Subseq(skip_left);
        if(len <= to_add.size()) {
            sb.append(to_add.Subseq(0, len));
            len = 0;
            skip_left = 0;
            break;
        } else {
            sb.append(to_add);
            skip_left = edge.end->size();
            len -= to_add.size();
        }
    }
    VERIFY(len == 0);
    return sb.BuildSequence();
}

struct Detour {
    size_t start;
    size_t end;
    std::vector<multigraph::Edge *> path;

    Detour(size_t start, size_t anEnd, std::vector<multigraph::Edge *> path) : start(start), end(anEnd),
                                                                                      path(std::move(path)) {}
};

class BulgeFinder {
private:
    const multigraph::MultiGraph *mg;
    std::unordered_map<int, std::unordered_map<int, int>> min_dist;
    size_t max_size;
    size_t max_diff;
    static int INF;

    size_t getMinDist(int v1, int v2) {
        if(min_dist.find(v1) == min_dist.end()) {
            typedef std::pair<int, int> StoredValue;
            std::priority_queue<StoredValue, std::vector<StoredValue>, std::greater<>> queue;
            std::unordered_map<int, size_t> res;
            queue.emplace(0, v1);
            while(!queue.empty()) {
                StoredValue next = queue.top();
                queue.pop();
                if(res.find(next.second) == res.end()) {
                    res[next.second] = next.first;
                    for(multigraph::Edge *edge : mg->vertices.at(next.second).outgoing) {
                        int new_len = next.first + int(edge->size()) - edge->end->size();
                        queue.emplace(new_len, edge->end->id);
                    }
                }
            }
        }
        const std::unordered_map<int, int> &m2 = min_dist.at(v1);
        if(m2.find(v2) == m2.end())
            return INF;
        else
            return m2.at(v2);
    }
public:
    BulgeFinder(multigraph::MultiGraph &mg, size_t max_size, size_t max_diff) : mg(&mg), max_size(max_size), max_diff(max_diff) {
    }

    std::vector<Detour> findSimpleBulges(const std::vector<multigraph::Edge *> &path) {
        std::vector<Detour> res;
        for(size_t i = 1; i + 1 < path.size(); i++) {
            multigraph::Vertex & start = *path[i]->start;
            multigraph::Vertex & end = *path[i]->end;
            for(multigraph::Edge * edge : start.outgoing) {
                if(edge != path[i] && edge->end == &end) {
                    res.emplace_back(Detour(i, i + 1, {edge}));
                }
            }
        }
        return std::move(res);
    }

    bool recoursiveFindBulges(std::vector<std::vector<multigraph::Edge *>> &bulges, std::vector<multigraph::Edge *> &bulge,
                          const multigraph::Edge &last_edge, int clen, int tlen) {
        if(bulge.back()->end == last_edge.end && bulge.back() != &last_edge && tlen - max_diff <= clen <= tlen + max_diff) {
            bulges.emplace_back(bulge);
            if(bulges.size() >= 20)
                return false;
        }
        for(multigraph::Edge *edge : bulge.back()->end->outgoing) {
            int new_len = clen + edge->size() - edge->end->size();
            int min_bulge_len = new_len + getMinDist(edge->end->id, last_edge.end->id);
            if(min_bulge_len <= tlen + max_diff) {
                bulge.push_back(edge);
                bool res = recoursiveFindBulges(bulges, bulge, last_edge, new_len, tlen);
                bulge.pop_back();
                if(!res)
                    return false;
            }
        }
        return true;
    }

    bool recoursiveFindBulge(std::vector<Detour> &res, const std::vector<multigraph::Edge *> &path, size_t start, size_t end, int path_len = -INF) {
        if(path_len == -INF) {
            int path_len = -int(path[start]->start->size());
            for (size_t i = start; i < end; i++) {
                path_len += int(path[i]->size()) - int(path[i]->end->size());
                if (path_len > max_size) {
                    return true;
                }
            }
        }
        std::vector<std::vector<multigraph::Edge *>> bulges;
        bool found_all = true;
        for(multigraph::Edge *edge : path[start]->start->outgoing) {
            if(edge == path[start])
                continue;
            std::vector<multigraph::Edge*> bulge = {edge};
            found_all &=!recoursiveFindBulges(bulges, bulge, *path[end - 1], int(edge->size() - edge->start->size() - edge->end->size()), path_len);
        }
        for(std::vector<multigraph::Edge *> &bulge : bulges) {
            res.emplace_back(start, end, std::move(bulge));
        }
        return found_all;
    }

    std::vector<Detour> findBulges(const std::vector<multigraph::Edge *> &path) {
        std::vector<Detour> res;
        bool found_all;
        for(size_t i = 1; i + 1 < path.size(); i++) {
            int len = -int(path[i]->start->size());
            for(size_t j = i + 1; j + 1 < path.size(); j++) {
                len += int(path[j]->size()) - int(path[j]->end->size());
                if(len > max_size) {
                    break;
                }
                found_all &= recoursiveFindBulge(res, path, i, j, len);
            }
        }
        std::cout << path.size() << " " << res.size() << " " << found_all << std::endl;
        return std::move(res);
    }
};

void AnalyseAndPrint(const Sequence &from_seq, const Sequence &to_seq) {
    std::vector<cigar_pair> cigar = defaultAlign(from_seq, to_seq);
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
                s1.emplace_back("ACTG"[from_seq[from_pos]]);
                from_pos++;
            }
            if(cp.type == 'I') {
                s2.emplace_back('-');
            } else {
                s2.emplace_back("ACTG"[to_seq[to_pos]]);
                to_pos++;
            }
        }
    }
    std::cout << std::string(s1.begin(), s1.end()) << std::endl;
    std::cout << std::string(m.begin(), m.end()) << std::endl;
    std::cout << std::string(s2.begin(), s2.end()) << std::endl;
}

bool CheckAndReplace(const Sequence &read_seq, const std::vector<std::pair<size_t, size_t>> &nails, std::vector<multigraph::Edge *> &path,
                     const Detour &detour){
    Sequence s = read_seq.Subseq(nails[detour.start - 1].first, nails[detour.end - 1].first + 1);
//    Sequence s1 = BuildSequence(std::next(path.begin(), detour.start), std::next(path.begin(), detour.end),
//                                nails[detour.start].second, path[detour.end]->end->size() - nails[detour.end].second - 1);
//    Sequence s2 = BuildSequence(detour.path.begin(), detour.path.end(), nails[detour.start].second, path[detour.end]->end->size() - nails[detour.end].second - 1);
    Sequence s1 = MPath({path.begin() + detour.start, path.begin() + detour.end}, nails[detour.start - 1].second, path[detour.end]->start->size() - nails[detour.end - 1].second - 1).getSeq();
    Sequence s2 = MPath(detour.path, nails[detour.start - 1].second, path[detour.end]->start->size() - nails[detour.end - 1].second - 1).getSeq();
    std::pair<size_t, size_t> scores = Score2(s, s1, s2);
    if(scores.second < scores.first) {
        std::cout << "Changed path " << scores.first << " " << scores.second << " " << Score(s, s1, {}) << " " << Score(s, s2, {}) << std::endl;
        for(multigraph::Edge *edge : path) {
            std::cout << edge->getId() << " ";
        }
        std::cout << std::endl;
        std::vector<multigraph::Edge *> res(path.begin(), path.begin() + detour.start);
        res.insert(res.end(), detour.path.begin(), detour.path.end());
        res.insert(res.end(), path.begin() + detour.end, path.end());
        path = std::move(res);
        for(multigraph::Edge *edge : path) {
            std::cout << edge->getId() << " ";
        }
        std::cout << std::endl;
        AnalyseAndPrint(s, s1);
        AnalyseAndPrint(s, s2);
        return true;
    }
    return false;
}

int BulgeFinder::INF = std::numeric_limits<int>::max() / 2;

void FixPath(const nano::GraphContig &graphContig, BulgeFinder &bulgeFinder, multigraph::MultiGraph &graph) {
    std::vector<int> ids = CodeToPath(graphContig.path);
    Sequence from_seq = graphContig.read_str.seq.Subseq(graphContig.qStart, graphContig.qEnd);
    MPath mpath(graphContig, graph);
    std::vector<cigar_pair> cigar = GAStringToCigar(graphContig.cigar);
    bool changed = true;
    while(changed) {
        changed = false;
        std::vector<std::pair<size_t, size_t>> nails = Nails(graph, from_seq, mpath, cigar);
        std::vector<multigraph::Edge *> path(mpath.begin(), mpath.end());
        for(const Detour &detour : bulgeFinder.findBulges(path)) {
            if(CheckAndReplace(from_seq, nails, path, detour)) {
                mpath = {std::move(path), mpath.getSkipLeft(), mpath.getSkipRight()};
                cigar = defaultAlign(from_seq, mpath.getSeq());
                changed = true;
                break;
            }
        }
    }
}


int main(int argc, char **argv) {
    AlgorithmParameters parameters({"threads=", "output-dir="}, {"nano", "graph"}, "");
    CLParser parser(parameters, {"o=output-dir", "t=threads"});
    AlgorithmParameterValues parameterValues = parser.parseCL(argc, argv);
    parameterValues.checkMissingValues();
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

    multigraph::MultiGraph mmg;
    mmg.LoadGFA(graph, false);
    multigraph::MultiGraph mg = mmg.DBG();

    std::unordered_map<std::string, std::vector<nano::GraphContig>> result = AlignOnt(logger, threads, dir, mg, reads);
    BulgeFinder bulgeFinder(mg, 5000, 1000);
    for(auto &it : result) {
        std::cout << it.first << " " << it.second.size() << std::endl;
        for(nano::GraphContig &al : it.second) {
            for(std::string & s: al.path) {
                std::cout << s << " ";
            }
            std::cout << std::endl;
            MPath mPath(al, mg);
//            AnalyseAndPrint(al.read_str.seq.Subseq(al.qStart, al.qEnd), GraphSeq(mg, al));
            FixPath(al, bulgeFinder, mg);
        }
    }
    return 0;
}