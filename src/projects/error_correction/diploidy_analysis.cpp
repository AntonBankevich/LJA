#include "diploidy_analysis.hpp"

using namespace ag;
using namespace dbg;

BulgePath BulgePathFinder::forwardPath(dbg::Vertex &start) {
    BulgePath res(start);
    dbg::Vertex * cur = &start;
    while(isBulgePathInner(*cur)) {
        res.extend(threshold);
        cur = &res.getFinish();
        if(cur == &start)
            return std::move(res);
    }
    return std::move(res);
}

BulgePathFinder::BulgePathFinder(ag::AssemblyGraph &dbg, double threshold) : dbg(dbg), threshold(threshold) {
    std::unordered_set<dbg::Vertex *> visited;
    for(auto &vertex : dbg.verticesUnique()) {
        if(visited.find(&vertex) != visited.end())
            continue;
        if(isBulgePathInner(vertex)) {
            BulgePath new_path = forwardPath(vertex);
            VERIFY(new_path.size() > 0);
            if(new_path.getStart() != new_path.getFinish()) {
                BulgePath p2 = forwardPath(vertex.rc());
                BulgePath p3 = p2.RC();
                new_path = p3 + new_path;
            }
            for(size_t i = 1; i + 1 <= new_path.size(); i++) {
                visited.emplace(&new_path.vertexAt(i));
                visited.emplace(&new_path.vertexAt(i).rc());
            }
            if(new_path.getStart() == new_path.getFinish()) {
                visited.emplace(&new_path.getStart());
                visited.emplace(&new_path.getStart().rc());
            }
            paths.emplace_back(new_path.RC());
            paths.emplace_back(std::move(new_path));
        }
    }
    for(dbg::Edge &edge : dbg.edges()) {
        if(visited.find(&edge.getFinish()) == visited.end() && visited.find(&edge.getStart()) == visited.end())
            paths.emplace_back(edge);
    }
}

SetUniquenessStorage BulgePathFinder::uniqueEdges(size_t min_len) const {
    std::vector<dbg::EdgeId> res;
    for(const BulgePath &bp : paths) {
        if(bp.size() == 1) {
            dbg::Edge &edge = *bp[0].first;
            if(edge.truncSize() > min_len || (
                    (edge.getStart().inDeg() == 0 || edge.getFinish().outDeg() == 0) &&
                    edge.truncSize() > min_len / 3 && edge.getCoverage() > 4)) {
                res.emplace_back(edge.getId());
            }
        } else {
            if(!bp.isBad(0) && (bp.conservativeLength() < bp.length() / 2 || bp.size() >= 4) && bp.length() > min_len && bp.conservativeLength() < bp.length() * 95 / 100) {
                for (auto &p : bp) {
                    if (p.first != p.second) {
                        res.emplace_back(p.first);
                        res.emplace_back(p.second);
                    }
                }
            }
        }
    }
    return {res.begin(), res.end()};
}

std::pair<ag::RAGraphPath, ag::RAGraphPath>
BulgePathCorrector::resolveBulgePath(const dbg::DBGAlignedReadStorage &reads, const BulgePath &path) const {
    ag::GraphPath p1; ag::GraphPath p2;
    ag::GraphPath repeat;
    for (size_t i = 0; i < path.size(); i++) {
        if (path.isBulge(i)) {
            if (p1.empty() && p2.empty()) {
                p1 = repeat + *path[i].first; p2 = repeat + *path[i].second;
            } else {
                dbg::Edge & e1 = *path[i].first;
                dbg::Edge & e2 = *path[i].second;
                const ag::SuffixRecord &rec1 = reads.getSuffixes().getSuffixRecord(p1.backEdge());
                const ag::SuffixRecord &rec2 = reads.getSuffixes().getSuffixRecord(p2.backEdge());
                size_t straight_score = rec1.countStartsWith(repeat + e1) + rec2.countStartsWith(repeat + e2);
                size_t switch_score = rec1.countStartsWith(repeat + e2) + rec2.countStartsWith(repeat + e1);
                p1 += repeat; p2 += repeat;
                if (straight_score >= switch_score) { p1 += e1; p2 += e2; }
                else { p1 += e2; p2 += e1; }
            }
            repeat = ag::GraphPath();
        } else {
            repeat += *path[i].first;
        }
    }
    p1 += repeat; p2 += repeat;
    return {p1.asRAPath(), p2.asRAPath()};
}

//    TODO: get rid of this or at least do alignment of ends.
std::string BulgePathCorrector::correctRead(const std::string &name, ag::GraphPath &read_path) {
    std::vector<Case> cases;
    std::vector<std::string> messages;
    std::vector<Segment<dbg::Edge>> path = oneline::initialize<Segment<dbg::Edge>>(read_path.begin(), read_path.end());
    for(size_t i = 0; i < path.size(); i++) {
        if(!cases.empty() && cases.back().read_to == i && cases.back().path_to != paths[cases.back().path_ind].size() &&
           (paths[cases.back().path_ind][cases.back().path_to].first == path[i].contig().getId() ||
           paths[cases.back().path_ind][cases.back().path_to].second == path[i].contig().getId())) {
            if(resolved[cases.back().path_ind].first[cases.back().path_to] != path[i].contig())
                cases.back().score1 += 1;
            if(resolved[cases.back().path_ind].second[cases.back().path_to] != path[i].contig())
                    cases.back().score2 += 1;
                cases.back().read_to += 1;
                cases.back().path_to += 1;
            } else {
                auto it = pathPoses.find(path[i].contig().getId());
                if (it != pathPoses.end()) {
                    size_t score1 = (resolved[it->second.path_ind].first[it->second.pos] != path[i].contig());
                    size_t score2 = (resolved[it->second.path_ind].second[it->second.pos] != path[i].contig());
                    cases.emplace_back(it->second.path_ind, it->second.pos, it->second.pos + 1, i, i + 1, score1,
                                       score2);
                }
            }
    }
    if(cases.empty())
        return "";
    ag::GraphPath res;
    for(Case & bp : cases) {
        for(size_t i = res.calculateSize(); i < bp.read_from; i++)
            res += path[i];
        if(bp.score1 == 0 || bp.score2 == 0) {
            for(size_t i = res.calculateSize(); i < bp.read_to; i++)
                res += path[i];
        } else {
            messages.emplace_back("bpc" + itos(std::min(bp.score1, bp.score2)));
            if(bp.score1 <= bp.score2) {
                for(size_t i = bp.read_from; i < bp.read_to; i++)
                    res += resolved[bp.path_ind].first[i - bp.read_from + bp.path_from];
            } else {
                for(size_t i = bp.read_from; i < bp.read_to; i++)
                    res += resolved[bp.path_ind].second[i - bp.read_from + bp.path_from];
            }
        }
    }
    for(size_t i = res.calculateSize(); i < path.size(); i++)
        res += path[i];
    if(res.front() != path.front() && res.front().size() > path.front().size()) {
        res.setCutLeft(res.leftCut() + res.front().size() - path.front().size());
        }
        if (res.back() != path.back() && res.back().size() > path.back().size()) {
            res.setCutRight(res.rightCut() + res.back().size() - path.back().size());
        }
        read_path = std::move(res);
    return join("_", messages);
}

void BulgePathCorrector::initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                                    dbg::DBGAlignedReadStorage &reads) {
    paths = BulgePathFinder(dbg, threshold).paths;
    for(BulgePath &path: BulgePathFinder(dbg, threshold).paths) {
        if (path.size() > 1 && path.length() > unique_length) {
            paths.emplace_back(path);
        }
    }
    for(size_t path_ind = 0; path_ind < paths.size(); path_ind++) {
        BulgePath &path = paths[path_ind];
        for(size_t i = 0; i < path.size(); i++) {
            std::pair<dbg::EdgeId, dbg::EdgeId> pair = path[i];
            pathPoses[pair.first] = {path_ind, i};
            pathPoses[pair.second] = {path_ind, i};
        }
        resolved.emplace_back(resolveBulgePath(reads, path));
    }
}

