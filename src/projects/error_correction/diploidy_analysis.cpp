#include "diploidy_analysis.hpp"

//BulgePath::BulgePath(std::vector<std::pair<dbg::Edge *, dbg::Edge *>> &&path_) : path(path_), start_(nullptr) {
//    VERIFY(path.size() > 0);
//    start_ = &path.front().first->getStart();
//}
//
//dbg::Vertex &BulgePath::getFinish() const {
//    if(path.empty())
//        return *start_;
//    return path.back().first->getFinish();
//}
//
//dbg::Vertex &BulgePath::getStart() const {
//    return *start_;
//}
//
//dbg::Vertex &BulgePath::getVertex(size_t ind) const {
//    VERIFY(ind <= size());
//    if(ind == size())
//        return getFinish();
//    return path[ind].first->getStart();
//}
//
//void BulgePath::extend(double threshold) {
//    dbg::Vertex &last = getFinish();
//    size_t deg = last.outDeg();
//    if(last.front().getCoverage() > threshold && last.back().getCoverage() > threshold)
//        path.emplace_back(&last.front(), &last.back());
//    else {
//        for(dbg::Edge &edge: last) {
//            if(edge.getCoverage() > threshold) {
//                path.emplace_back(&edge, &edge);
//                return;
//            }
//        }
//        VERIFY(last.outDeg() == 2 && last.front().getFinish() == last.back().getFinish());
//        path.emplace_back(&last.front(), &last.back());
//    }
//}
//
//BulgePath BulgePath::RC() {
//    if(path.empty()) {
//        return BulgePath(start_->rc());
//    }
//    std::vector<std::pair<dbg::Edge *, dbg::Edge *>> rc;
//    for(size_t i = 0; i < path.size(); i++) {
//        rc.emplace_back(&path[path.size() - 1 - i].first->rc(), &path[path.size() - 1 - i].second->rc());
//    }
//    return BulgePath(std::move(rc));
//}
//
//BulgePath BulgePath::operator+(const BulgePath &other) const {
//    VERIFY(getFinish() == other.getStart())
//    std::vector<std::pair<dbg::Edge *, dbg::Edge *>> sum(path);
//    sum.insert(sum.end(), other.path.begin(), other.path.end());
//    return BulgePath(std::move(sum));
//}
//
//dbg::Vertex &BulgePath::vertexAt(size_t ind) {
//    if(ind == 0)
//        return *start_;
//    return path[ind - 1].first->getFinish();
//}
//
//size_t BulgePath::length() const {
//    size_t res = 0;
//    for(auto & p : path) {
//        res += std::max(p.first->truncSize(), p.second->truncSize());
//    }
//    return res;
//}
//
//size_t BulgePath::bulgeLength() const {
//    size_t res = 0;
//    for(auto & p : path) {
//        if(p.first != p.second)
//            res += std::max(p.first->truncSize(), p.second->truncSize());
//    }
//    return res;
//}
//
//size_t BulgePath::conservativeLength() const {
//    size_t res = 0;
//    for(auto & p : path) {
//        if(p.first == p.second)
//            res += std::max(p.first->truncSize(), p.second->truncSize());
//    }
//    return res;
//}
//
//std::string BulgePath::str() const {
//    std::stringstream ss;
//    ss << getStart().getShortId();
//    for(const auto &p : path) {
//        if(p.first == p.second) {
//            ss << "-" << p.first->truncSize() << p.first->firstNucl() << "-" << p.first->getFinish().getShortId();
//        } else {
//            ss << "-(" << p.first->truncSize() << p.first->firstNucl() << "," <<
//               p.second->truncSize() << p.second->firstNucl() << ")-" << p.first->getFinish().getShortId();
//        }
//    }
//    return ss.str();
//}
//
//bool BulgePath::isBad(size_t bad_bulge_inner_size) const {
//    if(path.size() < 2)
//        return false;
//    for(const auto &p : path) {
//        if(p.first != p.second)
//            if(p.first->innerSize() > bad_bulge_inner_size || p.second->truncSize() > bad_bulge_inner_size) {
//                return false;
//            }
//    }
//    return true;
//}
//
//dbg::GraphPath BulgePath::randomPath() const {
//    dbg::GraphPath res(getStart());
//    for(const std::pair<dbg::Edge *, dbg::Edge *> &pair: path) {
//        res += *pair.first;
//    }
//    return std::move(res);
//}

BulgePath<dbg::DBGTraits> BulgePathFinder::forwardPath(dbg::Vertex &start) {
    BulgePath<dbg::DBGTraits> res(start);
    dbg::Vertex * cur = &start;
    while(isBulgePathInner(*cur)) {
        res.extend(threshold);
        cur = &res.getFinish();
        if(cur == &start)
            return std::move(res);
    }
    return std::move(res);
}

BulgePathFinder::BulgePathFinder(ag::AssemblyGraph<dbg::DBGTraits> &dbg, double threshold) : dbg(dbg), threshold(threshold) {
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
    for(const BulgePath<dbg::DBGTraits> &bp : paths) {
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

std::pair<std::vector<dbg::EdgeId>, std::vector<dbg::EdgeId>>
BulgePathCorrector::resolveBulgePath(const dbg::DBGAlignedReadStorage &reads, const BulgePath<dbg::DBGTraits> &path) const {
    dbg::GraphPath p1; dbg::GraphPath p2;
    dbg::GraphPath repeat;
    for (size_t i = 0; i < path.size(); i++) {
        if (path.isBulge(i)) {
            if (p1.empty() && p2.empty()) {
                p1 = repeat + *path[i].first; p2 = repeat + *path[i].second;
            } else {
                dbg::Edge & e1 = *path[i].first;
                dbg::Edge & e2 = *path[i].second;
                const ag::SuffixRecord<dbg::DBGTraits> &rec1 = reads.getSuffixes().getSuffixRecord(p1.backEdge());
                const ag::SuffixRecord<dbg::DBGTraits> &rec2 = reads.getSuffixes().getSuffixRecord(p2.backEdge());
                size_t straight_score = rec1.countStartsWith(repeat + e1) + rec2.countStartsWith(repeat + e2);
                size_t switch_score = rec1.countStartsWith(repeat + e2) + rec2.countStartsWith(repeat + e1);
                p1 += repeat; p2 += repeat;
                if (straight_score >= switch_score) { p1 += e1; p2 += e2; }
                else { p1 += e2; p2 += e1; }
            }
            repeat = dbg::GraphPath();
        } else {
            repeat += *path[i].first;
        }
    }
    p1 += repeat; p2 += repeat;
    return {p1.asEdgeIds(), p2.asEdgeIds()};
}

//    TODO: get rid of this or at least do alignment of ends.
std::string BulgePathCorrector::correctRead(const std::string &name, dbg::GraphPath &read_path) {
    std::vector<Case> cases;
    std::vector<std::string> messages;
    std::vector<Segment<dbg::Edge>> path = oneline::initialize<Segment<dbg::Edge>>(read_path.begin(), read_path.end());
    for(size_t i = 0; i < path.size(); i++) {
        if(!cases.empty() && cases.back().read_to == i && cases.back().path_to != paths[cases.back().path_ind].size() &&
           (paths[cases.back().path_ind][cases.back().path_to].first == path[i].contig().getId() ||
           paths[cases.back().path_ind][cases.back().path_to].second == path[i].contig().getId())) {
            if(resolved[cases.back().path_ind].first[cases.back().path_to] != path[i].contig().getId())
                cases.back().score1 += 1;
            if(resolved[cases.back().path_ind].second[cases.back().path_to] != path[i].contig().getId())
                    cases.back().score2 += 1;
                cases.back().read_to += 1;
                cases.back().path_to += 1;
            } else {
                auto it = pathPoses.find(path[i].contig().getId());
                if (it != pathPoses.end()) {
                    size_t score1 = (resolved[it->second.path_ind].first[it->second.pos] != path[i].contig().getId());
                    size_t score2 = (resolved[it->second.path_ind].second[it->second.pos] != path[i].contig().getId());
                    cases.emplace_back(it->second.path_ind, it->second.pos, it->second.pos + 1, i, i + 1, score1,
                                       score2);
                }
            }
    }
    if(cases.empty())
        return "";
    dbg::GraphPath res;
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
                    res += *resolved[bp.path_ind].first[i - bp.read_from + bp.path_from];
            } else {
                for(size_t i = bp.read_from; i < bp.read_to; i++)
                    res += *resolved[bp.path_ind].second[i - bp.read_from + bp.path_from];
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
    for(BulgePath<dbg::DBGTraits> &path: BulgePathFinder(dbg, threshold).paths) {
        if (path.size() > 1 && path.length() > unique_length) {
            paths.emplace_back(path);
        }
    }
    for(size_t path_ind = 0; path_ind < paths.size(); path_ind++) {
        BulgePath<dbg::DBGTraits> &path = paths[path_ind];
        for(size_t i = 0; i < path.size(); i++) {
            std::pair<dbg::EdgeId, dbg::EdgeId> pair = path[i];
            pathPoses[pair.first] = {path_ind, i};
            pathPoses[pair.second] = {path_ind, i};
        }
        resolved.emplace_back(resolveBulgePath(reads, path));
    }
}

