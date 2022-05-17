#include "diploidy_analysis.hpp"

BulgePath::BulgePath(std::vector<std::pair<dbg::Edge *, dbg::Edge *>> &&path_) : path(path_), start_(nullptr) {
    VERIFY(path.size() > 0);
    start_ = path.front().first->start();
}

dbg::Vertex &BulgePath::finish() const {
    if(path.empty())
        return *start_;
    return *path.back().first->end();
}

dbg::Vertex &BulgePath::start() const {
    return *start_;
}

dbg::Vertex &BulgePath::getVertex(size_t ind) const {
    VERIFY(ind <= size());
    if(ind == size())
        return finish();
    return *path[ind].first->start();
}

void BulgePath::extend(double threshold) {
    dbg::Vertex &last = finish();
    size_t deg = last.outDeg();
    if(last[0].getCoverage() > threshold && last[deg - 1].getCoverage() > threshold)
        path.emplace_back(&last[0], &last[deg - 1]);
    else {
        for(dbg::Edge &edge: last) {
            if(edge.getCoverage() > threshold) {
                path.emplace_back(&edge, &edge);
                return;
            }
        }
        VERIFY(last.outDeg() == 2 && last[0].end() == last[1].end());
        path.emplace_back(&last[0], &last[1]);
    }
}

BulgePath BulgePath::RC() {
    if(path.empty()) {
        return BulgePath(start_->rc());
    }
    std::vector<std::pair<dbg::Edge *, dbg::Edge *>> rc;
    for(size_t i = 0; i < path.size(); i++) {
        rc.emplace_back(&path[path.size() - 1 - i].first->rc(), &path[path.size() - 1 - i].second->rc());
    }
    return BulgePath(std::move(rc));
}

BulgePath BulgePath::operator+(const BulgePath &other) const {
    VERIFY(finish() == other.start())
    std::vector<std::pair<dbg::Edge *, dbg::Edge *>> sum(path);
    sum.insert(sum.end(), other.path.begin(), other.path.end());
    return BulgePath(std::move(sum));
}

dbg::Vertex &BulgePath::vertexAt(size_t ind) {
    if(ind == 0)
        return *start_;
    return *path[ind - 1].first->end();
}

size_t BulgePath::length() const {
    size_t res = 0;
    for(auto & p : path) {
        res += std::max(p.first->size(), p.second->size());
    }
    return res;
}

size_t BulgePath::bulgeLength() const {
    size_t res = 0;
    for(auto & p : path) {
        if(p.first != p.second)
            res += std::max(p.first->size(), p.second->size());
    }
    return res;
}

size_t BulgePath::conservativeLength() const {
    size_t res = 0;
    for(auto & p : path) {
        if(p.first == p.second)
            res += std::max(p.first->size(), p.second->size());
    }
    return res;
}

std::string BulgePath::str() const {
    std::stringstream ss;
    ss << start().getShortId();
    for(const auto &p : path) {
        if(p.first == p.second) {
            ss << "-" << p.first->size() << "ACGT"[p.first->seq[0]] << "-" << p.first->end()->getShortId();
        } else {
            ss << "-(" << p.first->size() << "ACGT"[p.first->seq[0]] << "," <<
               p.second->size() << "ACGT"[p.second->seq[0]] << ")-" << p.first->end()->getShortId();
        }
    }
    return ss.str();
}

bool BulgePath::isBad(size_t bad_bulge_size) const {
    if(path.size() < 2)
        return false;
    for(const auto &p : path) {
        if(p.first != p.second)
            if(p.first->size() > bad_bulge_size || p.second->size() > bad_bulge_size) {
                return false;
            }
    }
    return true;
}

dbg::Path BulgePath::randomPath() const {
    dbg::Path res(start());
    for(const std::pair<dbg::Edge *, dbg::Edge *> &pair: path) {
        res += *pair.first;
    }
    return std::move(res);
}

bool BulgePathFinder::checkVertexForward(const dbg::Vertex &v) {
    size_t cnt = 0;
    for(dbg::Edge &edge : v) {
        if(edge.getCoverage() > threshold)
            cnt++;
    }
    if(cnt > 1 && cnt != v.outDeg())
        return false;
    return cnt == 1 ||
           (v.outDeg() == 2 && v[0].end() == v[1].end() && v[0].size() < v[1].size() * 1.3 && v[1].size() < v[0].size() * 1.3);
}

bool BulgePathFinder::isInner(const dbg::Vertex &v) {
    return checkVertexForward(v) && checkVertexForward(v.rc());
}

BulgePath BulgePathFinder::forwardPath(dbg::Vertex &start) {
    BulgePath res(start);
    dbg::Vertex * cur = &start;
    while(isInner(*cur)) {
        res.extend(threshold);
        cur = &res.finish();
        if(cur == &start)
            return std::move(res);
    }
    return std::move(res);
}

BulgePathFinder::BulgePathFinder(dbg::SparseDBG &dbg, double threshold) : dbg(dbg), threshold(threshold) {
    std::unordered_set<dbg::Vertex *> visited;
    for(auto &it : dbg) {
        if(visited.find(&it.second) != visited.end())
            continue;
        if(isInner(it.second)) {
            BulgePath new_path = forwardPath(it.second);
            VERIFY(new_path.size() > 0);
            if(new_path.start() != new_path.finish()) {
                BulgePath p2 = forwardPath(it.second.rc());
                BulgePath p3 = p2.RC();
                new_path = p3 + new_path;
            }
            for(size_t i = 1; i + 1 <= new_path.size(); i++) {
                visited.emplace(&new_path.vertexAt(i));
                visited.emplace(&new_path.vertexAt(i).rc());
            }
            if(new_path.start() == new_path.finish()) {
                visited.emplace(&new_path.start());
                visited.emplace(&new_path.start().rc());
            }
            paths.emplace_back(new_path.RC());
            paths.emplace_back(std::move(new_path));
        }
    }
    for(dbg::Edge &edge : dbg.edges()) {
        if(visited.find(edge.end()) == visited.end() && visited.find((edge.start())) == visited.end())
            paths.emplace_back(edge);
    }
}

SetUniquenessStorage BulgePathFinder::uniqueEdges(size_t min_len) const {
    std::vector<dbg::Edge *> res;
    for(const BulgePath &bp : paths) {
        if(bp.size() == 1) {
            dbg::Edge &edge = *bp[0].first;
            if(edge.size() > min_len || (
                    (edge.start()->inDeg() == 0 || edge.end()->outDeg() == 0) &&
                    edge.size() > min_len / 3 && edge.getCoverage() > 4)) {
                res.emplace_back(&edge);
            }
        } else {
            if(!bp.isBad(dbg.hasher().getK()) && (bp.conservativeLength() < bp.length() / 2 || bp.size() >= 4) && bp.length() > min_len && bp.conservativeLength() < bp.length() * 95 / 100) {
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

std::pair<dbg::Path, dbg::Path>
BulgePathCorrector::resolveBulgePath(const RecordStorage &reads, const BulgePath &path) const {
    dbg::Path p1(path.start());
    dbg::Path p2(path.start());
    Sequence seq;
    size_t last = 0;
    Sequence s1;
    Sequence s2;
    for(size_t i = 0; i < path.size(); i++) {
        if(path.isBulge(i)) {
            if(s1 == s2) {
                p1 += *path[i].first;
                p2 += *path[i].second;
            } else {
                const VertexRecord &vr = reads.getRecord(path.getVertex(last));
                Sequence c1 = path[i].first->firstNucl();
                Sequence c2 = path[i].second->firstNucl();
                size_t score1 = vr.countStartsWith(s1 + seq + c1) + vr.countStartsWith(s2 + seq + c2);
                size_t score2 = vr.countStartsWith(s1 + seq + c2) + vr.countStartsWith(s2 + seq + c1);
                if(score1 >= score2) {
                    p1 += *path[i].first;
                    p2 += *path[i].second;
                } else {
                    p1 += *path[i].second;
                    p2 += *path[i].first;
                }
            }
            s1 = p1.back().firstNucl();
            s2 = p2.back().firstNucl();
            last = i;
            seq = Sequence();
        } else {
            p1 += *path[i].first;
            p2 += *path[i].first;
            seq = seq + path[i].first->firstNucl();
        }
    }
    return {std::move(p1), std::move(p2)};
}

std::string BulgePathCorrector::correctRead(dbg::GraphAlignment &path) {
    std::vector<Case> cases;
    std::vector<std::string> messages;
    for(size_t i = 0; i < path.size(); i++) {
        if(!cases.empty() && cases.back().read_to == i && cases.back().path_to < paths[cases.back().path_ind].size() &&
           (paths[cases.back().path_ind][cases.back().path_to].first == &path[i].contig() || paths[cases.back().path_ind][cases.back().path_to].second == &path[i].contig())) {
            if(resolved[cases.back().path_ind].first[cases.back().path_to] != path[i].contig())
                cases.back().score1 += 1;
            if(resolved[cases.back().path_ind].second[cases.back().path_to] != path[i].contig())
                cases.back().score2 += 1;
            cases.back().read_to += 1;
            cases.back().path_to += 1;
            continue;
        }
        auto it = pathPoses.find(&path[i].contig());
        if(it == pathPoses.end()) {
            continue;
        }
        size_t score1 = (resolved[it->second.path_ind].first[it->second.pos] != path[i].contig());
        size_t score2 = (resolved[it->second.path_ind].second[it->second.pos] != path[i].contig());
        cases.emplace_back(it->second.path_ind, it->second.pos, it->second.pos + 1, i, i + 1, score1, score2);
    }
    if(cases.empty())
        return "";
    dbg::GraphAlignment res;
    for(Case & bp : cases) {
        for(size_t i = res.size(); i < bp.read_from; i++)
            res += path[i];
        if(bp.score1 == 0 || bp.score2 == 0) {
            for(size_t i = res.size(); i < bp.read_to; i++)
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
    for(size_t i = res.size(); i < path.size(); i++)
        res += path[i];
    if(res.front() != path.front() && res.front().size() > path.front().size()) {
        res.front() = res.front().shrinkLeft(res.front().size() - path.front().size());
    }
    if(res.back() != path.back() && res.back().size() > path.back().size()) {
        res.back() = res.back().shrinkRight(res.back().size() - path.back().size());
    }
    path = std::move(res);
    return join("_", messages);
}

void BulgePathCorrector::initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads) {
    paths = BulgePathFinder(dbg, threshold).paths;
    for(BulgePath &path : BulgePathFinder(dbg, threshold).paths) {
        if(path.size() > 1 && path.length() > unique_length) {
            paths.emplace_back(path);
        }
    }
    for(size_t path_ind = 0; path_ind < paths.size(); path_ind++) {
        BulgePath &path = paths[path_ind];
        for(size_t i = 0; i < path.size(); i++) {
            std::pair<dbg::Edge *, dbg::Edge *> pair = path[i];
            pathPoses[pair.first] = {path_ind, i};
            pathPoses[pair.second] = {path_ind, i};
        }
        resolved.emplace_back(resolveBulgePath(reads, path));
    }
}

