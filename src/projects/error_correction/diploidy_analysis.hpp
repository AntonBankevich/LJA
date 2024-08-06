#pragma once

#include "uniqueness.hpp"
#include "dbg/sparse_dbg.hpp"
#include "error_correction.hpp"

template<class Traits>
class BulgePath {
public:
    typedef typename Traits::Edge::EdgeId EdgeId;
    typedef typename Traits::Vertex::VertexId VertexId;

private:
    typedef typename std::vector<std::pair<EdgeId, EdgeId>> storage_type;
    typedef typename storage_type::const_iterator iterator_type;
    VertexId start_;
    storage_type path;

public:
    typedef typename Traits::Vertex Vertex;
    typedef typename Traits::Edge Edge;

    explicit BulgePath(Vertex &_start) : start_(_start.getId()) {}

    explicit BulgePath(Edge &edge) : start_(edge.getStart().getId()) {
        path.emplace_back(edge.getId(), edge.getId());
    }

    explicit BulgePath(std::vector<std::pair<EdgeId, EdgeId>> &&path_) :
        start_({}), path(path_) {
        VERIFY(path.size() > 0);
        start_ = path.front().first->getStart().getId();
    }

    Vertex &finish() const{
        if (path.empty())
            return *start_;
        return path.back().first->getFinish();
    }

    Vertex &start() const{
        return *start_;
    }

    Vertex &getVertex(size_t ind) const{
        VERIFY(ind <= size());
        if (ind == size())
            return finish();
        return path[ind].first->getStart();
    }

    const std::pair<EdgeId, EdgeId> &operator[](size_t ind) const {
        return path[ind];
    }

    iterator_type begin() const { return path.begin(); }

    iterator_type end() const { return path.end(); }

    Vertex &vertexAt(size_t ind){
        if (ind == 0)
            return *start_;
        return path[ind - 1].first->getFinish();
    }

    BulgePath RC(){
        if (path.empty()) {
            return BulgePath(start_->rc());
        }
        storage_type rc;
        for (size_t i = 0; i < path.size(); i++) {
            rc.emplace_back(path[path.size() - 1 - i].first->rc().getId(), path[path.size() - 1 - i].second->rc().getId());
        }
        return BulgePath(std::move(rc));
    }

    void extend(double threshold){
        Vertex &last = finish();
        size_t deg = last.outDeg();
        if (last.front().getCoverage() > threshold && last.back().getCoverage() > threshold)
            path.emplace_back(last.front().getId(), last.back().getId());
        else {
            for (Edge &edge: last) {
                if (edge.getCoverage() > threshold) {
                    path.emplace_back(edge.getId(), edge.getId());
                    return;
                }
            }
            VERIFY(last.outDeg() == 2 && last.front().getFinish() == last.back().getFinish());
            path.emplace_back(last.front().getId(), last.back().getId());
        }
    }

    BulgePath operator+(const BulgePath &other) const{
        VERIFY(finish() == other.start())
        storage_type sum(path);
        sum.insert(sum.end(), other.path.begin(), other.path.end());
        return BulgePath(std::move(sum));
    }

    ag::GraphPath<Traits> randomPath() const{
        ag::GraphPath<Traits> res(start());
        for (const std::pair<EdgeId, EdgeId> &pair: path) {
            res += *pair.first;
        }
        return std::move(res);
    }

    bool isBulge(size_t ind) const { return path[ind].first != path[ind].second; }

    // QQ
    bool isBad(size_t bad_bulge_inner_size) const{
        if (path.size() < 2)
            return false;
        for (const auto &p: path) {
            if (p.first != p.second)
                if (p.first->innerSize() > bad_bulge_inner_size || p.second->truncSize() > bad_bulge_inner_size) {
                    return false;
                }
        }
        return true;
    }

    size_t size() const { return path.size(); }

    size_t length() const{
        size_t res = 0;
        for (auto &p: path) {
            res += std::max(p.first->truncSize(), p.second->truncSize());
        }
        return res;
    }

    size_t bulgeLength() const{
        size_t res = 0;
        for (auto &p: path) {
            if (p.first != p.second)
                res += std::max(p.first->truncSize(), p.second->truncSize());
        }
        return res;
    }

    size_t conservativeLength() const{
        size_t res = 0;
        for (auto &p: path) {
            if (p.first == p.second)
                res += std::max(p.first->truncSize(), p.second->truncSize());
        }
        return res;
    }

    std::string str() const{
        std::stringstream ss;
        ss << start().getId();
        for (const auto &p: path) {
            if (p.first == p.second) {
                ss << "-" << p.first->truncSize() << p.first->nuclLabel() << "-" << p.first->getFinish().getId();
            } else {
                ss << "-(" << p.first->truncSize() << p.first->nuclLabel() << "," <<
                   p.second->truncSize() << p.second->nuclLabel() << ")-" << p.first->getFinish().getId ();
            }
        }
        return ss.str();
    }
};

template<class Traits>
class BulgePathFinder {
public:
    typedef typename Traits::Edge::EdgeId EdgeId;
    typedef typename Traits::Vertex::VertexId VertexId;
    typedef typename Traits::Vertex Vertex;
    typedef typename Traits::Edge Edge;
private:
    bool checkVertexForward(const ag::BaseVertex<Traits> &v){
        size_t outgoing_edge_cnt = 0;
        for (Edge &edge: v) {
            if (edge.getCoverage() > threshold)
                outgoing_edge_cnt++;
        }
        //if (outgoing_edge_cnt > 1 && outgoing_edge_cnt != v.outDeg()) {return false;}
        if (outgoing_edge_cnt == 0) {return false;}
        if (outgoing_edge_cnt == 1) {return true;}
        if (v.outDeg() > 2) {return false;}
        if (v.front().getFinish() != v.back().getFinish()) {return false;}
        if (v.front().truncSize() > v.back().truncSize() * 1.3 ||
            v.back().truncSize() > v.front().truncSize() * 1.3) {return false;}
        return true;
    }

    bool isBulgePathInner(const ag::BaseVertex<Traits> &v) {
        return checkVertexForward(v) && checkVertexForward(v.rc());
    }

    BulgePath<Traits> forwardPath(Vertex &start){
        BulgePath<Traits> res(start);
        ag::BaseVertex<Traits> *cur = &start;
        while (isBulgePathInner(*cur)) {
            res.extend(threshold);
            cur = &res.finish();
            if (cur == &start)
                return std::move(res);
        }
        return std::move(res);
    }

    ag::AssemblyGraph<Traits> &dbg;
    double threshold;
public:
    std::vector<BulgePath<Traits>> paths;

    explicit BulgePathFinder(ag::AssemblyGraph<Traits> &dbg, double threshold = 0) :
    dbg(dbg), threshold(threshold) {
        std::unordered_set<ag::BaseVertex<Traits> *> visited;
        for (auto &vertex: dbg.verticesUnique()) {
            if (visited.find(&vertex) != visited.end())
                continue;
            if (isBulgePathInner(vertex)) {
                BulgePath<Traits> new_path = forwardPath(vertex);
                VERIFY(new_path.size() > 0);
                if (new_path.start() != new_path.finish()) {
                    BulgePath p2 = forwardPath(vertex.rc());
                    BulgePath p3 = p2.RC();
                    new_path = p3 + new_path;
                }
                for (size_t i = 1; i + 1 <= new_path.size(); i++) {
                    visited.emplace(&new_path.vertexAt(i));
                    visited.emplace(&new_path.vertexAt(i).rc());
                }
                if (new_path.start() == new_path.finish()) {
                    visited.emplace(&new_path.start());
                    visited.emplace(&new_path.start().rc());
                }
                //paths.emplace_back(new_path.RC());
                paths.emplace_back(std::move(new_path));
            }
        }
        /*
        for (Edge &edge: dbg.edges()) {
            if (visited.find(&edge.getFinish()) == visited.end() && visited.find(&edge.getStart()) == visited.end())
                paths.emplace_back(edge);
        }*/
    }

    SetUniquenessStorage uniqueEdges(size_t min_len) const{
        std::vector<EdgeId> res;
        for (const BulgePath<Traits> &bp: paths) {
            if (bp.size() == 1) {
                Edge &edge = *bp[0].first;
                if (edge.truncSize() > min_len || (
                        (edge.getStart().inDeg() == 0 || edge.getFinish().outDeg() == 0) &&
                        edge.truncSize() > min_len / 3 && edge.getCoverage() > 4)) {
                    res.emplace_back(edge.getId());
                        }
            } else {
                if (!bp.isBad(0) && (bp.conservativeLength() < bp.length() / 2 || bp.size() >= 4) &&
                    bp.length() > min_len && bp.conservativeLength() < bp.length() * 95 / 100) {
                    for (auto &p: bp) {
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
};


class BulgePathCorrector : public AbstractCorrectionAlgorithm {
private:
    struct PathPos {
        size_t path_ind;
        size_t pos;
    };

    struct Case {
        Case(size_t pathInd, size_t pathFrom, size_t pathTo, size_t readFrom, size_t readTo, size_t score1,
             size_t score2) : path_ind(pathInd),
                              path_from(pathFrom),
                              path_to(pathTo),
                              read_from(readFrom),
                              read_to(readTo),
                              score1(score1),
                              score2(score2) {}

        size_t path_ind;
        size_t path_from;
        size_t path_to;
        size_t read_from;
        size_t read_to;
        size_t score1;
        size_t score2;
    };

    std::pair<dbg::GraphPath, dbg::GraphPath>
    resolveBulgePath(const dbg::ReadAlignmentStorage &reads, const BulgePath<dbg::DBGTraits> &path) const {
        dbg::GraphPath p1(path.start());
        dbg::GraphPath p2(path.start());
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
                    const ag::VertexRecord<dbg::DBGTraits> &vr = reads.getRecord(path.getVertex(last));
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
                s1 = p1.backEdge().firstNucl();
                s2 = p2.backEdge().firstNucl();
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

    std::vector<BulgePath<dbg::DBGTraits>> paths;
    std::vector<std::pair<dbg::GraphPath, dbg::GraphPath>> resolved;
    std::unordered_map<dbg::EdgeId, PathPos> pathPoses;
    size_t unique_length;
    double threshold;
public:
    BulgePathCorrector(dbg::SparseDBG &dbg, dbg::ReadAlignmentStorage &reads, size_t unique_length,
                       double threshold) : AbstractCorrectionAlgorithm("BulgePathFixer"),
                                           unique_length(unique_length), threshold(threshold) {
    }

    void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                    dbg::ReadAlignmentStorage &reads) override {
        paths = BulgePathFinder(dbg, threshold).paths;
        for(BulgePath<dbg::DBGTraits> &path : BulgePathFinder(dbg, threshold).paths) {
            if(path.size() > 1 && path.length() > unique_length) {
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

    std::string correctRead(dbg::GraphPath &path) override{
        std::vector<Case> cases;
        std::vector<std::string> messages;
        for(size_t i = 0; i < path.size(); i++) {
            if(!cases.empty() && cases.back().read_to == i && cases.back().path_to < paths[cases.back().path_ind].size() &&
               (*paths[cases.back().path_ind][cases.back().path_to].first == path[i].contig() ||
                *paths[cases.back().path_ind][cases.back().path_to].second == path[i].contig())) {
                if(resolved[cases.back().path_ind].first[cases.back().path_to] != path[i].contig())
                    cases.back().score1 += 1;
                if(resolved[cases.back().path_ind].second[cases.back().path_to] != path[i].contig())
                    cases.back().score2 += 1;
                cases.back().read_to += 1;
                cases.back().path_to += 1;
                continue;
            }
            auto it = pathPoses.find(path[i].contig().getId());
            if(it == pathPoses.end()) {
                continue;
            }
            size_t score1 = (resolved[it->second.path_ind].first[it->second.pos] != path[i].contig());
            size_t score2 = (resolved[it->second.path_ind].second[it->second.pos] != path[i].contig());
            cases.emplace_back(it->second.path_ind, it->second.pos, it->second.pos + 1, i, i + 1, score1, score2);
        }
        if(cases.empty())
            return "";
        dbg::GraphPath res;
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
            res.front() = res.front().shrinkLeftBy(res.front().size() - path.front().size());
        }
        if(res.back() != path.back() && res.back().size() > path.back().size()) {
            res.back() = res.back().shrinkRightBy(res.back().size() - path.back().size());
        }
        path = std::move(res);
        return join("_", messages);
    }
};

