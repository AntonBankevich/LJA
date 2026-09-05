#pragma once

#include "uniqueness.hpp"
#include "dbg/sparse_dbg.hpp"
#include "dbg/dbg_read_alignment_storage.hpp"
#include "dbg_correction_algorithm.hpp"

namespace ag {
    class BulgePath {
    private:
        typedef typename std::vector<std::pair<EdgeId, EdgeId>> storage_type;
        typedef typename storage_type::const_iterator iterator_type;
        VertexId start_;
        storage_type path;

    public:
        explicit BulgePath(Vertex &_start) : start_(_start.getId()) {}

        explicit BulgePath(Edge &edge) : start_(edge.getStart().getId()) {
            path.emplace_back(edge.getId(), edge.getId());
        }

        explicit BulgePath(std::vector<std::pair<EdgeId, EdgeId>> &&path_) :
                start_({}), path(path_) {
            VERIFY(path.size() > 0);
            start_ = path.front().first->getStart().getId();
        }

        Vertex &getFinish() const {
            if (path.empty())
                return *start_;
            return path.back().first->getFinish();
        }

        Vertex &getStart() const {
            return *start_;
        }

        Vertex &getVertex(size_t ind) const {
            VERIFY(ind <= size());
            if (ind == size())
                return getFinish();
            return path[ind].first->getStart();
        }

        const std::pair<EdgeId, EdgeId> &operator[](size_t ind) const {
            return path[ind];
        }

        iterator_type begin() const { return path.begin(); }

        iterator_type end() const { return path.end(); }

        Vertex &vertexAt(size_t ind) {
            if (ind == 0)
                return *start_;
            return path[ind - 1].first->getFinish();
        }

        BulgePath RC() {
            if (path.empty()) {
                return BulgePath(start_->rc());
            }
            storage_type rc;
            for (size_t i = 0; i < path.size(); i++) {
                rc.emplace_back(path[path.size() - 1 - i].first->rc().getId(),
                                path[path.size() - 1 - i].second->rc().getId());
            }
            return BulgePath(std::move(rc));
        }

        void extend(double threshold) {
            Vertex &last = getFinish();
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
                VERIFY(last.outDeg() == 2 && last.frontVertex() == last.backVertex());
                path.emplace_back(last.front().getId(), last.back().getId());
            }
        }

        BulgePath operator+(const BulgePath &other) const {
            VERIFY(getFinish() == other.getStart())
            storage_type sum(path);
            sum.insert(sum.end(), other.path.begin(), other.path.end());
            return BulgePath(std::move(sum));
        }

        ag::GraphPath randomPath() const {
            ag::GraphPath res(getStart());
            for (const std::pair<EdgeId, EdgeId> &pair: path) {
                res += *pair.first;
            }
            return std::move(res);
        }

        bool isBulge(size_t ind) const { return path[ind].first != path[ind].second; }

        // QQ
        bool isBad(size_t bad_bulge_inner_size) const {
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

        size_t length() const {
            size_t res = 0;
            for (auto &p: path) {
                res += std::max(p.first->truncSize(), p.second->truncSize());
            }
            return res;
        }

        size_t bulgeLength() const {
            size_t res = 0;
            for (auto &p: path) {
                if (p.first != p.second)
                    res += std::max(p.first->truncSize(), p.second->truncSize());
            }
            return res;
        }

        size_t conservativeLength() const {
            size_t res = 0;
            for (auto &p: path) {
                if (p.first == p.second)
                    res += std::max(p.first->truncSize(), p.second->truncSize());
            }
            return res;
        }

        std::string str() const {
            std::stringstream ss;
            ss << getStart().getId();
            for (const auto &p: path) {
                if (p.first == p.second) {
                    ss << "-" << p.first->truncSize() << p.first->getCode() << "-" << p.first->getFinish().getId();
                } else {
                    ss << "-(" << p.first->truncSize() << p.first->getCode() << "," <<
                       p.second->truncSize() << p.second->getCode() << ")-" << p.first->getFinish().getId();
                }
            }
            return ss.str();
        }
    };

    class BulgePathFinder {
    private:
        bool checkVertexForward(const Vertex &v){
            size_t outgoing_edge_cnt = 0;
            for (Edge &edge: v) {
                if (edge.getCoverage() > threshold)
                    outgoing_edge_cnt++;
            }
            //if (outgoing_edge_cnt > 1 && outgoing_edge_cnt != v.outDeg()) {return false;}
            if (outgoing_edge_cnt == 0) {return false;}
            if (outgoing_edge_cnt == 1) {return true;}
            if (v.outDeg() > 2) {return false;}
            if (v.frontVertex() != v.backVertex()) {return false;}
            if (v.front().truncSize() > v.back().truncSize() * 1.3 ||
                v.back().truncSize() > v.front().truncSize() * 1.3) {return false;}
            return true;
        }

        bool isBulgePathInner(const Vertex &v) {
            return checkVertexForward(v) && checkVertexForward(v.rc());
        }

        BulgePath forwardPath(Vertex &start);

        ag::AssemblyGraph &dbg;
        double threshold;
    public:
        std::vector<BulgePath> paths;

        explicit BulgePathFinder(ag::AssemblyGraph &dbg, double threshold = 0);

        SetUniquenessStorage uniqueEdges(size_t min_len) const;
    };
}

namespace dbg {
    class BulgePathCorrector : public dbg::AbstractDBGCorrectionAlgorithm {
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

        std::pair<ag::RAGraphPath, ag::RAGraphPath>
        resolveBulgePath(const dbg::DBGAlignedReadStorage &reads, const ag::BulgePath &path) const;

        std::vector<ag::BulgePath> paths;
        std::vector<std::pair<ag::RAGraphPath, ag::RAGraphPath>> resolved;
        std::unordered_map<dbg::EdgeId, PathPos> pathPoses;
        size_t unique_length;
        double threshold;
    public:
        BulgePathCorrector(dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &reads, size_t unique_length,
                           double threshold) : dbg::AbstractDBGCorrectionAlgorithm("BulgePathFixer"),
                                               unique_length(unique_length), threshold(threshold) {
        }

        void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg,
                        dbg::DBGAlignedReadStorage &reads) override;

        std::string correctRead(const std::string &name, ag::GraphPath &path) override;
    };
}