#pragma once

#include "uniqueness.hpp"
#include "dbg/sparse_dbg.hpp"
#include "error_correction.hpp"

class BulgePath {
private:
    typedef typename std::vector<std::pair<dbg::Edge *, dbg::Edge *>> storage_type;
    typedef typename storage_type::const_iterator iterator_type;
    dbg::Vertex *start_;
    std::vector<std::pair<dbg::Edge *, dbg::Edge *>> path;
public:
    explicit BulgePath(dbg::Vertex &_start) : start_(&_start) {}
    explicit BulgePath(dbg::Edge &edge) : start_(edge.start()) {path.emplace_back(&edge, &edge);}
    explicit BulgePath(std::vector<std::pair<dbg::Edge *, dbg::Edge *>> &&path_);

    dbg::Vertex &finish() const;
    dbg::Vertex &start() const;
    dbg::Vertex &getVertex(size_t ind) const;
    const std::pair<dbg::Edge *, dbg::Edge *> &operator[](size_t ind) const {return path[ind];}
    iterator_type begin() const {return path.begin();}
    iterator_type end() const {return path.end();}
    dbg::Vertex &vertexAt(size_t ind);

    BulgePath RC();
    void extend(double threshold);
    BulgePath operator+(const BulgePath &other) const;
    dbg::Path randomPath() const;

    bool isBulge(size_t ind) const {return path[ind].first != path[ind].second;}
    bool isBad(size_t bad_bulge_size) const;
    size_t size() const {return path.size();}
    size_t length() const;
    size_t bulgeLength() const;
    size_t conservativeLength() const;
    std::string str() const;
};


class BulgePathFinder {
private:
    bool checkVertexForward(const dbg::Vertex &v);
    bool isInner(const dbg::Vertex &v);
    BulgePath forwardPath(dbg::Vertex &start);

    dbg::SparseDBG &dbg;
    double threshold;
public:
    std::vector<BulgePath> paths;

    explicit BulgePathFinder(dbg::SparseDBG &dbg, double threshold = 0);
    SetUniquenessStorage uniqueEdges(size_t min_len) const;
};

class BulgePathCorrector : public AbstractCorrectionAlgorithm {
private:
    struct PathPos {
        size_t path_ind;
        size_t pos;
    };

    struct Case {
        Case(size_t pathInd, size_t pathFrom, size_t pathTo, size_t readFrom, size_t readTo, size_t score1, size_t score2) : path_ind(pathInd),
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
    std::pair<dbg::Path, dbg::Path> resolveBulgePath(const RecordStorage &reads, const BulgePath &path) const;

    std::vector<BulgePath> paths;
    std::vector<std::pair<dbg::Path, dbg::Path>> resolved;
    std::unordered_map<dbg::Edge *, PathPos> pathPoses;
    size_t unique_length;
    double threshold;
public:
    BulgePathCorrector(dbg::SparseDBG &dbg, RecordStorage &reads, size_t unique_length,
                                           double threshold) : AbstractCorrectionAlgorithm("BulgePathFixer"),
                                           unique_length(unique_length), threshold(threshold) {
    }

    void initialize(logging::Logger &logger, size_t threads, dbg::SparseDBG &dbg, RecordStorage &reads) override;

    std::string correctRead(dbg::GraphAlignment &path) override;
};
