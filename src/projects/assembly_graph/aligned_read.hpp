#pragma once
#include "assembly_graph_base.hpp"
#include "graph_paths.hpp"

namespace ag {
    //TODO: get rid of read names entirely and replace with integer ids
    inline std::string encodeReadId(const std::string &str) {
        std::string res = str;
        std::replace(res.begin(), res.end(), ' ', '^');
        std::replace(res.begin(), res.end(), '\t', '$');
        return std::move(res);
    }

    inline std::string decodeReadId(const std::string &str) {
        std::string res = str;
        std::replace(res.begin(), res.end(), '^', ' ');
        std::replace(res.begin(), res.end(), '$', '\t');
        return std::move(res);
    }

    class AlignedReadDirection;

    class AlignedRead {
    private:
        std::string id = {};
        GraphPath  path = {};
        GraphPath  corrected_path = {};
        omp_lock_t writelock = {};
        bool corrected = false;
    public:
        AlignedRead() = default;
        AlignedRead(AlignedRead &&other) noexcept = default;
        AlignedRead &operator=(AlignedRead &&other) noexcept = default;
        explicit AlignedRead(std::string readId) : id(std::move(readId)), path(), corrected(false) {}
        AlignedRead(std::string readId, const ag::GraphPath &_path) : id(std::move(readId)), path(_path),
                                                                              corrected(false) {}

        void lock() { omp_set_lock(&writelock); }

        void unlock() { omp_unset_lock(&writelock); }

        bool operator<(const AlignedRead &other) const { return id < other.id; }

        AlignedReadDirection forward();

        AlignedReadDirection backward();

        const std::string &getId() const { return id; }

//        Hidden contract: a stored path is always the shortest path, by vertex count, that spells the
//        read's sequence — so it can never start with a prefix edge or end with a suffix edge (trimming
//        those loses no sequence but shortens the path). GraphPath::normalize() is what enforces this.
        GraphPath  &getPath() { return path; }

        const GraphPath  &getPath() const { return path; }

        GraphPath  &getCorrected() { return corrected_path; }

        const GraphPath  &getCorrected() const { return corrected_path; }

        bool checkCorrected() const { return corrected; }

        bool valid() const { return path.valid(); }

        void delayedInvalidate();

        void correct(GraphPath  cpath);

        void resetEdgeCodes();

        void setPath(GraphPath  new_path);

        void applyCorrection();

        static AlignedRead Load(std::istream &is, const IdIndex<Vertex> &index) {
            std::string id;
            is >> id;
            id = decodeReadId(id);
            return {id, GraphPath::Load(is, index)};
        }
    };


    inline std::ostream &operator<<(std::ostream &os, const AlignedRead &alignedRead) {
        return os << encodeReadId(alignedRead.getId()) << " " << alignedRead.getPath();
    }


    class AlignedReadDirection : public PathDirection {
    private:
        AlignedRead *read;
    public:
        AlignedReadDirection() = default;

        AlignedReadDirection(AlignedRead &read, bool rc) : PathDirection(read.getPath(), rc),
                                                                   read(&read) {
        }

        AlignedRead &getRead() const { return *read; }

        void setPath(GraphPath  path) const {
            this->isRC() ? read->setPath(path.RC()) : read->setPath(std::move(path));
        }

        AlignedReadDirection RC() const { return {*read, !this->isRC()}; }

        PathDirection  getCorrected() { return {read->getCorrected(), this->isRC()}; }
    };
}