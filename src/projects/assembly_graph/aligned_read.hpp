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

    template<class Traits>
    class AlignedReadDirection;

    template<class Traits>
    class AlignedRead {
    private:
        std::string id = {};
        GraphPath <Traits> path = {};
        GraphPath <Traits> corrected_path = {};
        omp_lock_t writelock = {};
        bool corrected = false;
    public:
        AlignedRead() = default;

        AlignedRead(AlignedRead &&other) noexcept = default;

        AlignedRead &operator=(AlignedRead &&other) noexcept = default;

        explicit AlignedRead(std::string readId) : id(std::move(readId)), path(), corrected(false) {}

        AlignedRead(std::string readId, const ag::GraphPath<Traits> &_path) : id(std::move(readId)), path(_path),
                                                                              corrected(false) {}

        void lock() { omp_set_lock(&writelock); }

        void unlock() { omp_unset_lock(&writelock); }

        bool operator<(const AlignedRead &other) const { return id < other.id; }

        AlignedReadDirection<Traits> forward();

        AlignedReadDirection<Traits> backward();

        const std::string &getId() const { return id; }

        GraphPath <Traits> &getPath() { return path; }

        const GraphPath <Traits> &getPath() const { return path; }

        GraphPath <Traits> &getCorrected() { return corrected_path; }

        const GraphPath <Traits> &getCorrected() const { return corrected_path; }

        bool checkCorrected() const { return corrected; }

        bool valid() const { return path.valid(); }

        void delayedInvalidate();

        void correct(GraphPath <Traits> cpath);

        void resetEdgeCodes();

        void setPath(GraphPath <Traits> new_path);

        void applyCorrection();

        static AlignedRead Load(std::istream &is, const IdIndex<typename Traits::Vertex> &index) {
            std::string id;
            is >> id;
            id = decodeReadId(id);
            return {id, GraphPath<Traits>::Load(is, index)};
        }
    };

    template<class Traits>
    inline std::ostream &operator<<(std::ostream &os, const AlignedRead<Traits> &alignedRead) {
        return os << encodeReadId(alignedRead.getId()) << " " << alignedRead.getPath();
    }

    template<class Traits>
    void AlignedRead<Traits>::applyCorrection() {
        if (corrected)
            setPath(std::move(corrected_path));
        corrected_path = {};
        corrected = false;
    }

    template<class Traits>
    void AlignedRead<Traits>::delayedInvalidate() {
        VERIFY(!corrected);
        corrected_path = {};
        corrected = true;
    }

    template<class Traits>
    class AlignedReadDirection : public PathDirection<Traits> {
    private:
        AlignedRead<Traits> *read;
    public:
        AlignedReadDirection() = default;

        AlignedReadDirection(AlignedRead<Traits> &read, bool rc) : PathDirection<Traits>(read.getPath(), rc),
                                                                   read(&read) {
        }

        AlignedRead<Traits> &getRead() const { return *read; }

        void setPath(GraphPath <Traits> path) const {
            this->isRC() ? read->setPath(path.RC()) : read->setPath(std::move(path));
        }

        AlignedReadDirection RC() const { return {*read, !this->isRC()}; }

        PathDirection <Traits> getCorrected() { return {read->getCorrected(), this->isRC()}; }
    };

    template<class Traits>
    AlignedReadDirection<Traits> AlignedRead<Traits>::forward() { return {*this, false}; }

    template<class Traits>
    AlignedReadDirection<Traits> AlignedRead<Traits>::backward() { return {*this, true}; }

    template<class Traits>
    void AlignedRead<Traits>::setPath(GraphPath <Traits> new_path) {
        path = std::move(new_path);
    }

    template<class Traits>
    void AlignedRead<Traits>::correct(GraphPath <Traits> cpath) {
        VERIFY_MSG(!corrected, "Attempt to correct path while previous correction was not yet applied");
        corrected_path = std::move(cpath);
        corrected = true;
    }

    template<class Traits>
    void AlignedRead<Traits>::resetEdgeCodes() {
        path.resetEdgeCodes();
    }

    template<class Traits>
    inline std::ostream &operator<<(std::ostream &os, const AlignedReadDirection<Traits> &dir) {
        os << dir.getRead().name << "_" << (dir.isForward() ? "F" : "B") << ":";
        bool first = true;
        for (auto &edge: dir) {
            if (!first)
                os << ",";
            else
                first = false;
            os << edge.getId();
        }
        return os;
    }
}