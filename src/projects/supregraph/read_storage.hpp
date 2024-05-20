#pragma once
#include "list_path.hpp"
#include "vertex_resolution.hpp"
#include "supregraph.hpp"
#include "unique_vertex_storage.hpp"

namespace spg {

    typedef typename ag::GraphPath<SPGTraits> GraphPath;
    typedef typename ag::PathIterator<SPGTraits> PathIterator;
    typedef typename ag::PathDirection<SPGTraits> PathDirection;
    class ReadDirection;

    class ReadRecord {
    public:
        std::string name;
        GraphPath path;
        ReadDirection forward();
        ReadDirection backward();
        explicit ReadRecord(std::string name, GraphPath path = {}) : name(std::move(name)), path(std::move(path)) {
        }
    };

    class ReadDirection : public PathDirection {
    private:
        ReadRecord *read;
    public:
        ReadDirection(ReadRecord &read, bool rc) : PathDirection(read.path, rc), read(&read) {
        }
        ReadRecord &getRead() const {return *read;}
    };

    inline std::ostream &operator<<(std::ostream &os, const ReadDirection &dir) {
        os << dir.getRead().name << "_" << (dir.isForward() ? "F" : "B") << ":";
        bool first = true;
        for(Edge &edge : dir) {
            if(!first)
                os << ",";
            else
                first = false;
            os << edge.getId();
        }
        return os;
    }


    class PathStorage {
    private:
        SupreGraph *spg;
        std::vector<ReadRecord> reads;

        void shrink(GraphPath &path) const {
            while(!path.empty() && path.frontEdge().isSuffix() && path.start().size() <= path.cutLeft()) {
                path.pop_front();
            }
            while(!path.empty() && path.backEdge().isPrefix() && path.finish().size() <= path.cutRight())
                path.pop_back();
        }
    public:

        explicit PathStorage(SupreGraph &spg);
        PathStorage(PathStorage &&other) noexcept = default;
//        PathStorage(PathStorage &&other) noexcept : ResolutionListener(*other.spg), spg(other.spg), read_index(std::move(other.read_index)),
//                    outgoing_index(std::move(other.outgoing_index)), inner_index(std::move(other.inner_index)),
//                    reads(std::move(other.reads)) {
//
//        }
        PathStorage(const PathStorage &other) = delete;

        std::vector<ReadRecord>::iterator begin() {return reads.begin();}
        std::vector<ReadRecord>::iterator end() {return reads.end();}
        std::vector<ReadRecord>::const_iterator begin() const {return reads.begin();}
        std::vector<ReadRecord>::const_iterator end() const {return reads.end();}

        ReadRecord &operator[](size_t ind) {return reads[ind];}
        const ReadRecord &operator[](size_t ind) const {return reads[ind];}

//        auto passingReads(Vertex &v) const {
//            for(auto &it : getReadPositions(inc.rc())) {
//                auto pos = it.second;
//                auto &path = storage[it.first].path;
//                if(pos == path.path.begin())
//                    continue;
//                auto next = pos;
//                --next;
//                res.add(inc, (**next).rc());
//                has_covering = true;
//            }
//        }


        ReadRecord &addRead(std::string name, GraphPath path);

    };

    class PathIndex : ResolutionListener {
    private:
        std::unordered_map<EdgeId, std::vector<std::pair<ReadDirection, PathIterator>>> read_index;
        std::unordered_map<VertexId, std::vector<ReadDirection>> outgoing_index;
        std::unordered_map<VertexId, std::vector<Segment<Vertex>>> inner_index;

        void processPassing(Edge &edge, const spg::VertexResolutionResult &resolution);
    public:
        struct PassingRead {
            EdgePair edges;
            ReadDirection direction;
            PathIterator position;

            PassingRead(const ReadDirection &direction, const PathIterator &position) :
                    edges(*position, *(position + 1)), direction(direction), position(position) {}
        };

        PathIndex(SupreGraph &spg, PathStorage &storage);

        ComplexIterableStorage<Generator<std::vector<std::pair<ReadDirection, PathIterator>>::iterator, PassingRead>> getPassing(Vertex &v) &;
        const std::vector<std::pair<spg::ReadDirection, PathIterator>> &getReadPositions(Edge &edge) const;
        const std::vector<spg::ReadDirection> &getOutgoingReads(Vertex &v) const;
        const std::vector<Segment<Vertex>> &getInnerReads(Vertex &v) const;

        bool checkReadIndexConsistency() const;

        void fireAddVertex(Vertex &vertex);
        void fireAddEdge(Edge &edge);
        void fireDeleteVertex(Vertex &vertex) override;
        void fireDeleteEdge(Edge &edge) override;

        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;

    };
}