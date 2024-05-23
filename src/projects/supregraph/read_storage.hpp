#pragma once
#include "list_path.hpp"
#include "vertex_resolution.hpp"
#include "supregraph.hpp"
#include "unique_vertex_storage.hpp"
#include "listeners.hpp"

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
        ReadDirection RC() const {return {*read, !isRC()};}
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

    template<class ContigId>
    class IdSegment {
    public:
        ContigId id;
        size_t cut_left;
        size_t cut_right;
        IdSegment() : id(), cut_left(0), cut_right(0) {};
        IdSegment(ContigId id, size_t cut_left, size_t cut_right) : id(id), cut_left(cut_left), cut_right(cut_right) {}
        Segment<typename ContigId::base> asSegment() {return {*id, cut_left, id->truncSize() - cut_right};}
        IdSegment nest(const IdSegment &other) {return {other.id, other.cut_left + cut_left, other.cut_right + cut_right};}
        size_t truncSize() {return id->truncSize() - cut_left - cut_right;}
    };

    class Embedding {
    private:
        std::unordered_map<EdgeId, IdSegment<EdgeId>> edge_embedding;
        std::unordered_map<VertexId, IdSegment<VertexId>> vertex_embedding;

        Segment<Edge> embed(EdgeId eid, size_t cut_left = 0, size_t cut_right = 0) const;
        Segment<Vertex> embed(VertexId vid, size_t cut_left = 0, size_t cut_right = 0) const;
    public:
        void remap(const GraphPath &path, Vertex &v);
        GraphPath embed(const GraphPath &path);
        void clear() {edge_embedding.clear(); vertex_embedding.clear();}
    };


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
        ReadRecord &Oppa();

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
//        When reads_ready==false these two, as well as the path storage itself are invalidated
        mutable std::unordered_map<EdgeId, std::vector<std::pair<ReadDirection, PathIterator>>> read_index;
        mutable std::unordered_map<EdgeId, std::vector<ReadDirection>> outgoing_index;

//        inner index always remains valid
        mutable std::unordered_map<VertexId, std::vector<Segment<Vertex>>> inner_index;
//        remapping is used only when reads_ready=false and allows to reconstruct read paths when necessary
        mutable Embedding embedding;
        PathStorage *storage;
        mutable bool reads_ready = false;

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

        void unprepare() const {
            if(reads_ready) {
                for (auto &it: read_index)
                    it.second.clear();
                reads_ready = false;
            }
        }

        void restorePaths() const {
            for (ReadRecord &rr: *this->storage) {
                rr.path = embedding.embed(rr.path);
            }
        }
        void prepareIndex() const;
        void prepare() const {
            if(!reads_ready) {
                restorePaths();
                prepareIndex();
                embedding.clear();
                reads_ready = true;
            }
        }


        ComplexIterableStorage<Generator<std::vector<std::pair<ReadDirection, PathIterator>>::iterator, PassingRead>> getPassing(Vertex &v) const &;
        const std::vector<std::pair<spg::ReadDirection, PathIterator>> &getReadPositions(Edge &edge) const;
        const std::vector<Segment<Vertex>> &getInnerReads(Vertex &v) const;

        bool checkReadIndexConsistency() const;

        void fireAddVertex(Vertex &vertex) override;
        void fireAddEdge(Edge &edge) override;
        void fireDeleteVertex(Vertex &vertex) override;
        void fireDeleteEdge(Edge &edge) override;

        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;
        void fireMergePath(const GraphPath &path, Vertex &new_vertex) override;
        void fireMergeLoop(const GraphPath &path, Vertex &new_vertex) override;
    };
}