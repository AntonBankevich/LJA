#pragma once
#include <utility>

#include "list_path.hpp"
#include "vertex_resolution.hpp"
#include "supregraph.hpp"
#include "unique_vertex_storage.hpp"
#include "assembly_graph/listeners.hpp"

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

    class LinkedSegments {
    private:
        size_t k;
        std::vector<Segment<Vertex>> segments;
        void normalize() {
            if(segments.empty())
                return;
            std::sort(segments.begin(), segments.end());
            size_t cur = 1;
            for(size_t i = 1; i < segments.size(); i++) {
                if(segments[i].left + k > segments[cur - 1].right) {
                    segments[cur] = segments[i];
                    cur++;
                } else {
                    segments[cur - 1].right = segments[i].right;
                }
            }
            segments.erase(segments.begin() + cur, segments.end());
        }
        LinkedSegments(size_t k, std::vector<Segment<Vertex>> segments) : k(k), segments(std::move(segments)) {
        }
    public:
        LinkedSegments(size_t k) : k(k) {
        }
        void operator+=(const Segment<Vertex> &seg) {
            segments.emplace_back(seg);
            normalize();
        }
        void operator+=(const LinkedSegments &other) {
            segments.insert(segments.end(), other.segments.begin(), other.segments.end());
            normalize();
        }
        LinkedSegments embed(const Segment<Vertex> &other) {
            std::vector<Segment<Vertex>> res;
            for(const auto &seg : segments) {
                res.emplace_back(seg.nest(other));
            }
            return {k, std::move(res)};
        }
        size_t longestExtension(size_t dive) const {
            size_t res = dive;
            for(const Segment<Vertex> &seg : segments) {
                if(seg.left + k <= dive)
                    res = std::max(res, seg.right);
                else
                    break;
            }
            return res;
        }
    };
    class PathIndex : public ag::ResolutionListener<SPGTraits> {
    private:
//        When reads_ready==false these two, as well as the path storage itself are invalidated
        mutable std::unordered_map<EdgeId, std::vector<std::pair<ReadDirection, PathIterator>>> read_index;
//        This is required only for chain rule and theoretical paper. Need to move it into a separate listener.
        mutable std::unordered_map<VertexId, size_t> vertex_dive;

        mutable std::unordered_map<VertexId, LinkedSegments> inner_index;
//        remapping is used only when reads_ready=false and allows to reconstruct read paths when necessary
        mutable Embedding embedding;
        PathStorage *storage;
        mutable bool reads_ready = false;

        void processPassing(Edge &edge, const spg::VertexResolutionResult &resolution);
        void prepareIndex() const;
    public:
        struct PassingRead {
            InOutEdgePair edges;
            ReadDirection direction;
            PathIterator position;

            PassingRead(const ReadDirection &direction, const PathIterator &position) :
                    edges(*position, *(position + 1)), direction(direction), position(position) {}
        };

//        Contract! No read from storage is fully contained within more than one core vertex.
//        This is always true if all reads have length at least K and spg was converted from DBG with K.
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
        void prepare() const {
            if(!reads_ready) {
                restorePaths();
                prepareIndex();
                embedding.clear();
                reads_ready = true;
            }
        }

//        void addInnerSegment(const Segment<Vertex> &seg) const {
//            Vertex &v = seg.contig();
////            if(v.outDeg() == 1 && v.front().isSuffix())
////                VERIFY(seg.cutRight() >= v.front().getFinish().size());
////            if(v.inDeg() == 1 && v.rc().front().isSuffix())
////                VERIFY(seg.cutLeft() >= v.rc().front().getFinish().size());
//            inner_index[v.getId()] += seg;
//        }


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