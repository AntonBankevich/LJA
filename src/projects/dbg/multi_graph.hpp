#pragma once

#include <sequences/sequence.hpp>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <experimental/filesystem>
#include <fstream>
#include <common/string_utils.hpp>
#include <sequences/contigs.hpp>
#include <utility>
#include <common/iterator_utils.hpp>
#include <common/object_id.hpp>

namespace multigraph {

    class MultiGraph;
    class Vertex;
    class Edge;
    typedef Position<Edge> EdgePosition;
    typedef Segment<Edge> EdgeSegment;

    typedef ObjectId<Vertex, int> VertexId;
    typedef ObjectId<Edge, int> EdgeId;
    typedef ObjectId<const Vertex, int> ConstVertexId;
    typedef ObjectId<const Edge, int> ConstEdgeId;
    class Vertex {
        friend class MultiGraph;
    private:
        Sequence seq;
        int id;
        std::vector<Edge *> outgoing;
        std::string label;
        Vertex *_rc = nullptr;

        void setRC(Vertex &other) {_rc = &other;other._rc = this;}
        void addOutgoing(Edge &edge) {outgoing.emplace_back(&edge);}
        void removeOutgoing(Edge &edge) {outgoing.erase(std::remove(outgoing.begin(), outgoing.end(), &edge), outgoing.end());}
    public:
        explicit Vertex(Sequence seq, int id, std::string label = "") : seq(std::move(seq)), id(id), label(std::move(label)) {outgoing.reserve(4);}
        Vertex(): seq(""), id(0), label("") {VERIFY(false);}
        Vertex(const Vertex &) = delete;

        Vertex(Vertex && v) noexcept : seq(std::move(v.seq)), id(v.id), label(std::move(v.label)) {VERIFY (v.outgoing.empty() && v._rc == nullptr);}

        bool isCanonical() const {return seq <= !seq;}
        size_t inDeg() const {return _rc->outgoing.size();}
        size_t outDeg() const {return outgoing.size();}
        size_t size() const {return seq.size();};

        VertexId getId() {return {id, this};}
        int getInnerId() const {return id;}
        ConstVertexId getId() const {return {id, this};}
        const Sequence &getSeq() const {return seq;}
        const std::string &getLabel() const {return label;}
        Vertex &rc() {return *_rc;}
        const Vertex &rc() const {return *_rc;}
        TransformingIterator<std::vector<Edge *>::iterator, Edge> begin();
        TransformingIterator<std::vector<Edge *>::iterator, Edge> end();
        TransformingIterator<std::vector<Edge *>::const_iterator, const Edge> begin() const;
        TransformingIterator<std::vector<Edge *>::const_iterator, const Edge> end() const;

        Edge &operator[](size_t index) {return *outgoing[index];}
        const Edge &operator[](size_t index) const {return *outgoing[index];}
        Edge &front() {return *outgoing[0];}
        const Edge &front() const {return *outgoing[0];}
        Edge &back() {return *outgoing[outgoing.size() - 1];}
        const Edge &back() const {return *outgoing[outgoing.size() - 1];}

        bool operator==(const Vertex &other) const;
        bool operator!=(const Vertex &other) const {return !(*this == other);}
        bool operator<(const Vertex &other) const {return this->id < other.id;}
        bool operator>(const Vertex &other) const {return this->id < other.id;}
    };


    class EdgeData {
    private:
        Edge *edge;
        size_t cov = 0;
    public:
        explicit EdgeData(Edge &edge) : edge(&edge) {}
        void incCov(int delta) {cov += delta;}
        size_t intCov() const {return cov;}
        double getCoverage() const;
    };

    struct Edge {
        friend class MultiGraph;
    private:
        Sequence seq;
        int id{};
        size_t sz{};
        bool canonical{};
        std::string label;
        Vertex *_start = nullptr;
        Vertex *_end = nullptr;
        Edge *_rc = nullptr;
        mutable EdgeData data;
        void setRC(Edge &other) {_rc = &other;}
    public:
        typedef int id_type;
        explicit Edge(Vertex &start, Vertex &end, Sequence _seq, int id = 0, std::string label = "") :
                    seq(std::move(_seq)), id(id), label(std::move(label)),
                    _start(&start), _end(&end), data(*this) {
            sz = seq.size();
            canonical = seq <= !seq;
        }
        Edge(const Edge &) = delete;
        Edge(Edge && e) noexcept : seq(std::move(e.seq)), id(e.id), sz(e.sz), canonical(e.canonical), label(std::move(e.label)), data(*this) {
            VERIFY(e._start == nullptr && e._end == nullptr && e._rc == nullptr);
        }
        Edge() : data(*this) {VERIFY(false);}

        Sequence getSeq() const;
        Sequence truncSeq() const;
        Sequence fullSubseq(size_t from, size_t to) const {return seq.Subseq(from, to + getStartSize());}
        size_t truncSize() const {return sz - getStart().size();}
        size_t getStartSize() const {return getStart().size();}
        EdgeId getId() {return {id, this};}
        ConstEdgeId getId() const {return {id, this};}
        int getInnerId() const {return id;}
        const std::string &getLabel() const {return label;}
        size_t size() const {return sz;}
        Vertex &getStart() {return *_start;}
        Vertex &getFinish() {return *_end;}
        const Vertex &getStart() const {return *_start;}
        const Vertex &getFinish() const {return *_end;}
//        TODO: create reasonable coverage for multiplex graph
        EdgeData &getData() const {return data;}
        Edge &rc() {return *_rc;}
        const Edge &rc() const {return *_rc;}

        std::string getReverseLabel() const;

//TODO: get rid of the function below
        size_t overlap() const;
        bool isCanonical() const;
        bool isTip() const;
        bool isSimpleBridge();

//        We use convention that only edges of the same graph are going to be compared
        bool operator==(const Edge &other) const;
        bool operator!=(const Edge &other) const {return !(*this ==other);}
        bool operator<(const Edge &other) const {return this->id < other.id;}
        bool operator>(const Edge &other) const {return this->id < other.id;}
    };
    typedef std::unordered_map<std::string, std::vector<std::string>> deleted_edges_map;
    
    class MultiGraph {
    public:
        typedef multigraph::Vertex Vertex;
        typedef multigraph::Edge Edge;
    private:
        int maxVId = 0;
        int maxEId = 0;
        std::unordered_map<int, Vertex> vertices_map;
        std::unordered_map<int, Edge> edges_map;

    public:
        MultiGraph() = default;
        MultiGraph(MultiGraph &&other) = default;
        MultiGraph &operator=(MultiGraph &&other) = default;
        MultiGraph(const MultiGraph &) = delete;

        size_t edgeNumber() const {return edges_map.size();}
        size_t size() const {return vertices_map.size();}
        VertexId getVertexById(int id);
        ConstVertexId getVertexById(int id) const;
        EdgeId getEdgeById(int id);
        ConstEdgeId getEdgeById(int id) const;

        IterableStorage<TransformingIterator<std::unordered_map<int, Vertex>::iterator, Vertex>> vertices();
        IterableStorage<TransformingIterator<std::unordered_map<int, Vertex>::const_iterator, const Vertex>> vertices() const;
        IterableStorage<TransformingIterator<std::unordered_map<int, Edge>::iterator, Edge>> edges();
        IterableStorage<TransformingIterator<std::unordered_map<int, Edge>::const_iterator, const Edge>> edges() const;
        void printStats(std::ostream &os);
        void checkConsistency() const;

        Vertex &addVertex(const Sequence &seq, int id = 0, std::string label = "");
        Edge &addEdge(Vertex &from, Vertex &to, Sequence seq, int id = 0, std::string label = "");

        void internalRemoveEdge(Edge &edge);
        void internalRemoveIsolatedVertex (Vertex &vertex);
        deleted_edges_map attemptCompressVertex(Vertex &v);

//        deleted_edges_map deleteAndCompress(Edge &edge);

    };

    class MultiGraphHelper {
    public:
        MultiGraphHelper() = default;


        static MultiGraph LoadGFA(const std::experimental::filesystem::path &gfa_file, bool int_ids);
        static MultiGraph TransformToEdgeGraph(const MultiGraph &mg, size_t tip_size = 4001);
        static MultiGraph Delete(const MultiGraph &mg, const std::unordered_set<ConstEdgeId> &to_delete, const std::unordered_set<ConstVertexId> &to_delete_vertices = {});
        static MultiGraph MergeAllPaths(const MultiGraph &mg, bool verbose);

        static std::vector<EdgeId> uniquePathForward(Edge &edge);
        static std::vector<ConstEdgeId> uniquePathForward(const Edge &edge);
        static std::vector<EdgeId> uniquePath(Edge &edge);
        static std::vector<ConstEdgeId> uniquePath(const Edge &edge);

        static std::vector<Contig> extractContigs(const MultiGraph &mg, bool cut_overlaps);
        static void printExtractedContigs(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool cut_overlaps);
        static void printDot(const MultiGraph &mg, const std::experimental::filesystem::path &f);
//This is ugly duplication of code. It could be avoided using templates but it is ugly too. No viable solution for that in C++
        static void printEdgeGFA(const std::experimental::filesystem::path &f, const std::vector<ConstVertexId> &component, bool labels = false);
        static void printEdgeGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool labels = false);
        static void printVertexGFA(const std::experimental::filesystem::path &f, const std::vector<ConstVertexId> &component);
        static void printVertexGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f);
        static std::vector<std::vector<ConstVertexId>> split(const MultiGraph &mg);
    };

}
