#pragma once

#include "assembly_graph/paths.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/compact_path.hpp"
#include <sequences/sequence.hpp>
#include <sequences/contigs.hpp>
#include <common/string_utils.hpp>
#include <common/iterator_utils.hpp>
#include <common/object_id.hpp>
#include <experimental/filesystem>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <fstream>
#include <utility>

namespace multigraph {

    class MultiGraph;
    class MGVertex;
    class MGEdge;
    typedef Position<MGEdge> EdgePosition;
    typedef Segment<MGEdge> EdgeSegment;

    class MGVertexData {
    protected:
        std::string label;
    public:
        MGVertexData(std::string label = "") : label(std::move(label)) {}
        MGVertexData RC() const {
            return *this;
        }
    };

    class MGEdgeData {
    protected:
        size_t cov = 0;
        std::string label;
    public:
        MGEdgeData(std::string label) : label(std::move(label)) {}
        void incCov(int delta) {
#pragma omp atomic
            cov += delta;
        }
        size_t intCov() const {return cov;}
        MGEdgeData RC() const {return {getReverseLabel()};}
        std::string getReverseLabel() const {
            if(label.empty())
                return {};
            std::vector<std::string> tokens = ::split(label, "_");
            std::string res;
            for (size_t i = tokens.size() -1; i > 0; i --)
                res += tokens[i] + "_";
            res += tokens[0];
            return res;
        }

        template<class I>
        static MGEdgeData Merge(I begin, I end) {
            std::vector<std::string> labels;
            size_t cov = 0;
            for(;begin != end; ++begin) {
                labels.emplace_back(begin->label);
                cov += begin->cov;
            }
            MGEdgeData res(join("_", labels));
            res.cov = cov;
            return std::move(res);
        }

    };


    struct MGTraits {
        typedef MGVertex Vertex;
        typedef MGEdge Edge;
        typedef MGEdgeData EdgeData;
        typedef MGVertexData VertexData;
    };

    class MGVertex : public ag::BaseVertex<MGTraits>, public MGVertexData {
    public:
        explicit MGVertex(id_type id, Sequence seq, VertexData data) : ag::BaseVertex<MGTraits>(id, std::move(seq)), MGVertexData(std::move(data)) {}
        explicit MGVertex(id_type id, bool canonical, VertexData data) : ag::BaseVertex<MGTraits>(id, canonical), MGVertexData(std::move(data)) {}
//        MGVertex(): seq(""), id(0), label("") {VERIFY(false);}
        MGVertex(const MGVertex &) = delete;

//        MGVertex(MGVertex && v) noexcept : seq(std::move(v.seq)), id(v.id), label(std::move(v.label)) {VERIFY (v.outgoing.empty() && v._rc == nullptr);}

        const std::string &getLabel() const {return label;}
    };


    class MGEdge : public ag::BaseEdge<MGTraits>, public MGEdgeData {
    public:
        explicit MGEdge(id_type id, MGVertex &start, MGVertex &end, Sequence _seq, MGEdgeData data) :
                ag::BaseEdge<MGTraits>(id, start, end, std::move(_seq)), MGEdgeData(std::move(data)) {}
        MGEdge(const MGEdge &) = delete;
//        MGEdge(MGEdge && e) noexcept : seq(std::move(e.seq)), id(e.id), sz(e.sz), canonical(e.canonical), label(std::move(e.label)), MGEdgeData(*this) {
//            VERIFY(e._start == nullptr && e._end == nullptr && e._rc == nullptr);
//        }
//        MGEdge() : MGEdgeData(*this) {VERIFY(false);}

        const std::string &getLabel() const {return label;}
        double getCoverage() const {return double(cov) / truncSize();}
//        TODO: create reasonable coverage for multiplex graph

//        bool isSimpleBridge();
    };
    typedef std::unordered_map<std::string, std::vector<std::string>> deleted_edges_map;

    typedef MGEdge Edge;
    typedef MGVertex Vertex;
    typedef MGEdge::EdgeId EdgeId;
    typedef MGVertex::VertexId VertexId;
    typedef MGEdge::ConstEdgeId ConstEdgeId;
    typedef MGVertex::ConstVertexId ConstVertexId;
    typedef ag::GraphPath<multigraph::MGTraits> GraphPath;
    typedef ag::CompactPath<multigraph::MGTraits> CompactPath;

    class MultiGraph : public ag::AssemblyGraph<MGTraits> {
    public:
        MultiGraph() = default;
        MultiGraph(MultiGraph &&other) = default;
        MultiGraph &operator=(MultiGraph &&other) = default;
        MultiGraph(const MultiGraph &) = delete;

//        Vertex &addVertex(const Sequence &seq, int id = 0, std::string label = "");
//        Edge &addEdge(Vertex &from, Vertex &to, Sequence seq, int id = 0, std::string label = "");

//        deleted_edges_map deleteAndCompress(Edge &edge);

    };

    class MultiGraphHelper {
    public:
        MultiGraphHelper() = default;

        static MultiGraph LoadGFA(const std::experimental::filesystem::path &gfa_file, bool int_ids);
        static MultiGraph LoadEdgeGFA(const std::experimental::filesystem::path &gfa_file, size_t K);
        static MultiGraph TransformToEdgeGraph(const MultiGraph &mg, size_t tip_size = 4001);
        static MultiGraph Delete(const MultiGraph &mg, const std::unordered_set<ConstEdgeId> &to_delete, const std::unordered_set<ConstVertexId> &to_delete_vertices = {});

        static std::vector<EdgeId> uniquePathForward(MGEdge &edge);
        static std::vector<ConstEdgeId> uniquePathForward(const MGEdge &edge);
        static std::vector<EdgeId> uniquePath(MGEdge &edge);
        static std::vector<ConstEdgeId> uniquePath(const MGEdge &edge);

        static std::vector<Contig> extractContigs(const MultiGraph &mg, bool cut_overlaps);
        static void printExtractedContigs(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool cut_overlaps);
        static void printDot(const MultiGraph &mg, const std::experimental::filesystem::path &f);
//This is ugly duplication of code. It could be avoided using templates but it is ugly too. No viable solution for that in C++
        static void printEdgeGFA(const std::experimental::filesystem::path &f, const std::vector<ConstVertexId> &component, bool labels = false);
        static void printEdgeGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool labels = false);
        static void printVertexGFA(const std::experimental::filesystem::path &f, const std::vector<ConstVertexId> &component);
        static void printVertexGFA(const MultiGraph &mg, const std::experimental::filesystem::path &f);
        static std::vector<std::vector<ConstVertexId>> split(const MultiGraph &mg);
        static void checkConsistency(multigraph::MultiGraph &mg) {
            std::unordered_set<ConstEdgeId> eset;
            std::unordered_set<ConstVertexId> vset;
            for(const MGEdge &edge: mg.edges()) {
                eset.emplace(edge.getId());
                VERIFY(edge.rc().getStart() == edge.getFinish().rc());
                VERIFY(edge.rc().rc() == edge);
            }
            for(const MGEdge &edge: mg.edges()) {
                VERIFY(eset.find(edge.rc().getId()) != eset.end());
            }
            for(const MGVertex &v : mg.vertices()) {
                vset.emplace(v.getId());
                VERIFY(v.rc().rc() == v);
                for(const MGEdge &edge : v) {
                    VERIFY(eset.find(edge.getId()) != eset.end());
                    VERIFY(edge.getStart() == v);
                }
            }
            for(const MGVertex &v : mg.vertices()) {
                VERIFY(vset.find(v.getId()) != vset.end());
            }
        }

    };

}
