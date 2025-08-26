#pragma once

#include "assembly_graph/random_access_paths.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include <sequences/sequence.hpp>
#include <sequences/contigs.hpp>
#include <common/string_utils.hpp>
#include <common/iterator_utils.hpp>
#include <common/object_id.hpp>
#include <common/logging.hpp>
#include <experimental/filesystem>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <fstream>
#include <utility>
#include <assembly_graph/data_structures/component.hpp>

namespace multigraph {
    using ag::Vertex;
    using ag::Edge;
    using ag::VertexId;
    using ag::EdgeId;
    using ag::ConstVertexId;
    using ag::ConstEdgeId;

    typedef ag::AssemblyGraph MultiGraph;
    typedef Position<Edge> EdgePosition;
    typedef Segment<Edge> EdgeSegment;

    class LabelStorage : public ag::ResolutionListener {
    private:
        std::unordered_map<ConstEdgeId, std::vector<EdgeId>> labels;
    public:
        explicit LabelStorage(ag::AssemblyGraph &fire);

        void fireAddEdge(Edge &e) override {labels[e.getId()] = {e.getId()};}
        void fireDeleteEdge(Edge &e) override {labels.erase(e.getId());}

        void fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) override {VERIFY(false);}
        void fireMergeLoop(const ag::GraphPath  &path, Vertex &new_vertex) override {VERIFY(false);}
        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) override {
            std::vector<EdgeId> res;
            for(Edge &e: path.edges()) {
                std::vector<EdgeId> &tmp = labels.at(e.getId());
                res.insert(res.end(), tmp.begin(), tmp.end());
            }
            labels[new_edge.getId()] = std::move(res);
        }
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                             const AlignmentForm &left_al, const AlignmentForm &right_al) override {VERIFY(false);}
        void fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) override {VERIFY(false);}
        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, ag::AssemblyGraph &graph) override {VERIFY(false);}
        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override {VERIFY(false);};
        const std::vector<EdgeId> &getLabel(const Edge &edge) const {return labels.at(edge.getId());}
        std::string stringLabel(const Edge &edge) const {
            const std::vector<EdgeId> &label = labels.at(edge.getId());
            if(label.empty())
                return "";
            std::stringstream ss;
            ss << label.front();
            for(size_t i = 1; i < label.size(); i++) {
                ss << "_" << label[i];
            }
            return ss.str();
        }
        std::function<std::string(const Edge&)> getLabeler() const {
            return [this](const Edge &edge)->std::string{return stringLabel(edge);};
        }
    };




    class MultiGraphHelper {
    public:
        MultiGraphHelper() = default;

        static MultiGraph LoadGFA(const std::experimental::filesystem::path &gfa_file, bool int_ids);
        static MultiGraph LoadEdgeGFA(const std::experimental::filesystem::path &gfa_file, size_t K);
        static MultiGraph TransformToEdgeGraph(logging::Logger &logger, const MultiGraph &mg, size_t tip_size = 4001);
        static MultiGraph Delete(const MultiGraph &mg, const std::unordered_set<ConstEdgeId> &to_delete, const std::unordered_set<ConstVertexId> &to_delete_vertices = {});

//        static std::vector<EdgeId> uniquePathForward(Edge &edge);
//        static std::vector<ConstEdgeId> uniquePathForward(const Edge &edge);
//        static std::vector<EdgeId> uniquePath(Edge &edge);
//        static std::vector<ConstEdgeId> uniquePath(const Edge &edge);

        static std::vector<Contig> extractContigs(const MultiGraph &mg, bool cut_overlaps);
        static void printExtractedContigs(const MultiGraph &mg, const std::experimental::filesystem::path &f, bool cut_overlaps);
        static void printDot(const MultiGraph &mg, const std::experimental::filesystem::path &f);
        static void printDot2(const MultiGraph &mg, const std::experimental::filesystem::path &f);
//This is ugly duplication of code. It could be avoided using templates but it is ugly too. No viable solution for that in C++
    };

}
