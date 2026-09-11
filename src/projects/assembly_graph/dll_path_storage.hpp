#pragma once
#include "alignment_chain.hpp"
#include "dll_path.hpp"

namespace ag {
    class DLLAlignmentPath : public ds::DoublyLinkedList<AlignmentFragment>{
    public:
        DLLAlignmentPath() = default;
        DLLAlignmentPath(const std::vector<AlignmentChain<Contig, Edge>> &als);
        DLLAlignmentPath(Contig &contig, const GraphPath &path);
    };

    class DLLAlignmentStorage : public ResolutionListener {
    public:
        typedef DLLAlignmentPath::iterator DLLPathPosition;
    private:
        std::list<Contig> contigs;
        std::unordered_map<std::string, DLLAlignmentPath> alignments;
        std::unordered_map<ConstVertexId, std::vector<DLLPathPosition>> vertex_map;
        std::unordered_map<ConstEdgeId, std::vector<DLLPathPosition>> edge_map;

        void addContig(Contig & contig, DLLAlignmentPath && p);

        DLLPathPosition insertBefore(Edge &edge, DLLPathPosition pos, AlignmentFragment fragment);
        DLLPathPosition insertAfter(Edge &edge, DLLPathPosition pos, AlignmentFragment fragment);
        DLLPathPosition insertBefore(Vertex &vertex, DLLPathPosition pos, AlignmentFragment fragment);
        DLLPathPosition insertAfter(Vertex &vertex, DLLPathPosition pos, AlignmentFragment fragment);

        bool hasRecords(const Vertex &vertex) const {return vertex_map.find(vertex.getId()) != vertex_map.end();}
        bool hasRecords(const Edge &edge) const {return edge_map.find(edge.getId()) != edge_map.end();}

    public:
        DLLAlignmentStorage(ResolutionFire &fire) : ResolutionListener(fire, "DLLAlignmentStorage"){}

        //TODO: change interface to avoid direct usage of dbg related AlignmentChain objects
        DLLAlignmentPath &addContig(Contig new_contig, std::vector<AlignmentChain<Contig, Edge>> &als);
        void fireDeleteVertex(Vertex &v) override;
        void fireDeleteEdge(Edge &e) override;
        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) override;
        void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {VERIFY(false);}
        void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireSplitEdge(Edge &edge, const RAGraphPath &split) override;
        void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) override {}
        void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) override;

//        Names (contig ids, still RC-oriented e.g. "-name") of all contigs with a fragment
//        anchored at vertex as an inner fragment, deduplicated. Empty if vertex is untracked.
        std::vector<std::string> passingContigs(const Vertex &vertex) const;
        std::vector<std::string> passingForwardContigs(const Vertex &vertex) const;

        std::function<std::string(const Vertex &)> getVertexTooltipper() const;
        std::function<std::string(const Edge &)> getEdgeTooltipper() const;
        std::function<std::string(const Edge &)> getEdgeColorer() const;
        std::function<std::string(const Vertex &)> getVertexColorer() const;
        VertexInfo getVertexInfo() const;
        EdgeInfo getEdgeInfo() const;
        Printer getPrinter() const;

        void print(std::ostream &out);
    };
}
