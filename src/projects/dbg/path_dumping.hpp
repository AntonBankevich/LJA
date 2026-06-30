#pragma once
#include "sparse_dbg.hpp"
#include "sequences/contigs.hpp"
#include "dbg_read_alignment_storage.hpp"
#include <unordered_map>
#include <utility>
#include <assembly_graph/data_structures/splitters.hpp>
#include "assembly_graph/data_structures/component.hpp"
#include "dbg_graph_aligner.hpp"
#include "assembly_graph/visualization.hpp"

namespace ag {
    //TODO rewrite for AssemblyGraph
    class AlignedContigStorage : ResolutionListener {
    private:
        std::unordered_map<ag::ConstEdgeId, std::vector<ag::AlignmentChain<Contig, Edge>>> edge_alignments;
        std::unordered_map<ag::ConstVertexId, std::vector<ag::AlignmentChain<Contig, Vertex>>> vertex_alignments;
        std::vector<Contig> contigs;
        omp_lock_t writelock = {};

        void lock() {omp_set_lock(&writelock);}
        void unlock() {omp_unset_lock(&writelock);}
    public:
        std::vector<Contig>::iterator begin() {return contigs.begin();}
        std::vector<Contig>::iterator end() {return contigs.end();}
        explicit AlignedContigStorage(AssemblyGraph & graph) : ResolutionListener(graph, "AlignedContigsStorage") {
        }

        AlignedContigStorage(const AlignedContigStorage &) = delete;

        AlignedContigStorage(AlignedContigStorage &&other)  noexcept = default;

        void addContig(Contig &&contig);
        void print(std::ostream &os);
        std::function<std::string(const dbg::Edge &edge)> pathInfo() const;
        std::function<std::string(const dbg::Edge &edge)> colorer(const std::string &color = "brown") const;
        ag::EdgeInfo edgeInfo() const {return ag::EdgeInfo::Tooltiper(pathInfo()) + ag::EdgeInfo::Colorer(colorer());}

        template<class T>
        std::unordered_map<typename T::const_pointer_type, std::vector<ag::AlignmentChain<Contig, T>>> GroupByContig(const std::vector<ag::AlignmentChain<Contig, T>> &rec_list);
        void Fill(logging::Logger &logger, size_t threads, dbg::KmerIndex &index);

        void fireAddVertex(Vertex &v) override {lock(); vertex_alignments[v.getId()] = {}; unlock();}
        void fireAddEdge(Edge &e) override {lock(); edge_alignments[e.getId()] = {}; unlock();}
        void fireDeleteVertex(Vertex &v) override {lock(); vertex_alignments.erase(v.getId()); unlock();}
        void fireDeleteEdge(Edge &e) override {lock(); edge_alignments.erase(e.getId()); unlock();}

        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;

        // virtual void fireMergePath(const RAGraphPath &path, Vertex &new_vertex) {
        //     lock();
        //     std::vector<AlignmentChain<Contig, Vertex>> &new_als = vertex_alignments.at(new_vertex.getId());
        //     size_t shift = 0;
        //     size_t last_size = 0;
        //     for (Edge &edge: path.edges()) {
        //         vertex_alignments[edge.getFinish().getId()].emplace_back(edge.getAlignment());
        //     }
        //     unlock();
        // }
        virtual void fireMergeLoop(const ag::GraphPath &path, Vertex &new_vertex) {VERIFY(false);}
        virtual void fireMergePathToEdge(const RAGraphPath &path, Edge &new_edge) {}
        virtual void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                         const AlignmentForm &left_al, const AlignmentForm &right_al) {}
        virtual void fireSplitEdge(Edge &edge, const RAGraphPath &split) {VERIFY(false);}

        virtual void fireResetEdgeCodes(logging::Logger &logger, size_t threads, AssemblyGraph &graph) {}

        virtual void fireResolveVertex(Vertex &core, const VertexResolutionResult &resolution) {};

    };

    template<class T>
    std::unordered_map<typename T::const_pointer_type, std::vector<AlignmentChain<Contig, T>>> AlignedContigStorage::
    GroupByContig(const std::vector<AlignmentChain<Contig, T>> &rec_list) {
        std::vector<std::pair<typename T::const_pointer_type, std::vector<ag::AlignmentChain<Contig, T>>>> res;
        std::vector<ag::AlignmentChain<Contig, T> > next;
        for(ag::AlignmentChain<Contig, T> rec : rec_list) {
            if(!next.empty() && (next[0].seg_to.contig() != rec.seg_to.contig())) {
                res.emplace_back(next[0].seg_to.contig().getId(), std::move(next));
                next.clear();
            }
            next.emplace_back(rec);
        }
        if(!next.empty()) {
            res.emplace_back(next[0].seg_to.contig().getId(), std::move(next));
        }
        return {res.begin(), res.end()};
    }
}

namespace dbg {
    void FillAlignmentsStorage(logging::Logger &logger, size_t threads, ag::AlignedContigStorage &storage,dbg::KmerIndex &index);
    void PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const std::string &stage,
                dbg::SparseDBG &dbg, dbg::DBGAlignedReadStorage &readStorage, const io::Library &paths_lib, const io::Library &references_lib, bool small);
}

// inline void printEdge(std::ostream &os, dbg::Edge &edge, const std::string &extra_label = "",
//                const std::string &color = "black") {
//     dbg:: Vertex &end = edge.getFinish();
//     os << "\"" << edge.getStart().getId() << "\" -> \"" << end.getId() <<
//        "\" [label=\"" << edge.getInnerId() << " " << edge.getCode() << " " << edge.truncSize() << "(" << edge.getCoverage() << ")\"";
//     if(!extra_label.empty()) {
//         os << " labeltooltip=\"" << extra_label << "\"";
// //        os << "\\n"<<extra_label;
//     }
//     os << " color=\"" + color + "\"]\n";
// }

// namespace std {
//     inline std::function<std::string(const dbg::Edge &)>
//     operator+(const std::function<std::string(const dbg::Edge &)> &l1,
//               const std::function<std::string(const dbg::Edge &)> &l2) {
//         return [l1, l2](const dbg::Edge &edge) ->std::string {
//             std::string s1 = l1(edge);
//             std::string s2 = l2(edge);
//             if(s1.empty())
//                 return s2;
//             return s1 + "\\n" + s2;
//         };
//     }
// }

/*
inline void printDot(std::ostream &os, const ag::Component &component, const std::function<std::string(dbg::Edge &)> &labeler,
              const std::function<std::string(dbg::Edge &)> &edge_colorer) {
    os << "digraph {\nnodesep = 0.5;\n";
    std::unordered_set<dbg::VertexId> extended;
    for(dbg::Edge &edge : component.edges()) {
        extended.emplace(edge.getFinish().getId());
        extended.emplace(edge.getStart().getId());
    }
    for(dbg::VertexId vid : extended) {
        dbg::Vertex &vert = *vid;
        std::string color = component.covers(vert) ? "white" : "yellow";
        os << vert.getId();
        os << " [style=filled fillcolor=\"" + color + "\"";
        if(vert.size() < 10)
            os << " label=" << vert.getSeq();
        os << "]\n";
    }
    for(dbg::Edge &edge : component.edges()) {
        printEdge(os, edge, labeler(edge), edge_colorer(edge));
    }
    os << "}\n";
}


inline void printDot(std::ostream &os, const ag::Component &component) {
    const std::function<std::string(dbg::Edge &)> labeler = [](dbg::Edge &) {return "";};
    const std::function<std::string(dbg::Edge &)> colorer = [](dbg::Edge &) {return "black";};
    printDot(os, component, labeler, colorer);
}

inline void printDot(std::ostream &os, const ag::Component &component, const std::function<std::string(dbg::Edge &)> &labeler) {
    const std::function<std::string(dbg::Edge &)> colorer = [](dbg::Edge &) {return "black";};
    printDot(os, component, labeler, colorer);
}

inline void printDot(const std::experimental::filesystem::path &f, const ag::Component &component, const std::function<std::string(dbg::Edge &)> &labeler,
                     const std::function<std::string(dbg::Edge &)> &edge_colorer) {
    std::ofstream os;
    os.open(f);
    printDot(os, component, labeler, edge_colorer);
    os.close();
}


inline void printDot(const std::experimental::filesystem::path &f, const ag::Component &component) {
    std::ofstream os;
    os.open(f);
    printDot(os, component);
    os.close();
}

inline void printDot(const std::experimental::filesystem::path &f, const ag::Component &component, const std::function<std::string(dbg::Edge &)> &labeler) {
    std::ofstream os;
    os.open(f);
    printDot(os, component, labeler);
    os.close();
}
*/
// inline void DrawSplit(const ag::Component &component, const std::experimental::filesystem::path &dir,
//                const std::function<std::string(const dbg::Edge &)> &labeler, const std::function<std::string(const dbg::Edge &)> &colorer,
//                size_t len = 100000) {
//     ag::Printer printer;
//     printer.setEdgeInfo(ag::ObjInfo<dbg::Edge>({labeler}, {colorer}, {}));
//     ensure_dir_existance(dir);
//     std::vector<ag::Component> split = ag::LengthSplitter(len).split(component);
//     for(size_t i = 0; i < split.size(); i++) {
//         std::experimental::filesystem::path f = dir / (std::to_string(i) + ".dot");
//         std::ofstream os;
//         os.open(f);
//         printer.printDot(os, split[i]);
//         //printDot(os, split[i], labeler, colorer);
//         os.close();
//     }
// }
//
// inline void DrawSplit(const ag::Component &component, const std::experimental::filesystem::path &dir,
//                       const std::function<std::string(const dbg::Edge &)> &labeler, size_t len = 100000) {
//     std::function<std::string(const dbg::Edge &)> colorer = [](const dbg::Edge &){return "black";};
//     DrawSplit(component, dir, labeler, colorer, len);
// }
//
// inline void DrawSplit(const ag::Component &component, const std::experimental::filesystem::path &dir, size_t len = 100000) {
//     std::function<std::string(const dbg::Edge &)> labeler = [](const dbg::Edge &){return "";};
//     std::function<std::string(const dbg::Edge &)> colorer = [](const dbg::Edge &){return "black";};
//     DrawSplit(component, dir, labeler, colorer, len);
// }

