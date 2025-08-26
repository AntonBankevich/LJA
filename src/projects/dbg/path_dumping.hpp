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

namespace dbg {
    //TODO rewrite for AssemblyGraph
    class AlignedContigStorage {
    private:
        std::unordered_map<ag::ConstEdgeId, std::vector<ag::AlignmentChain<Contig, dbg::Edge>>> alignments;
        std::vector<Contig*> stored_contigs;
        const dbg::SparseDBG * dbg;

    public:
        explicit AlignedContigStorage(dbg::SparseDBG & dbg_) : dbg(&dbg_) {
        }

        AlignedContigStorage(const AlignedContigStorage &) = delete;

        AlignedContigStorage(AlignedContigStorage &&other)  noexcept = default;

        ~AlignedContigStorage() {
            for(Contig * contig : stored_contigs) {
                delete contig;
            }
        }

        void addContig(const Contig &contig);
        void Fill(size_t threads, dbg::KmerIndex &index);
        void print(std::ostream &os);
        std::function<std::string(const dbg::Edge &edge)> pathInfo() const;
        std::function<std::string(const dbg::Edge &edge)> colorer(const std::string &color = "brown") const;
        ag::EdgeInfo edgeInfo() const {return ag::EdgeInfo::Tooltiper(pathInfo()) + ag::EdgeInfo::Colorer(colorer());}
    };

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

