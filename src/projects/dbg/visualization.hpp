#pragma once
#include "sparse_dbg.hpp"
#include "sequences/contigs.hpp"
#include <unordered_map>
#include <utility>
#include "component.hpp"
using namespace dbg;

class GraphAlignmentStorage {
private:
    std::unordered_map<const Edge *, std::vector<PerfectAlignment<Contig, Edge>>> alignments;
    std::vector<Contig*> stored_contigs;
    SparseDBG & dbg;

    void innerFill(const Contig &old_contig) {
        stored_contigs.emplace_back(new Contig(old_contig));
        Contig &contig = *stored_contigs.back();
        std::vector<PerfectAlignment<Contig, Edge>> path = GraphAligner(dbg).carefulAlign(contig);
        for(PerfectAlignment<Contig, Edge> &al : path) {
            alignments[&al.seg_to.contig()].emplace_back(al);
        }
    }

    void printEdge(std::ostream &os, Vertex & start, Edge &edge) {
        Vertex &end = *edge.end();
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 100000 << " -> ";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 100000 << "\n";
        std::vector<PerfectAlignment<Contig, Edge>> &als = alignments[&edge];
        for(auto & al : als) {
            os << "\n" << al.seg_from << "->" << al.seg_to;
        }
        os << "\n";
    }

public:
    explicit GraphAlignmentStorage(SparseDBG & dbg_) : dbg(dbg_) {
    }

    GraphAlignmentStorage(const GraphAlignmentStorage &) = delete;

    GraphAlignmentStorage(GraphAlignmentStorage &&other)  noexcept : dbg(other.dbg) {
        std::swap(other.alignments, alignments);
        std::swap(other.stored_contigs, stored_contigs);
    }

    ~GraphAlignmentStorage() {
        for(Contig * contig : stored_contigs) {
            delete contig;
        }
    }

    void fill(const Contig &contig) {
        innerFill(contig);
        innerFill(contig.RC());
    }

    void print(std::ostream &os) {
        for(auto &it : alignments) {
            const Edge &edge = *it.first;
            os << edge.getId() << "\n";
            if (alignments.find(&edge) == alignments.end())
                return;
            const std::vector<PerfectAlignment<Contig, Edge>> &als = alignments.find(&edge)->second;
            if (als.empty()) {
                return;
            }
            os << als[0].seg_from << "->" << als[0].seg_to;
            for (size_t i = 1; i < als.size(); i++) {
                const PerfectAlignment<Contig, Edge> &al = als[i];
                os << "\n" << al.seg_from << "->" << al.seg_to;
            }
        }
    }

    std::function<std::string(const Edge &edge)> labeler() const {
        std::function<std::string(const Edge &edge)> res = [this](const Edge &edge) {
            if (alignments.find(&edge) == alignments.end())
                return std::string("");
            std::stringstream ss;
            const std::vector<PerfectAlignment<Contig, Edge>> &als = alignments.find(&edge)->second;
            if (als.empty()) {
                return std::string("");
            }
            size_t num = std::min<size_t>(10, als.size());
            ss << als[0].seg_from << "->" << als[0].seg_to.coordinaresStr();
            for (size_t i = 1; i < num; i++) {
                const PerfectAlignment<Contig, Edge> &al = als[i];
                ss << "\\n" << al.seg_from << "->" << al.seg_to.coordinaresStr() << "\n";
            }
            return ss.str();
        };
        return res;
    }
};

inline void printEdge(std::ostream &os, Edge &edge, const std::string &extra_label = "",
               const std::string &color = "black") {
    Vertex &end = *edge.end();
    os << "\"" << edge.start()->getShortId() << "\" -> \"" << end.getShortId() <<
       "\" [label=\"" << "ACGT"[edge.seq[0]] << " " << edge.size() << "(" << edge.getCoverage() << ")";
    if(!extra_label.empty()) {
        os << "\\n"<<extra_label;
    }
    os << "\", color=\"" + color + "\"]\n";
}

namespace std {
    inline std::function<std::string(Edge &)> operator+(const std::function<std::string(Edge &)> &l1, const std::function<std::string(Edge &)> &l2) {
        return [l1, l2](Edge &edge) ->std::string {
            return l1(edge) + "\n" + l2(edge);
        };
    }
}

inline void printDot(std::ostream &os, const Component &component, const std::function<std::string(Edge &)> &labeler,
              const std::function<std::string(Edge &)> &edge_colorer) {
    os << "digraph {\nnodesep = 0.5;\n";
    std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> extended;
    for(Edge &edge : component.edgesUnique()) {
        extended.emplace(edge.end()->hash());
        extended.emplace(edge.start()->hash());
    }
    for(hashing::htype vid : extended) {
        for(Vertex * vit : component.graph().getVertices(vid)) {
            Vertex &vert = *vit;
            std::string color = component.covers(vert) ? "white" : "yellow";
            os << vert.getShortId() << " [style=filled fillcolor=\"" + color + "\"]\n";
        }
    }
    for(Edge &edge : component.edges()) {
        printEdge(os, edge, labeler(edge), edge_colorer(edge));
    }
    os << "}\n";
}


inline void printDot(std::ostream &os, const Component &component) {
    const std::function<std::string(Edge &)> labeler = [](Edge &) {return "";};
    const std::function<std::string(Edge &)> colorer = [](Edge &) {return "black";};
    printDot(os, component, labeler, colorer);
}

inline void printDot(std::ostream &os, const Component &component, const std::function<std::string(Edge &)> &labeler) {
    const std::function<std::string(Edge &)> colorer = [](Edge &) {return "black";};
    printDot(os, component, labeler, colorer);
}

inline void printDot(const std::experimental::filesystem::path &f, const Component &component, const std::function<std::string(Edge &)> &labeler,
                     const std::function<std::string(Edge &)> &edge_colorer) {
    std::ofstream os;
    os.open(f);
    printDot(os, component, labeler, edge_colorer);
    os.close();
}


inline void printDot(const std::experimental::filesystem::path &f, const Component &component) {
    std::ofstream os;
    os.open(f);
    printDot(os, component);
    os.close();
}

inline void printDot(const std::experimental::filesystem::path &f, const Component &component, const std::function<std::string(Edge &)> &labeler) {
    std::ofstream os;
    os.open(f);
    printDot(os, component, labeler);
    os.close();
}

inline void DrawSplit(const Component &component, const std::experimental::filesystem::path &dir,
               const std::function<std::string(Edge &)> &labeler, const std::function<std::string(Edge &)> &colorer,
               size_t len = 100000) {
    ensure_dir_existance(dir);
    std::vector<Component> split = LengthSplitter(len).split(component);
    for(size_t i = 0; i < split.size(); i++) {
        std::experimental::filesystem::path f = dir / (std::to_string(i) + ".dot");
        std::ofstream os;
        os.open(f);
        printDot(os, split[i], labeler, colorer);
        os.close();
    }
}

inline void DrawSplit(const Component &component, const std::experimental::filesystem::path &dir,
                      const std::function<std::string(Edge &)> &labeler, size_t len = 100000) {
    std::function<std::string(Edge &)> colorer = [](Edge &){return "black";};
    DrawSplit(component, dir, labeler, colorer, len);
}
inline void DrawSplit(const Component &component, const std::experimental::filesystem::path &dir, size_t len = 100000) {
    std::function<std::string(Edge &)> labeler = [](Edge &){return "";};
    std::function<std::string(Edge &)> colorer = [](Edge &){return "black";};
    DrawSplit(component, dir, labeler, colorer, len);
}