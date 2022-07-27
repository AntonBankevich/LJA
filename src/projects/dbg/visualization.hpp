#pragma once
#include "sparse_dbg.hpp"
#include "sequences/contigs.hpp"
#include <unordered_map>
#include <utility>
#include "component.hpp"

class GraphAlignmentStorage {
private:
    std::unordered_map<const dbg::Edge *, std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>>> alignments;
    std::vector<Contig*> stored_contigs;
    dbg::SparseDBG & dbg;

    void innerFill(const Contig &old_contig) {
        stored_contigs.emplace_back(new Contig(old_contig));
        Contig &contig = *stored_contigs.back();
        std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> path = dbg::GraphAligner(dbg).carefulAlign(contig);
        for(dbg::PerfectAlignment<Contig, dbg::Edge> &al : path) {
            alignments[&al.seg_to.contig()].emplace_back(al);
        }
    }

    void printEdge(std::ostream &os, dbg::Vertex & start, dbg::Edge &edge) {
        dbg::Vertex &end = *edge.end();
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 100000 << " -> ";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 100000 << "\n";
        std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> &als = alignments[&edge];
        for(auto & al : als) {
            os << "\n" << al.seg_from << "->" << al.seg_to;
        }
        os << "\n";
    }

public:
    explicit GraphAlignmentStorage(dbg::SparseDBG & dbg_) : dbg(dbg_) {
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

    void addContig(const Contig &contig) {
        stored_contigs.emplace_back(new Contig(contig));
        stored_contigs.emplace_back(new Contig(contig.RC()));
    }

    void Fill(size_t threads) {
        ParallelRecordCollector<dbg::PerfectAlignment<Contig, dbg::Edge>> records(threads);
#pragma omp parallel for default(none) shared(stored_contigs, records)
        for(size_t i = 0; i < stored_contigs.size(); i++) {
            Contig &contig = *stored_contigs[i];
            std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> path = dbg::GraphAligner(dbg).carefulAlign(contig);
            for(dbg::PerfectAlignment<Contig, dbg::Edge> &al : path) {
                records.emplace_back(al);
            }
        }
        std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> rec_list = records.collect();
        __gnu_parallel::sort(rec_list.begin(), rec_list.end());
        std::vector<std::pair<const dbg::Edge *, std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>>>> res;
        std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> next;
        for(dbg::PerfectAlignment<Contig, dbg::Edge> rec : rec_list) {
            if(next.empty() || next[0].seg_to.contig() == rec.seg_to.contig()) {
                next.emplace_back(rec);
            } else {
                dbg::Edge *nedge = &next[0].seg_to.contig();
                res.emplace_back(nedge, std::move(next));
                next.clear();
            }
        }
        if(!next.empty()) {
            dbg::Edge *nedge = &next[0].seg_to.contig();
            res.emplace_back(nedge, std::move(next));
        }
        alignments = {res.begin(), res.end()};
    }

    void print(std::ostream &os) {
        for(auto &it : alignments) {
            const dbg::Edge &edge = *it.first;
            os << edge.getId() << "\n";
            if (alignments.find(&edge) == alignments.end())
                return;
            const std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> &als = alignments.find(&edge)->second;
            if (als.empty()) {
                return;
            }
            os << als[0].seg_from << "->" << als[0].seg_to;
            for (size_t i = 1; i < als.size(); i++) {
                const dbg::PerfectAlignment<Contig, dbg::Edge> &al = als[i];
                os << "\n" << al.seg_from << "->" << al.seg_to;
            }
        }
    }

    std::function<std::string(dbg::Edge &edge)> labeler() const {
        std::function<std::string(const dbg::Edge &edge)> res = [this](const dbg::Edge &edge) {
            if (alignments.find(&edge) == alignments.end())
                return std::string("");
            std::stringstream ss;
            const std::vector<dbg::PerfectAlignment<Contig, dbg::Edge>> &als = alignments.find(&edge)->second;
            if (als.empty()) {
                return std::string("");
            }
            size_t num = std::min<size_t>(10, als.size());
            ss << als[0].seg_from << "->" << als[0].seg_to.coordinaresStr();
            for (size_t i = 1; i < num; i++) {
                const dbg::PerfectAlignment<Contig, dbg::Edge> &al = als[i];
                ss << "\\n" << al.seg_from << "->" << al.seg_to.coordinaresStr();
            }
            return ss.str();
        };
        return res;
    }
};

inline void printEdge(std::ostream &os, dbg::Edge &edge, const std::string &extra_label = "",
               const std::string &color = "black") {
    dbg:: Vertex &end = *edge.end();
    os << "\"" << edge.start()->getShortId() << "\" -> \"" << end.getShortId() <<
       "\" [label=\"" << "ACGT"[edge.seq[0]] << " " << edge.size() << "(" << edge.getCoverage() << ")\"";
    if(!extra_label.empty()) {
        os << " labeltooltip=\"" << extra_label << "\"";
//        os << "\\n"<<extra_label;
    }
    os << " color=\"" + color + "\"]\n";
}

namespace std {
    inline std::function<std::string(dbg::Edge &)> operator+(const std::function<std::string(dbg::Edge &)> &l1, const std::function<std::string(dbg::Edge &)> &l2) {
        return [l1, l2](dbg::Edge &edge) ->std::string {
            std::string s1 = l1(edge);
            std::string s2 = l2(edge);
            if(s1.empty())
                return s2;
            return s1 + "\\n" + s2;
        };
    }
}

inline void printDot(std::ostream &os, const dbg::Component &component, const std::function<std::string(dbg::Edge &)> &labeler,
              const std::function<std::string(dbg::Edge &)> &edge_colorer) {
    os << "digraph {\nnodesep = 0.5;\n";
    std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> extended;
    for(dbg::Edge &edge : component.edgesUnique()) {
        extended.emplace(edge.end()->hash());
        extended.emplace(edge.start()->hash());
    }
    for(hashing::htype vid : extended) {
        for(dbg::Vertex * vit : component.graph().getVertices(vid)) {
            dbg::Vertex &vert = *vit;
            std::string color = component.covers(vert) ? "white" : "yellow";
            os << vert.getShortId() << " [style=filled fillcolor=\"" + color + "\"]\n";
        }
    }
    for(dbg::Edge &edge : component.edges()) {
        printEdge(os, edge, labeler(edge), edge_colorer(edge));
    }
    os << "}\n";
}


inline void printDot(std::ostream &os, const dbg::Component &component) {
    const std::function<std::string(dbg::Edge &)> labeler = [](dbg::Edge &) {return "";};
    const std::function<std::string(dbg::Edge &)> colorer = [](dbg::Edge &) {return "black";};
    printDot(os, component, labeler, colorer);
}

inline void printDot(std::ostream &os, const dbg::Component &component, const std::function<std::string(dbg::Edge &)> &labeler) {
    const std::function<std::string(dbg::Edge &)> colorer = [](dbg::Edge &) {return "black";};
    printDot(os, component, labeler, colorer);
}

inline void printDot(const std::experimental::filesystem::path &f, const dbg::Component &component, const std::function<std::string(dbg::Edge &)> &labeler,
                     const std::function<std::string(dbg::Edge &)> &edge_colorer) {
    std::ofstream os;
    os.open(f);
    printDot(os, component, labeler, edge_colorer);
    os.close();
}


inline void printDot(const std::experimental::filesystem::path &f, const dbg::Component &component) {
    std::ofstream os;
    os.open(f);
    printDot(os, component);
    os.close();
}

inline void printDot(const std::experimental::filesystem::path &f, const dbg::Component &component, const std::function<std::string(dbg::Edge &)> &labeler) {
    std::ofstream os;
    os.open(f);
    printDot(os, component, labeler);
    os.close();
}

inline void DrawSplit(const dbg::Component &component, const std::experimental::filesystem::path &dir,
               const std::function<std::string(dbg::Edge &)> &labeler, const std::function<std::string(dbg::Edge &)> &colorer,
               size_t len = 100000) {
    ensure_dir_existance(dir);
    std::vector<dbg::Component> split = dbg::LengthSplitter(len).split(component);
    for(size_t i = 0; i < split.size(); i++) {
        std::experimental::filesystem::path f = dir / (std::to_string(i) + ".dot");
        std::ofstream os;
        os.open(f);
        printDot(os, split[i], labeler, colorer);
        os.close();
    }
}

inline void DrawSplit(const dbg::Component &component, const std::experimental::filesystem::path &dir,
                      const std::function<std::string(dbg::Edge &)> &labeler, size_t len = 100000) {
    std::function<std::string(dbg::Edge &)> colorer = [](dbg::Edge &){return "black";};
    DrawSplit(component, dir, labeler, colorer, len);
}
inline void DrawSplit(const dbg::Component &component, const std::experimental::filesystem::path &dir, size_t len = 100000) {
    std::function<std::string(dbg::Edge &)> labeler = [](dbg::Edge &){return "";};
    std::function<std::string(dbg::Edge &)> colorer = [](dbg::Edge &){return "black";};
    DrawSplit(component, dir, labeler, colorer, len);
}