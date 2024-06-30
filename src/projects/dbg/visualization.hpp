#pragma once
#include "sparse_dbg.hpp"
#include "sequences/contigs.hpp"
#include "graph_alignment_storage.hpp"
#include <unordered_map>
#include <utility>
#include <assembly_graph/splitters.hpp>
#include "assembly_graph/component.hpp"
#include "dbg_graph_aligner.hpp"

#ifdef USE_LIBTORCH
#include "torch/torch.h"
#endif // USE_LIBTORCH

//TODO rewrite for AssemblyGraph
class GraphPathStorage {
private:
    std::unordered_map<dbg::ConstEdgeId, std::vector<ag::AlignmentChain<Contig, dbg::Edge>>> alignments;
    std::vector<Contig*> stored_contigs;
    const dbg::SparseDBG * dbg;

public:
    explicit GraphPathStorage(dbg::SparseDBG & dbg_) : dbg(&dbg_) {
    }

    GraphPathStorage(const GraphPathStorage &) = delete;

    GraphPathStorage(GraphPathStorage &&other)  noexcept = default;

    ~GraphPathStorage() {
        for(Contig * contig : stored_contigs) {
            delete contig;
        }
    }

    void addContig(const Contig &contig) {
        stored_contigs.emplace_back(new Contig(contig));
        stored_contigs.emplace_back(new Contig(contig.RC()));
    }

    void Fill(size_t threads, dbg::KmerIndex &index) {
        ParallelRecordCollector<ag::AlignmentChain<Contig, dbg::Edge>> records(threads);
#pragma omp parallel for schedule(dynamic, 10) default(none) shared(stored_contigs, records, index)
        for(size_t i = 0; i < stored_contigs.size(); i++) {
            Contig &contig = *stored_contigs[i];
            std::vector<ag::AlignmentChain<Contig, dbg::Edge>> path = index.carefulAlign(contig);
            for(ag::AlignmentChain<Contig, dbg::Edge> &al : path) {
                records.emplace_back(al);
            }
        }
        std::vector<ag::AlignmentChain<Contig, dbg::Edge>> rec_list = records.collect();
        __gnu_parallel::sort(rec_list.begin(), rec_list.end());
        std::vector<std::pair<dbg::ConstEdgeId , std::vector<ag::AlignmentChain<Contig, dbg::Edge>>>> res;
        std::vector<ag::AlignmentChain<Contig, dbg::Edge>> next;
        for(ag::AlignmentChain<Contig, dbg::Edge> rec : rec_list) {
            if(!next.empty() && next[0].seg_to.contig() != rec.seg_to.contig()) {
                res.emplace_back(next[0].seg_to.contig().getId(), std::move(next));
                next.clear();
            }
            next.emplace_back(rec);
        }
        if(!next.empty()) {
            res.emplace_back(next[0].seg_to.contig().getId(), std::move(next));
        }
        alignments = {res.begin(), res.end()};
    }

    void print(std::ostream &os) {
        for(auto &it : alignments) {
            const dbg::Edge &edge = *it.first;
            os << edge.getInnerId() << "\n";
            if (alignments.find(edge.getId()) == alignments.end())
                return;
            const std::vector<ag::AlignmentChain<Contig, dbg::Edge>> &als = alignments.find(edge.getId())->second;
            if (als.empty()) {
                return;
            }
            os << als[0].seg_from << "->" << als[0].seg_to;
            for (size_t i = 1; i < als.size(); i++) {
                const ag::AlignmentChain<Contig, dbg::Edge> &al = als[i];
                os << "\n" << al.seg_from << "->" << al.seg_to;
            }
        }
    }

    std::function<std::string(dbg::Edge &edge)> labeler() const {
        std::function<std::string(const dbg::Edge &edge)> res = [this](const dbg::Edge &edge) {
            if (alignments.find(edge.getId()) == alignments.end())
                return std::string("");
            std::stringstream ss;
            const std::vector<ag::AlignmentChain<Contig, dbg::Edge>> &als = alignments.find(edge.getId())->second;
            if (als.empty()) {
                return std::string("");
            }
            size_t num = std::min<size_t>(10, als.size());
            ss << als[0].seg_from << "->" << als[0].seg_to.coordinaresStr();
            for (size_t i = 1; i < num; i++) {
                const ag::AlignmentChain<Contig, dbg::Edge> &al = als[i];
                ss << "\\n" << al.seg_from << "->" << al.seg_to.coordinaresStr();
            }
            return ss.str();
        };
        return res;
    }
};

inline void printEdge(std::ostream &os, dbg::Edge &edge, const std::string &extra_label = "",
               const std::string &color = "black") {
    dbg:: Vertex &end = edge.getFinish();
    os << "\"" << edge.getStart().getId() << "\" -> \"" << end.getId() <<
       "\" [label=\"" << edge.getInnerId() << " " << edge.nuclLabel() << " " << edge.truncSize() << "(" << edge.getCoverage() << ")\"";
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

#ifdef USE_LIBTORCH
inline void printPT(const std::experimental::filesystem::path &dir, const dbg::Component &component) {
    std::unordered_set<dbg::VertexId> extended;
    for(dbg::Edge &edge : component.edges()) {
        extended.emplace(edge.getFinish().getId());
        extended.emplace(edge.getStart().getId());
    }
    std::unordered_map<int64_t, int64_t> node_map;
    std::vector<int64_t> node_data;
    int new_id = 0;
    for(dbg::VertexId vid : extended) {
        dbg::Vertex &vert = *vid;
        // TODO: implement direct converstion to integer ids
        std::stringstream stream;
        stream << vert.getId();
        int node;
        stream >> node;
        node_map[node] = new_id++;
        node_data.push_back(node);
    }
    std::vector<int64_t> edge_index_data;
    std::vector<float> edge_attr_data;
    for(dbg::Edge &edge : component.edges()) {
        dbg:: Vertex &end = edge.getFinish();
        // TODO: implement direct converstion to integer ids
        std::stringstream ss_src, ss_dst;
        int64_t src, dst;
        ss_src << edge.getStart().getId();
        ss_src >> src;
        ss_dst << end.getId();
        ss_dst >> dst;
        edge_index_data.push_back(node_map[src]);
        edge_index_data.push_back(node_map[dst]);
        auto cov = edge.getCoverage();
        auto len = edge.truncSize();
        edge_attr_data.push_back(cov);
        edge_attr_data.push_back(len);
    }

    int64_t edge_index_len = edge_index_data.size() / 2;
    auto options_int = torch::TensorOptions().dtype(torch::kInt64);
    torch::Tensor edge_index = torch::from_blob(edge_index_data.data(), {edge_index_len, 2}, options_int).clone();
    // TODO: implement serialization as a function
    auto pickled_ei = torch::pickle_save(edge_index);
    std::ofstream ei_out(dir / "edge_index.pt", std::ios::out | std::ios::binary);
    ei_out.write(pickled_ei.data(), pickled_ei.size());
    ei_out.close();
    int64_t nodes_len = node_data.size();
    torch::Tensor nodes = torch::from_blob(node_data.data(), {1, nodes_len}, options_int).clone();
    auto pickled_n = torch::pickle_save(node_data);
    std::ofstream n_out(dir / "nodes.pt", std::ios::out | std::ios::binary);
    n_out.write(pickled_n.data(), pickled_n.size());
    int64_t edge_attr_len = edge_attr_data.size() / 2;
    auto options_float = torch::TensorOptions().dtype(torch::kFloat32);
    torch::Tensor edge_attrs = torch::from_blob(edge_attr_data.data(), {edge_attr_len, 2}, options_float).clone();
    auto pickled_ea = torch::pickle_save(edge_attrs);
    std::ofstream ea_out(dir / "edge_attrs.pt", std::ios::out | std::ios::binary);
    ea_out.write(pickled_ea.data(), pickled_ea.size());
    ea_out.close();
}
#endif // USE_LIBTORCH

inline void printDot(std::ostream &os, const dbg::Component &component, const std::function<std::string(dbg::Edge &)> &labeler,
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
    std::vector<dbg::Component> split = ag::LengthSplitter<dbg::DBGTraits>(len).split(component);
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
void PrintPaths(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, const std::string &stage,
               dbg::SparseDBG &dbg, dbg::ReadAlignmentStorage &readStorage, const io::Library &paths_lib, bool small);
