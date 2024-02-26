#pragma once

#include "component.hpp"
#include "sparse_dbg.hpp"
namespace dbg {
    inline std::string DefaultEdgeName(Edge &edge) {
        return edge.getInnerId().str();
    }

    inline std::string SaveEdgeName(Edge &edge) {
        return edge.getInnerId().str() + "_" + edge.rc().getInnerId().str();
    }

    inline void printFasta(std::ostream &out, const Component &component,
                           const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        for(Edge &edge : component.edgesUnique()) {
            out << ">" << name(edge) << "\n" << edge.getSeq() << "\n";
        }
//        size_t cnt = 0;
//        size_t masked_cnt = 1;
//        for (Edge &edge : component.edges()) {
//            Sequence edge_seq = edge.getStart().getSeq() + edge.truncSeq();
//            Vertex &end = edge.getFinish();
//            out << ">" << cnt << "_";
//            if(mask && !component.contains(edge.getStart())) {
//                out << masked_cnt << "0000";
//                masked_cnt++;
//            }
//            out << edge.getStart().getInnerId();
//            out << "_";
//            if(mask && !component.contains(edge.getFinish())) {
//                out << masked_cnt << "0000";
//                masked_cnt++;
//            }
//            out << end.getInnerId();
//            out << "_" << edge.truncSize() << "_" << edge.getCoverage() << "\n";
//            out << edge_seq << "\n";
//            cnt++;
//        }
    }

    inline void printAssembly(std::ostream &out, const Component &component,
                              const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        size_t cnt = 0;
        for (Edge &edge : component.edgesUnique()) {
            Sequence edge_seq = edge.getStart().getSeq() + edge.truncSeq();
            Vertex &end = edge.getFinish();
            out << ">" << cnt << "_" << name(edge) << "_" << edge.truncSize()
                << "_" << edge.getCoverage() << "\n";
            out << edge_seq << "\n";
            cnt++;
        }
    }

    inline Sequence cheatingCutStart(Sequence seq, unsigned char c, size_t min_size, size_t k) {
        size_t pos = seq.size() - min_size;
        while(pos > 0 && seq[pos + k] != c)
            pos--;
        return seq.Subseq(pos, seq.size());
    }

    inline void printFasta(const std::experimental::filesystem::path &outf, const Component &component,
                           const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        std::ofstream out;
        out.open(outf);
        printFasta(out, component, name);
        out.close();
    }

    inline void printAssembly(const std::experimental::filesystem::path &outf, const Component &component,
                              const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        std::ofstream out;
        out.open(outf);
        printAssembly(out, component, name);
        out.close();
    }

    inline void printFasta(const std::experimental::filesystem::path &outf, SparseDBG &dbg,
                           const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        std::ofstream out;
        out.open(outf);
        printFasta(out, Component(dbg), name);
        out.close();
    }

    inline void printAssembly(const std::experimental::filesystem::path &outf, SparseDBG &dbg,
                              const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        std::ofstream out;
        out.open(outf);
        printAssembly(out, Component(dbg), name);
        out.close();
    }

    inline void printGFA(std::ostream &out, const Component &component, bool calculate_coverage,
                         const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        std::unordered_map<const Edge *, std::string> eids;
        for (Edge &edge : component.edgesUnique()) {
            VERIFY(edge.isCanonical());
            eids[&edge] = name(edge);
            eids[&edge.rc()] = name(edge);
            if (calculate_coverage)
                out << "S\t" << name(edge) << "\t" << edge.getStart().getSeq() << edge.truncSeq()
                    << "\tKC:i:" << edge.intCov() << "\n";
            else
                out << "S\t" << name(edge) << "\t" << edge.getStart().getSeq() << edge.truncSeq() << "\n";
        }
        for (Vertex &vertex : component.verticesUnique()) {
            for (const Edge &out_edge : vertex) {
                std::string outid = eids[&out_edge];
                bool outsign = out_edge.isCanonical();
                for (const Edge &inc_edge : vertex.rc()) {
                    std::string incid = eids[&inc_edge];
                    bool incsign = inc_edge.isCanonical();
                    out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                        << (outsign ? "+" : "-") << "\t" << vertex.size() << "M" << "\n";
                }
            }
        }
    }

    inline void printGFA(const std::experimental::filesystem::path &outf, const Component &component, bool calculate_coverage,
                         const std::function<std::string(Edge &)> &name = &DefaultEdgeName) {
        std::ofstream out;
        out.open(outf);
        printGFA(out, component, calculate_coverage, name);
        out.close();
    }

    inline void printGraphAlignments(std::ostream &out, const std::vector<dbg::GraphPath> &als) {
        for(size_t i = 0; i < als.size(); i++) {
            out << ">" << i <<"\n" << als[i].Seq() << "\n";
        }
    }

    inline void printGraphAlignments(const std::experimental::filesystem::path &f, const std::vector<dbg::GraphPath> &als) {
        std::ofstream out;
        out.open(f);
        printGraphAlignments(out, als);
        out.close();
    }
}