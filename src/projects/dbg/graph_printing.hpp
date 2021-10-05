#pragma once

#include "component.hpp"
#include "sparse_dbg.hpp"
namespace dbg {
    inline void printFasta(std::ostream &out, const Component &component, bool mask = false) {
        size_t cnt = 0;
        size_t masked_cnt = 1;
        for (Edge &edge : component.edges()) {
            Sequence edge_seq = edge.start()->seq + edge.seq;
            Vertex &end = *edge.end();
            out << ">" << cnt << "_";
            if(mask && !component.contains(*edge.start())) {
                out << masked_cnt << "0000";
                masked_cnt++;
            }
            out << edge.start()->hash() << int(edge.start()->isCanonical());
            out << "_";
            if(mask && !component.contains(*edge.end())) {
                out << masked_cnt << "0000";
                masked_cnt++;
            }
            out << end.hash() << int(end.isCanonical());
            out << "_" << edge.size() << "_" << edge.getCoverage() << "\n";
            out << edge_seq << "\n";
            cnt++;
        }
    }

    inline void printAssembly(std::ostream &out, const Component &component) {
        size_t cnt = 0;
        for (Edge &edge : component.edgesUnique()) {
            Sequence edge_seq = edge.start()->seq + edge.seq;
            Vertex &end = *edge.end();
            out << ">" << cnt << "_" << edge.start()->hash() << int(edge.start()->isCanonical()) <<
                        "_" << end.hash() << int(end.isCanonical()) << "_" << edge.size()
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

    inline void cheatingFasta(const std::experimental::filesystem::path &outf, const Component &component, size_t cut) {
        std::ofstream out;
        out.open(outf);
        size_t cnt = 0;
        size_t k = component.graph().hasher().getK();
        for (Edge &edge : component.edges()) {
            Sequence edge_seq = edge.start()->seq + edge.seq;
            if(edge.size() > cut) {
                if(!component.contains(*edge.start())) {
                    edge_seq = cheatingCutStart(edge_seq, edge.seq[0], cut, k);
                } else if(!component.contains(*edge.end())) {
                    edge_seq = !cheatingCutStart(!edge_seq, edge.rc().seq[0], cut, k);
                }
            }
            Vertex &end = *edge.end();
            out << ">" << cnt << "_" << edge.start()->hash() << int(edge.start()->isCanonical()) <<
                "_" << end.hash() << int(end.isCanonical()) << "_" << edge_seq.size() - edge.start()->seq.size()
                << "_" << edge.getCoverage() << "\n";
            out << edge_seq << "\n";
            cnt++;
        }
        out.close();
    }

    inline void printFasta(const std::experimental::filesystem::path &outf, const Component &component, bool mask = false) {
        std::ofstream out;
        out.open(outf);
        printFasta(out, component, mask);
        out.close();
    }

    inline void printAssembly(const std::experimental::filesystem::path &outf, const Component &component) {
        std::ofstream out;
        out.open(outf);
        printAssembly(out, component);
        out.close();
    }

    inline void printGFA(std::ostream &out, const Component &component, bool calculate_coverage) {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        for (Edge &edge : component.edges()) {
            if (edge.start()->isCanonical(edge)) {
                if (calculate_coverage)
                    out << "S\t" << edge.oldId() << "\t" << edge.start()->seq << edge.seq
                        << "\tKC:i:" << edge.intCov() << "\n";
                else
                    out << "S\t" << edge.oldId() << "\t" << edge.start()->seq << edge.seq << "\n";
            }
        }
        for (Vertex &vertex : component.verticesUnique()) {
            for (const Edge &out_edge : vertex) {
                std::string outid = out_edge.oldId();
                bool outsign = vertex.isCanonical(out_edge);
                for (const Edge &inc_edge : vertex.rc()) {
                    std::string incid = inc_edge.oldId();
                    bool incsign = !vertex.rc().isCanonical(inc_edge);
                    out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                        << (outsign ? "+" : "-") << "\t" << component.graph().hasher().getK() << "M" << "\n";
                }
            }
        }
    }

    inline void printGFA(const std::experimental::filesystem::path &outf, const Component &component, bool calculate_coverage) {
        std::ofstream out;
        out.open(outf);
        printGFA(out, component, calculate_coverage);
        out.close();
    }
}