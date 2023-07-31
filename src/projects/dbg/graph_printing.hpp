#pragma once

#include "component.hpp"
#include "sparse_dbg.hpp"
namespace dbg {
    inline void printFasta(std::ostream &out, const Component &component, bool mask = false) {
        size_t cnt = 0;
        size_t masked_cnt = 1;
        for (Edge &edge : component.edges()) {
            Sequence edge_seq = edge.getStart().getSeq() + edge.truncSeq();
            Vertex &end = edge.getFinish();
            out << ">" << cnt << "_";
            if(mask && !component.contains(edge.getStart())) {
                out << masked_cnt << "0000";
                masked_cnt++;
            }
            out << edge.getStart().hash() << int(edge.getStart().isCanonical());
            out << "_";
            if(mask && !component.contains(edge.getFinish())) {
                out << masked_cnt << "0000";
                masked_cnt++;
            }
            out << end.hash() << int(end.isCanonical());
            out << "_" << edge.truncSize() << "_" << edge.getData().getCoverage() << "\n";
            out << edge_seq << "\n";
            cnt++;
        }
    }

    inline void printAssembly(std::ostream &out, const Component &component) {
        size_t cnt = 0;
        for (Edge &edge : component.edgesUnique()) {
            Sequence edge_seq = edge.getStart().getSeq() + edge.truncSeq();
            Vertex &end = edge.getFinish();
            out << ">" << cnt << "_" << edge.getStart().hash() << int(edge.getStart().isCanonical()) <<
                "_" << end.hash() << int(end.isCanonical()) << "_" << edge.truncSize()
                << "_" << edge.getData().getCoverage() << "\n";
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
            Sequence edge_seq = edge.getStart().getSeq() + edge.truncSeq();
            if(edge.truncSize() > cut) {
                if(!component.contains(edge.getStart())) {
                    edge_seq = cheatingCutStart(edge_seq, edge.truncSeq()[0], cut, k);
                } else if(!component.contains(edge.getFinish())) {
                    edge_seq = !cheatingCutStart(!edge_seq, edge.rc().truncSeq()[0], cut, k);
                }
            }
            Vertex &end = edge.getFinish();
            out << ">" << cnt << "_" << edge.getStart().hash() << int(edge.getStart().isCanonical()) <<
                "_" << end.hash() << int(end.isCanonical()) << "_" << edge_seq.size() - edge.getStart().getSeq().size()
                << "_" << edge.getData().getCoverage() << "\n";
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
        std::unordered_map<const Edge *, std::string> eids;
        for (Edge &edge : component.edges()) {
            if (edge.isCanonical()) {
                eids[&edge] = edge.oldId();
                eids[&edge.rc()] = edge.oldId();
                if (calculate_coverage)
                    out << "S\t" << edge.oldId() << "\t" << edge.getStart().getSeq() << edge.truncSeq()
                        << "\tKC:i:" << edge.getData().intCov() << "\n";
                else
                    out << "S\t" << edge.oldId() << "\t" << edge.getStart().getSeq() << edge.truncSeq() << "\n";
            }
        }
        for (Vertex &vertex : component.verticesUnique()) {
            for (const Edge &out_edge : vertex) {
                std::string outid = eids[&out_edge];
                bool outsign = out_edge.isCanonical();
                for (const Edge &inc_edge : vertex.rc()) {
                    std::string incid = eids[&inc_edge];
                    bool incsign = inc_edge.isCanonical();
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

    inline void printGraphAlignments(std::ostream &out, const std::vector<DBGGraphPath> &als) {
        for(size_t i = 0; i < als.size(); i++) {
            out << ">" << i <<"\n" << als[i].Seq() << "\n";
        }
    }

    inline void printGraphAlignments(const std::experimental::filesystem::path &f, const std::vector<DBGGraphPath> &als) {
        std::ofstream out;
        out.open(f);
        printGraphAlignments(out, als);
        out.close();
    }
}