#pragma once

#include "assembly_graph/data_structures/component.hpp"
#include "sparse_dbg.hpp"

namespace dbg {
//TODO: move to Printer functionality and make it universal for all graphs
    inline void printFasta(std::ostream &out, const ag::Component &component,
                           const std::function<std::string(const ag::Edge &)> &name = &ag::DefaultEdgeName) {
        for(Edge &edge : component.edgesUnique()) {
            out << ">" << name(edge) << "\n" << edge.getSeq() << "\n";
        }
    }

    inline void printFasta(const std::experimental::filesystem::path &outf, const ag::Component &component,
                           const std::function<std::string(const ag::Edge &)> &name = &ag::DefaultEdgeName) {
        std::ofstream out;
        out.open(outf);
        printFasta(out, component, name);
        out.close();
    }

    inline void printFasta(const std::experimental::filesystem::path &outf, SparseDBG &dbg,
                           const std::function<std::string(const ag::Edge &)> &name = &ag::DefaultEdgeName) {
        std::ofstream out;
        out.open(outf);
        printFasta(out, ag::Component(dbg), name);
        out.close();
    }
}
