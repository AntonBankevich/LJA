#pragma once

#include "dbg_read_alignment_storage.hpp"
#include "graph_printing.hpp"
#include "assembly_graph/visualization.hpp"

namespace dbg {
    struct Subdataset {
        Subdataset(ag::Component component) : component(std::move(component)) {}

        ag::Component component;
        std::vector<ag::AlignedRead *> reads;
        std::string id = "";

        void Save(const std::experimental::filesystem::path &dir,
                  const ag::Printer &printer) const {
            recreate_dir(dir);
            std::experimental::filesystem::path reads_file = dir / "reads.fasta";
            std::experimental::filesystem::path graph_file = dir / "graph.gfa";
            std::experimental::filesystem::path dot_file = dir / "graph.dot";
            printer.printDot(dot_file, component);
            printer.printGFA(graph_file, component);
            //dbg::printGFA(graph_file, component, true);
            //printDot(dot_file, component, labeler);
            std::ofstream os;
            os.open(reads_file);
            for (ag::AlignedRead *read: reads) {
                os << ">" << read->getId() << "\n" << read->getPath().Seq() << "\n";
            }
            os.close();
        }
    };

    inline void
    FillSubdatasets(std::vector<Subdataset> &result, const std::vector<dbg::DBGAlignedReadStorage *> &storages,
                    bool add_out_edges = true) {
        std::unordered_map<dbg::Vertex *, std::vector<size_t>> cmap;
        for (size_t i = 0; i < result.size(); i++) {
            for (dbg::Vertex &vert: result[i].component.vertices()) {
                cmap[&vert].emplace_back(i);
            }
        }
        for (dbg::DBGAlignedReadStorage *recordStorage: storages)
            for (ag::AlignedRead &read: recordStorage->getReads()) {
                if (!read.valid())
                    continue;
                ag::GraphPath al = read.getPath();
                std::vector<size_t> cids;
                for (Vertex & vertex : al.innerVertices()) {
                    if (cmap.find(&vertex) != cmap.end())
                        cids.insert(cids.end(), cmap[&vertex].begin(), cmap[&vertex].end());
                }
                const std::vector<size_t> &other = cmap[&al.getStart()];
                if (al.isSingleton()) {
                    if (add_out_edges) {
                        for (Vertex &vertex : al.vertices())
                            cids.insert(cids.end(), cmap[&vertex].begin(), cmap[&vertex].end());
                    } else {
                        size_t p2 = 0;
                        for (size_t c1: other) {
                            while (p2 < other.size() && other[p2])
                                p2++;
                            if (p2 < other.size() && c1 == other[p2]) {
                                cids.push_back(c1);
                            }
                        }
                    }
                }
                std::sort(cids.begin(), cids.end());
                cids.erase(std::unique(cids.begin(), cids.end()), cids.end());
                for (size_t cid: cids) {
                    result[cid].reads.emplace_back(&read);
                }
            }
    }
}
