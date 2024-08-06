#pragma once

#include "graph_alignment_storage.hpp"
#include "visualization.hpp"
#include "graph_printing.hpp"
#include "assembly_graph/visualization.hpp"

namespace dbg {
    struct Subdataset {
        Subdataset(dbg::Component component) : component(std::move(component)) {}

        dbg::Component component;
        std::vector<ag::AlignedRead<DBGTraits> *> reads;
        std::string id = "";

        void Save(const std::experimental::filesystem::path &dir,
                  const std::function<std::string(const dbg::Edge &edge)> &label_func) const {
            recreate_dir(dir);
            std::experimental::filesystem::path reads_file = dir / "reads.fasta";
            std::experimental::filesystem::path graph_file = dir / "graph.gfa";
            std::experimental::filesystem::path dot_file = dir / "graph.dot";
            Printer<DBGTraits> printer;
            printer.addEdgeInfo(ObjInfo<Edge>({label_func},{}, {}));
            printer.printDot(dot_file, component);
            printer.printGFA(graph_file, component);
            //dbg::printGFA(graph_file, component, true);
            //printDot(dot_file, component, labeler);
            std::ofstream os;
            os.open(reads_file);
            for (ag::AlignedRead<DBGTraits> *read: reads) {
                os << ">" << read->id << "\n" << read->path.unpack().Seq() << "\n";
            }
            os.close();
        }
    };

    inline void
    FillSubdatasets(std::vector<Subdataset> &result, const std::vector<dbg::ReadAlignmentStorage *> &storages,
                    bool add_out_edges = true) {
        std::unordered_map<dbg::Vertex *, std::vector<size_t>> cmap;
        for (size_t i = 0; i < result.size(); i++) {
            for (dbg::Vertex &vert: result[i].component.vertices()) {
                cmap[&vert].emplace_back(i);
            }
        }
        for (dbg::ReadAlignmentStorage *recordStorage: storages)
            for (ag::AlignedRead<DBGTraits> &read: *recordStorage) {
                if (!read.valid())
                    continue;
                dbg::GraphPath al = read.path.unpack();
                std::vector<size_t> cids;
                for (size_t i = 1; i < al.size(); i++) {
                    if (cmap.find(&al.getVertex(i)) != cmap.end())
                        cids.insert(cids.end(), cmap[&al.getVertex(i)].begin(), cmap[&al.getVertex(i)].end());
                }
                const std::vector<size_t> &other = cmap[&al.getVertex(0)];
                if (al.size() == 1) {
                    if (add_out_edges) {
                        for (size_t i = 0; i < 2; i++)
                            cids.insert(cids.end(), cmap[&al.getVertex(i)].begin(), cmap[&al.getVertex(i)].end());
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
