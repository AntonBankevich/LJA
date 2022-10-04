#pragma once

#include "graph_alignment_storage.hpp"
#include "visualization.hpp"
#include "graph_printing.hpp"

struct Subdataset {
    Subdataset(dbg::Component component) : component(std::move(component)) {}
    dbg::Component component;
    std::vector<AlignedRead *> reads;
    std::string id = "";

    void Save(const std::experimental::filesystem::path &dir, const std::function<std::string(dbg::Edge &edge)>&labeler) const {
        recreate_dir(dir);
        std::experimental::filesystem::path reads_file = dir / "reads.fasta";
        std::experimental::filesystem::path graph_file = dir / "graph.fasta";
        std::experimental::filesystem::path dot_file = dir / "graph.dot";
        dbg::printFasta(graph_file, component);
        printDot(dot_file, component, labeler);
        std::ofstream os;
        os.open(reads_file);
        for(AlignedRead * read : reads) {
            os << ">" << read->id << "\n" << read->path.getAlignment().Seq() << "\n";
        }
        os.close();
    }
};

inline void FillSubdatasets(std::vector<Subdataset> &result, const std::vector<RecordStorage *> &storages, bool add_out_edges = true) {
    std::unordered_map<dbg::Vertex *, std::vector<size_t>> cmap;
    for(size_t i = 0; i < result.size(); i++) {
        for(dbg::Vertex &vert : result[i].component.vertices()) {
            cmap[&vert].emplace_back(i);
        }
    }
    for(RecordStorage *recordStorage : storages)
        for(AlignedRead &read : *recordStorage) {
            if(!read.valid())
                continue;
            dbg::GraphAlignment al = read.path.getAlignment();
            std::vector<size_t> cids;
            for(size_t i = 1; i < al.size(); i++) {
                if(cmap.find(&al.getVertex(i)) != cmap.end())
                    cids.insert(cids.end(), cmap[&al.getVertex(i)].begin(), cmap[&al.getVertex(i)].end());
            }
            const std::vector<size_t> &other = cmap[&al.getVertex(0)];
            if(al.size() == 1) {
                if(add_out_edges) {
                    for (size_t i = 0; i < 2; i++)
                        cids.insert(cids.end(), cmap[&al.getVertex(i)].begin(), cmap[&al.getVertex(i)].end());
                } else {
                    size_t p2 = 0;
                    for(size_t c1 : other) {
                        while(p2 < other.size() && other[p2])
                            p2++;
                        if(p2 < other.size() && c1 == other[p2]) {
                            cids.push_back(c1);
                        }
                    }
                }
            }
            std::sort(cids.begin(), cids.end());
            cids.erase(std::unique(cids.begin(), cids.end()), cids.end());
            for(size_t cid : cids) {
                result[cid].reads.emplace_back(&read);
            }
        }
}
