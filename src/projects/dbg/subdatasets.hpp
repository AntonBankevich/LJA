#pragma once

#include "graph_alignment_storage.hpp"
#include "visualization.hpp"
#include "graph_printing.hpp"

struct Subdataset {
    Subdataset(dbg::Component component) : component(std::move(component)) {}
    dbg::Component component;
    std::vector<AlignedRead *> reads;

    void Save(const std::experimental::filesystem::path &dir) const {
        recreate_dir(dir);
        std::experimental::filesystem::path reads_file = dir / "reads.fasta";
        std::experimental::filesystem::path graph_file = dir / "graph.fasta";
        std::experimental::filesystem::path dot_file = dir / "graph.dot";
        dbg::printFasta(graph_file, component);
        printDot(dot_file, component);
        std::ofstream os;
        os.open(reads_file);
        for(AlignedRead * read : reads) {
            os << ">" << read->id << "\n" << read->path.getAlignment().Seq() << "\n";
        }
        os.close();
    }
};

std::vector<Subdataset> SubdatasetSplit(dbg::SparseDBG &dbg, const std::vector<RecordStorage *> &storages, size_t min_len, bool add_out_edges = true) {
    std::vector<dbg::Component> components = dbg::LengthSplitter(min_len).splitGraph(dbg);
    std::vector<Subdataset> result;
    std::unordered_map<dbg::Vertex *, size_t> cmap;
    for(size_t i = 0; i < components.size(); i++) {
        result.emplace_back(components[i]);
        for(dbg::Vertex &vert : result.back().component.vertices()) {
            cmap[&vert] = i;
        }
    }
    for(RecordStorage *recordStorage : storages)
        for(AlignedRead &read : *recordStorage) {
            if(!read.valid())
                continue;
            dbg::GraphAlignment al = read.path.getAlignment();
            std::vector<size_t> cids;
            for(size_t i = 1; i < al.size(); i++) {
                cids.emplace_back(cmap[&al.getVertex(i)]);
            }
            if(al.size() == 1 && (add_out_edges || cmap[&al.getVertex(0)] == cmap[&al.getVertex(1)])) {
                cids.emplace_back(cmap[&al.getVertex(0)]);
                cids.emplace_back(cmap[&al.getVertex(1)]);
            }
            std::sort(cids.begin(), cids.end());
            cids.erase(std::unique(cids.begin(), cids.end()), cids.end());
            for(size_t cid : cids) {
                result[cid].reads.emplace_back(&read);
            }
        }
    return std::move(result);
}