#pragma once

#include <dbg/paths.hpp>
#include "GraphContig.hpp"

struct Detour {
    size_t start;
    size_t end;
    MGGraphPath path;

    Detour(size_t start, size_t anEnd, MGGraphPath path) : start(start), end(anEnd),
                                                           path(std::move(path)) {}
};


class BulgeFinder {
private:
    const multigraph::MultiGraph *mg;
    std::unordered_map<multigraph::ConstVertexId, std::unordered_map<multigraph::ConstVertexId, size_t>> min_dist;
    size_t max_size;
    size_t max_diff;
    static int INF;

    size_t getMinDist(const multigraph::Vertex &v1, const multigraph::Vertex &v2);
    bool recursiveFindBulges(std::vector<MGGraphPath> &bulges, MGGraphPath &bulge, const multigraph::Edge &last_edge, size_t clen, size_t tlen);
public:
    BulgeFinder(multigraph::MultiGraph &mg, size_t max_size, size_t max_diff) : mg(&mg), max_size(max_size), max_diff(max_diff) {
    }
    std::vector<Detour> findSimpleBulges(const MGGraphPath &path);
    bool recursiveFindBulge(std::vector<Detour> &res, const MGGraphPath &path, size_t start, size_t end, size_t path_len = INF);
    std::vector<Detour> findBulges(const MGGraphPath &path);
};

std::unordered_map<std::string, std::vector<nano::GraphContig>> AlignOnt(logging::Logger &logger, const size_t threads,
             const std::experimental::filesystem::path &dir, const multigraph::MultiGraph &mg, const io::Library &ont_reads, bool reuse_alignment = false);

MGGraphPath ContigToPath(const nano::GraphContig &al, multigraph::MultiGraph &graph);

MGGraphPath FixPath(const nano::GraphContig &graphContig, BulgeFinder &bulgeFinder, multigraph::MultiGraph &graph);
