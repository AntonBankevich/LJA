#pragma once

#include "GraphContig.hpp"
#include <dbg/multi_graph.hpp>
#include <dbg/paths.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>

struct Detour {
    size_t start;
    size_t end;
    multigraph::GraphPath path;

    Detour(size_t start, size_t anEnd, multigraph::GraphPath path) : start(start), end(anEnd),
                                                           path(std::move(path)) {}
};


class BulgeFinder {
private:
    const multigraph::MultiGraph *mg;
    std::unordered_map<multigraph::ConstVertexId, std::unordered_map<multigraph::ConstVertexId, size_t>> min_dist;
    size_t max_size;
    size_t max_diff;
    static int INF;

    size_t getMinDist(const multigraph::MGVertex &v1, const multigraph::MGVertex &v2);
    bool recursiveFindBulges(std::vector<multigraph::GraphPath> &bulges, multigraph::GraphPath &bulge, const multigraph::MGEdge &last_edge, size_t clen, size_t tlen);
public:
    BulgeFinder(multigraph::MultiGraph &mg, size_t max_size, size_t max_diff) : mg(&mg), max_size(max_size), max_diff(max_diff) {
    }
    std::vector<Detour> findSimpleBulges(const multigraph::GraphPath &path);
    bool recursiveFindBulge(std::vector<Detour> &res, const multigraph::GraphPath &path, size_t start, size_t end, size_t path_len = INF);
    std::vector<Detour> findBulges(const multigraph::GraphPath &path);
};

std::unordered_map<std::string, std::vector<nano::GraphContig>> AlignOnt(logging::Logger &logger, const size_t threads,
             const std::experimental::filesystem::path &dir, const multigraph::MultiGraph &mg, const io::Library &ont_reads, bool reuse_alignment = false);

multigraph::GraphPath ContigToPath(const nano::GraphContig &al, multigraph::MultiGraph &graph);

multigraph::GraphPath FixPath(const nano::GraphContig &graphContig, BulgeFinder &bulgeFinder, multigraph::MultiGraph &graph);
