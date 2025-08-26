#pragma once

#include "GraphContig.hpp"
#include "ReadsAligner.h"
#include <dbg/multi_graph.hpp>
#include <assembly_graph/random_access_paths.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <dbg/aln_reads_reader.hpp>
#include <assembly_graph/visualization.hpp>

struct Detour {
    size_t start;
    size_t end;
    ag::PathPosition startPos;
    ag::PathPosition endPos;
    ag::GraphPath path;

    Detour(size_t start, size_t end, ag::PathPosition startPos,
           ag::PathPosition endPos, ag::GraphPath path) : start(start), end(end),
            startPos(startPos), endPos(endPos), path(std::move(path)) {}
};


class BulgeFinder {
private:
    const multigraph::MultiGraph *mg;
    std::unordered_map<multigraph::ConstVertexId, std::unordered_map<multigraph::ConstVertexId, size_t>> min_dist;
    size_t max_size;
    size_t max_diff;
    static int INF;

    size_t getMinDist(const multigraph::Vertex &v1, const multigraph::Vertex &v2);
    bool recursiveFindBulges(std::vector<ag::GraphPath> &bulges, ag::GraphPath &bulge, const multigraph::Edge &last_edge, size_t clen, size_t tlen);
public:
    BulgeFinder(multigraph::MultiGraph &mg, size_t max_size, size_t max_diff) : mg(&mg), max_size(max_size), max_diff(max_diff) {
    }
    std::vector<Detour> findSimpleBulges(const ag::GraphPath &path);
    bool recursiveFindBulge(std::vector<Detour> &res, const ag::GraphPath &path, size_t from, size_t to,
                            ag::PathPosition start, ag::PathPosition end, size_t path_len = INF);
    std::vector<Detour> findBulges(const ag::GraphPath &path);
};

std::unordered_map<std::string, std::vector<nano::GraphContig>>
AlignOnt(logging::Logger &logger, const size_t threads, const std::experimental::filesystem::path &dir,
         multigraph::MultiGraph &mg, const io::Library &ont_reads, bool reuse_alignment);
ag::GraphPath ContigToPath(const nano::GraphContig &al, multigraph::MultiGraph &graph);

ag::GraphPath FixPath(const nano::GraphContig &graphContig, BulgeFinder &bulgeFinder, multigraph::MultiGraph &graph);
