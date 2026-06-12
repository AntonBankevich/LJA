#pragma once

#include <common/logging.hpp>
#include "dbg/multi_graph.hpp"
#include "sequences/seqio.hpp"

class DecompressingManager {
private:
    ag::AssemblyGraph *g;
    std::experimental::filesystem::path dir;
    size_t min_overlap;
    size_t max_repeat;
    bool debug;
    std::unordered_map<ag::ConstVertexId, Segment<ag::Vertex>> segs;
    std::unordered_map<ag::VertexId, Sequence> uncompressed;
    std::unordered_map<ag::EdgeId, AlignmentForm> overlap_alignment;
public:
    DecompressingManager(ag::AssemblyGraph &assemblyGraph, std::experimental::filesystem::path dir,
                         size_t min_overlap, size_t max_repeat, bool debug) :
        g(&assemblyGraph), dir(dir), min_overlap(min_overlap), max_repeat(max_repeat), debug(debug) {}

    ag::AssemblyGraph &graph() const {return *g;}

    void ReduceAndUncompress(logging::Logger &logger, size_t threads, const io::Library &corrected_reads, const io::Library &reads);

    void printReduction(std::experimental::filesystem::path path) const;

    void calculateOverlaps(logging::Logger &logger, size_t threads);

    void printUncompressedGraph(logging::Logger &logger, size_t threads, std::experimental::filesystem::path path);
    std::vector<Contig> printAssembly(logging::Logger &logger, size_t threads);
};

// void printUncompressedResults(logging::Logger &logger, size_t threads, multigraph::MultiGraph &graph,
//                               const std::unordered_map<ag::VertexId, Segment<ag::Vertex>> &segs,
//                               const std::unordered_map<ag::VertexId , Sequence> &uncompressed_results,
//                               const std::experimental::filesystem::path &out_dir, bool debug);