#pragma once

// Shared helpers for building small, hand-constructed ag::AssemblyGraph instances in tests,
// without going through k-mer hashing/SparseDBG construction.
//
// Every vertex sequence is generated from a monotonically increasing counter encoded in the
// A/C alphabet only. Since a vertex's reverse complement can only contain T/G bases, no vertex
// sequence can ever equal another vertex's (or its own) reverse complement - so all vertices
// are guaranteed distinct from their own and each other's RC, and no self-RC (palindromic)
// vertices are ever created by accident.

#include "assembly_graph/assembly_graph.hpp"
#include "assembly_graph/ag_algorithms.hpp"
#include "assembly_graph/data_structures/read_alignment_storage.hpp"
#include "assembly_graph/data_structures/suffix_tracker.hpp"
#include "common/logging.hpp"

namespace ag_test {

    inline Sequence MakeSeq(size_t counter, size_t k) {
        std::string s(k, 'A');
        for (size_t i = 0; i < k && i < sizeof(size_t) * 8; i++) {
            if ((counter >> i) & 1u)
                s[i] = 'C';
        }
        return Sequence(s);
    }

    struct GraphBuilder {
        ag::AssemblyGraph graph;
        size_t k;
        size_t counter = 0;

        explicit GraphBuilder(size_t k) : k(k) {}

        ag::Vertex &newVertex() {
            return graph.addVertex(MakeSeq(counter++, k));
        }

        // Any two fresh vertices can always be connected this way regardless of any real overlap:
        // full_seq starts with a's own label and ends with b's own label by construction.
        ag::Edge &connect(ag::Vertex &a, ag::Vertex &b) {
            return graph.addEdge(a, b, a.getSeq() + b.getSeq());
        }

        // Builds a linear chain of `length` fresh vertices and `length - 1` edges connecting them.
        // The chain is bounded by junctions on both ends (first vertex has inDeg 0, last has outDeg 0),
        // so AllUnbranchingPaths/MergeAllToEdges will treat it as a single unbranching run.
        std::vector<ag::EdgeId> buildChain(size_t length) {
            VERIFY(length >= 2);
            std::vector<ag::EdgeId> res;
            ag::Vertex *prev = &newVertex();
            for (size_t i = 1; i < length; i++) {
                ag::Vertex &next = newVertex();
                res.emplace_back(connect(*prev, next).getId());
                prev = &next;
            }
            return res;
        }

        // Builds `branchCount` fully independent (disjoint) linear chains of `branchLength` edges
        // each. Every vertex has out-degree/in-degree at most 1 (AssemblyGraph caps out-degree at 4,
        // matching the 4-nucleotide alphabet, so this deliberately avoids any shared branch point),
        // giving `branchCount` unbranching paths that MergeAllToEdges can merge fully in parallel
        // with no shared vertex between them.
        std::vector<std::vector<ag::EdgeId>> buildStar(size_t branchCount, size_t branchLength) {
            std::vector<std::vector<ag::EdgeId>> branches;
            for (size_t b = 0; b < branchCount; b++)
                branches.emplace_back(buildChain(branchLength + 1));
            return branches;
        }
    };

    // Builds a GraphPath spanning edges[begin..end) (a contiguous sub-range of a chain built by
    // buildChain/buildStar), with no left/right nucleotide cuts.
    inline ag::GraphPath MakePath(const std::vector<ag::EdgeId> &edges, size_t begin, size_t end) {
        VERIFY(begin < end && end <= edges.size());
        ag::GraphPath path(*edges[begin]);
        for (size_t i = begin + 1; i < end; i++)
            path += *edges[i];
        return path;
    }

    inline logging::Logger &TestLogger() {
        static logging::Logger logger(false);
        return logger;
    }
}
