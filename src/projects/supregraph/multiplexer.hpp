#pragma once

#include "vertex_resolution.hpp"
#include "supregraph.hpp"
#include "read_storage.hpp"
#include <unordered_set>
#include <common/logging.hpp>
#include <assembly_graph/data_structures/suffix_tracker.hpp>

namespace spg {

    class Multiplexer {
    private:
        SupreGraph &graph;
        DecisionRule &rule;
        size_t max_core_length;
        std::unordered_set<VertexId> core_queue;//Store only canonical vertices
//        TODO: remove reads parameter
        ag::AlignedReadStorage<SPGTraits> &reads;
    public:
        Multiplexer(SupreGraph &graph, ag::AlignedReadStorage<SPGTraits> &reads, DecisionRule &rule, size_t max_core_length);
        Multiplexer(Multiplexer &&) = delete;
        Multiplexer(Multiplexer &) = delete;

        VertexResolutionResult multiplex(logging::Logger &logger, size_t threads, Vertex &vertex);

        VertexResolutionResult multiplex(logging::Logger &logger, size_t threads);

        bool finished() const {return core_queue.empty();}

        // No outer edges permitted
        void fullMultiplex(logging::Logger &logger, size_t threads);
    };
}