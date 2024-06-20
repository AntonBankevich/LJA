#pragma once

#include "vertex_resolution.hpp"
#include "supregraph.hpp"
#include "read_storage.hpp"
#include <unordered_set>
#include <common/logging.hpp>

namespace spg {

    class Multiplexer {
    private:
        SupreGraph &graph;
        DecisionRule &rule;
        size_t max_core_length;
        std::unordered_set<VertexId> core_queue;//Store only canonical vertices
    public:
        Multiplexer(SupreGraph &graph, DecisionRule &rule, size_t max_core_length);
        Multiplexer(Multiplexer &&) = delete;
        Multiplexer(Multiplexer &) = delete;

        VertexResolutionResult multiplex(logging::Logger &logger, size_t threads, Vertex &vertex, PathIndex &index);

        VertexResolutionResult multiplex(logging::Logger &logger, size_t threads, PathIndex &index);

        bool finished() const {return core_queue.empty();}

        // No outer edges permitted
        void fullMultiplex(logging::Logger &logger, size_t threads, PathIndex &index);
    };
}