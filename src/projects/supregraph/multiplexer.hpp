#pragma once

#include "abstract_decision_rule.hpp"

#include "read_storage.hpp"
#include <unordered_set>
#include <common/logging.hpp>
#include <assembly_graph/data_structures/suffix_tracker.hpp>
#include <set>

namespace spg {

    class Multiplexer {
    private:
        ag::AssemblyGraph &graph;
        DecisionRule &rule;
        size_t max_core_length;
        std::set<std::pair<size_t, VertexId>> core_queue;//Store only canonical vertices
        std::deque<VertexId> merge_queue;
//        TODO: remove reads parameter
        ag::AlignedReadStorage &reads;
    public:
        Multiplexer(ag::AssemblyGraph &graph, ag::AlignedReadStorage &reads, DecisionRule &rule, size_t max_core_length);
        Multiplexer(Multiplexer &&) = delete;
        Multiplexer(Multiplexer &) = delete;

        void pushCore(Vertex &vertex);
        Vertex &popCore();

//        Lets callers drive multiplexing through a progressively loosened vertex-size threshold:
//        core_queue is already ordered by vertex size, so raising this and checking hasReadyCore()
//        is enough to gate which cores get resolved without touching multiplex()/merge().
        void setMaxCoreLength(size_t new_max_core_length) {max_core_length = new_max_core_length;}
        size_t getMaxCoreLength() const {return max_core_length;}
        bool hasReadyCore() const {return !core_queue.empty() && core_queue.begin()->first <= max_core_length;}
        bool hasPendingMerge() const {return !merge_queue.empty();}

        std::vector<VertexId> multiplex(logging::Logger &logger, size_t threads, Vertex &vertex);
        std::vector<VertexId> merge(logging::Logger &logger, size_t threads, Vertex &vertex);

        std::vector<VertexId> process(logging::Logger &logger, size_t threads);

        bool finished() const {return core_queue.empty() && merge_queue.empty();}

        // No outer edges permitted
        void fullMultiplex(logging::Logger &logger, size_t threads);
    };
}