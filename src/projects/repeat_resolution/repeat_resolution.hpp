//
// Created by Andrey Bzikadze on 10/13/21.
//

#pragma once

#include "../dbg/sparse_dbg.hpp"
#include "../dbg/graph_alignment_storage.hpp"
#include "rr_graph.hpp"
#include <graphlite/serialize.hpp>
#include <vector>

namespace repeat_resolution {

    class RepeatResolver {
        dbg::SparseDBG &dbg;
        std::vector<RecordStorage *> storages;
        std::experimental::filesystem::path dir;
        bool debug;

        [[nodiscard]] std::vector<const Sequence *> get_edgeid2seq() const {
            std::vector<const Sequence *> edgeid2seq;
            for (const Edge &edge : dbg.edges()) {
                edgeid2seq.emplace_back(&edge.seq);
            }
            return edgeid2seq;
        }


    public:
        RepeatResolver(dbg::SparseDBG &dbg,
                       std::vector<RecordStorage *> storages,
                       std::experimental::filesystem::path dir,
                       bool debug) : dbg{dbg}, storages{std::move(storages)}, dir{std::move(dir)}, debug{debug} {
            std::experimental::filesystem::create_directory(this->dir);
        }

        void resolve_repeats(logging::Logger & logger) {
            std::vector<const Sequence*> edge2seq { get_edgeid2seq() };
            RRGraph rr_graph { dbg };
            rr_graph.serialize_to_dot(dir / "init_graph.dot");
            rr_graph.resolve_graph();
            rr_graph.serialize_to_dot(dir / "resolved_graph.dot");
        }
    };

} // End namespace repeat_resolution