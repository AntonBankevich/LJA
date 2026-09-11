#pragma once

#include "unique_vertex_storage.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include "assembly_graph/visualization.hpp"
#include "assembly_graph/data_structures/component.hpp"
#include "common/dir_utils.hpp"
#include "common/string_utils.hpp"
#include <dbg/dbg_read_alignment_storage.hpp>
#include <experimental/filesystem>
#include <unordered_map>
#include <string>

namespace spg {

//    Reads sequences from paths, aligns them to spg and registers them with path_tracker
//    (an ag::DLLAlignmentStorage already attached to spg as a listener).
    void PrepareDLLPathTracker(logging::Logger &logger, size_t threads, dbg::SparseDBG &spg, size_t w,
                const io::Library &paths, ag::DLLAlignmentStorage &path_tracker);

//    Listens to the same graph-editing events ag::DLLAlignmentStorage already tracks and,
//    whenever such an event leaves a tracked path anchored at a vertex (storage.contigNames(vertex)
//    non-empty), dumps a dot-file snapshot of that vertex's neighbourhood into a directory
//    dedicated to that path, numbered in the order snapshots are taken for that path.
//    Must be constructed (and therefore registered as a listener) after `storage`, so that by
//    the time this class's own fireXxx handler runs for a given edit, storage's vertex_map
//    already reflects the post-edit state for that same edit.
    class PathTracker : public ag::ResolutionListener {
    private:
//        Keeps track of the files already written for one tracked path: its output directory
//        and the monotonic counter used to number the next snapshot. Printing itself stays in
//        PathTracker; this struct only ever hands out where the next file should go.
        struct PathFigures {
            std::experimental::filesystem::path dir;
            size_t cnt = 0;

            explicit PathFigures(std::experimental::filesystem::path dir) : dir(std::move(dir)) {
                ensure_dir_existance(this->dir);
            }

            std::experimental::filesystem::path nextFile(const std::string &event_tag) {
                std::experimental::filesystem::path fname = dir / (itos(cnt, 4) + "_" + event_tag + ".dot");
                cnt++;
                return fname;
            }
        };

        const ag::DLLAlignmentStorage *storage;
        ag::Printer printer;
        std::experimental::filesystem::path dir;
        size_t radius;
        size_t max_size;
        std::unordered_map<std::string, PathFigures> paths;

        PathFigures &figuresFor(const std::string &contig_name);
        void drawFiguresForVertex(Vertex &vertex, const std::string &event_tag);

    public:
        PathTracker(ag::ResolutionFire &fire, const ag::DLLAlignmentStorage &storage,
                    const ag::Printer &printer, std::experimental::filesystem::path dir,
                    size_t radius = 10000, size_t max_size = 100);

        void fireEdgeToSupreVertex(Vertex &v, Edge &e) override;
        void fireMergePath(const ag::RAGraphPath &path, Vertex &new_vertex) override;
        void fireMergePathToEdge(const ag::RAGraphPath &path, Edge &new_edge) override;
        void fireMergeTipsToEdge(Edge &new_edge, Edge &left, Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) override;
        void fireSplitEdge(Edge &edge, const ag::RAGraphPath &split) override;
        void fireResolveVertex(Vertex &core, const ag::VertexResolutionResult &resolution) override;
    };
}
