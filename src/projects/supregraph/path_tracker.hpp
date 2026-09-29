#pragma once

#include "unique_vertex_storage.hpp"
#include "assembly_graph/dll_path_storage.hpp"
#include <dbg/dbg_read_alignment_storage.hpp>

namespace spg {

//    Reads sequences from paths, aligns them to spg and registers them with path_tracker
//    (an ag::DLLAlignmentStorage already attached to spg as a listener). Drawing snapshots of
//    the tracked paths as the graph changes is ag::DLLAlignmentStorage's own job now (see
//    startLogging/startDrawing and ag::PathDrawer in dll_path_storage.hpp) -- there is no
//    separate listener class here anymore.
    void PrepareDLLPathTracker(logging::Logger &logger, size_t threads, dbg::SparseDBG &spg, size_t w,
                const io::Library &paths, ag::DLLAlignmentStorage &path_tracker);
}
