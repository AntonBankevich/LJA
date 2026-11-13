#pragma once

#include "assembly_graph/graph_listeners.hpp"

namespace spg {
//    Base for listeners that only work on a supregraph. The DBG-only operations are prohibited: they abort and cannot
//    be overridden. Such listeners must be attached after the graph stops being processed as a DBG.
    class SPGListener : public ag::ResolutionListener {
    public:
        using ag::ResolutionListener::ResolutionListener;

        void fireMergePathToEdge(const ag::RAGraphPath &path, ag::Edge &new_edge) final {VERIFY(false);}
        void fireMergeTipsToEdge(ag::Edge &new_edge, ag::Edge &left, ag::Edge &right,
                                 const AlignmentForm &left_al, const AlignmentForm &right_al) final {VERIFY(false);}
        void fireSplitEdge(ag::Edge &edge, const ag::RAGraphPath &split) final {VERIFY(false);}
    };
}
