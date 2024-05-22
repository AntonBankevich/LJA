#pragma once
#include "listeners.hpp"
#include "vertex_resolution.hpp"
#include "supregraph_base.hpp"
#include "assembly_graph/assembly_graph.hpp"
#include "sequences/contigs.hpp"
#include "assembly_graph/paths.hpp"
#include "assembly_graph/compact_path.hpp"

namespace spg {

    class SupreGraph : public ag::AssemblyGraph<SPGTraits>, public ResolutionFire {
    public:
        SupreGraph() = default;
        SupreGraph(SupreGraph &&other) = default;
        SupreGraph &operator=(SupreGraph &&other) = default;
        SupreGraph(const SupreGraph &) = delete;

        Vertex &addSPGVertex(Sequence seq, bool cyclic, bool inf_left, bool inf_right, Vertex::id_type id = Vertex::id_type());
//        Depricated. Now graph should always be in normal form even though it creates 1-in-1-out vertices.
        Vertex &outerEdgeToVertex(Edge &edge);
        void IsolateAndMark(Vertex &v);

//        Graph should be in normal form. ResolutionPlan should connect incoming edges with outgoing (not reverse-complement)
        VertexResolutionResult resolveVertex(Vertex &core, const VertexResolutionPlan &resolution);
        Vertex &mergePath(const GraphPath &path);
    };

}