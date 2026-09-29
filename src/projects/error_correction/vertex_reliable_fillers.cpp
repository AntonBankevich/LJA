#include "vertex_reliable_fillers.hpp"

namespace ag {
    size_t VertexCoverageReliableFiller::normalizeReliability(AssemblyGraph &graph) const {
        size_t res = 0;
        std::vector<ag::EdgeId> eids = oneline::map(graph.edgesUnique().begin(), graph.edgesUnique().end(), IdTransformer<Edge>());
        for (EdgeId &eid : eids) {if (eid->isPrefix()) eid=eid->rc().getId();}
        std::sort(eids.begin(), eids.end(), [](EdgeId a, EdgeId b) { return a->fullSize() > b->fullSize(); });
        for (EdgeId &eid : eids) {
            if (eid->isSuffix() && eid->getStart().reliability == VertexReliability::reliable) {
                if (eid->getFinish().reliability != VertexReliability::reliable)
                    res++;
                eid->getFinish().reliability = eid->getFinish().rc().reliability =  VertexReliability::reliable;
            }
        }
        return res;
    }

    size_t VertexCoverageReliableFiller::fill(AssemblyGraph &graph) const {
        size_t cnt = 0;
        for (Vertex &v: graph.vertices()) {
            if (!v.hasCoverageInfo())
                v.reliability = VertexReliability::unknown;
            else if (v.getSPGCoverage() < threshold)
                v.reliability = VertexReliability::unreliable;
            else {
                v.reliability = VertexReliability::reliable;
                if (v.isCanonical())
                    ++cnt;
            }
        }
        cnt += normalizeReliability(graph);
        return cnt;
    }

    size_t VertexCoverageReliableFiller::refill(AssemblyGraph &graph) const {
        for (Vertex &v: graph.vertices())
            v.reliability = VertexReliability::unknown;
        return fill(graph);
    }

    size_t VertexCoverageReliableFiller::loggedRefill(logging::Logger &logger, AssemblyGraph &graph) const {
        logger.info() << "Running VertexCoverageReliableFiller" << std::endl;
        size_t res = refill(graph);
        logger.info() << "Marked " << res << " vertices as reliable" << std::endl;
        return res;
    }
}
